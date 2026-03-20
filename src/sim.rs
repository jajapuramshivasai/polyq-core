use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::f64::consts::PI;

pub use crate::qc::{read_qasm_file, write_qasm_file, write_qasm_string, QasmError};

/// Phase resolution parameter.
///
/// We store phases in units of π / 2^PHASE_BITS.
///
/// We make the remainder phase ring **2π-periodic** by using modulus 2^(PHASE_BITS+1).
/// This allows representing θ = π (i.e. -1) exactly, which is required for MCZ.
pub const PHASE_BITS: u32 = 16;

/// Phase integer in units of π / 2^PHASE_BITS, reduced modulo 2^(PHASE_BITS+1).
pub type Phase = u32;

#[inline(always)]
pub const fn phase_modulus_u32() -> u32 {
    1u32 << (PHASE_BITS + 1)
}

#[inline(always)]
pub const fn phase_mask_u32() -> u32 {
    (1u32 << (PHASE_BITS + 1)) - 1
}

#[inline(always)]
fn phase_add(a: Phase, b: Phase) -> Phase {
    (a + b) & (phase_mask_u32() as Phase)
}

#[inline(always)]
fn phase_sub(a: Phase, b: Phase) -> Phase {
    (a.wrapping_sub(b)) & (phase_mask_u32() as Phase)
}

#[inline(always)]
fn phase_from_u32(x: u32) -> Phase {
    (x & phase_mask_u32()) as Phase
}

#[inline(always)]
fn phase_to_angle_rad(p: Phase) -> f64 {
    // θ = π * p / 2^PHASE_BITS
    // With p reduced mod 2^(PHASE_BITS+1), this is periodic mod 2π.
    PI * (p as f64) / ((1u32 << PHASE_BITS) as f64)
}

/// Addition in Z4 ring (Clifford group arithmetic).
#[inline(always)]
fn z4_add(a: u8, b: u8) -> u8 {
    (a + b) & 3
}

/// Subtraction in Z4 ring.
#[inline(always)]
fn z4_sub(a: u8, b: u8) -> u8 {
    (a + 4 - b) & 3
}

/// Quantum gate representation for {Clifford} + {T,Diadic-RZ} circuits.
///
/// - H(q): Hadamard gate on qubit q
/// - Z(q): Pauli-Z gate (phase π)
/// - S(q): Phase gate (phase π/2)
/// - T(q): T gate (phase π/4)
/// - RZ(q, phase): Arbitrary Z rotation by phase units (mod 2π)
/// - CZ(a, b): Controlled-Z gate between qubits a and b
/// - MCZ(ctrls): Multi-controlled Z: multiply by -1 iff all controls are 1
#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),
    Z(usize),
    S(usize),
    T(usize),
    RZ(usize, Phase),
    CZ(usize, usize),
    MCZ(Vec<usize>),
}

/// Remainder term for non-Clifford / general phase tracking.
///
/// Each term is a weighted AND-monomial over Boolean variables.
#[derive(Clone, Debug)]
pub struct PhaseTerm {
    pub weight: Phase,
    pub vars: Vec<usize>,
}

/// Compiled phase polynomial for Clifford+T simulation.
///
/// Contains all information needed for fast amplitude and statevector computation.
#[derive(Clone, Debug)]
pub struct CompiledPhasePoly {
    pub num_qubits: usize,
    pub num_vars: usize,
    pub num_h: usize,
    pub output_vars: Vec<usize>,
    pub b4: Vec<FixedBitSet>,
    pub v4: Vec<u8>,
    pub eps4: u8,
    pub rem: Vec<PhaseTerm>,
}

/// In-memory representation of a quantum circuit.
pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self {
        Self { num_qubits, gates: Vec::new() }
    }

    pub fn h(&mut self, q: usize) {
        self.gates.push(Gate::H(q));
    }
    pub fn z(&mut self, q: usize) {
        self.gates.push(Gate::Z(q));
    }
    pub fn s(&mut self, q: usize) {
        self.gates.push(Gate::S(q));
    }
    pub fn t(&mut self, q: usize) {
        self.gates.push(Gate::T(q));
    }
    pub fn rz(&mut self, q: usize, phase: Phase) {
        self.gates.push(Gate::RZ(q, phase));
    }
    pub fn cz(&mut self, q1: usize, q2: usize) {
        self.gates.push(Gate::CZ(q1, q2));
    }
    pub fn mcz(&mut self, ctrls: &[usize]) {
        let mut v = ctrls.to_vec();
        v.sort_unstable();
        v.dedup();
        self.gates.push(Gate::MCZ(v));
    }

    pub fn compile(&self) -> CompiledPhasePoly {
        let optimized_gates = Transpiler::optimize(self.num_qubits, self.gates.clone());
        compile_clifford_t(self.num_qubits, &optimized_gates)
    }
}

// --- 2. MULTI-PASS TRANSPILER (PHASE TELEPORTATION) ---
//
// This pass coalesces Z/S/T/RZ phases on a wire until they hit an H or CZ/MCZ.

pub struct Transpiler;

impl Transpiler {
    pub fn optimize(num_qubits: usize, gates: Vec<Gate>) -> Vec<Gate> {
        let mut current_gates = gates;
        loop {
            let start_count = current_gates.len();
            current_gates = Self::phase_teleportation_pass(num_qubits, current_gates);
            if current_gates.len() >= start_count {
                break;
            }
        }
        current_gates
    }

    fn phase_teleportation_pass(num_qubits: usize, gates: Vec<Gate>) -> Vec<Gate> {
        let mut optimized = Vec::with_capacity(gates.len());
        let mut pending_rz: Vec<i64> = vec![0; num_qubits]; // units π/2^PHASE_BITS, 2π periodic by exp

        for gate in gates {
            match gate {
                Gate::T(q) => pending_rz[q] += 1i64 << (PHASE_BITS - 2), // π/4
                Gate::S(q) => pending_rz[q] += 1i64 << (PHASE_BITS - 1), // π/2
                Gate::Z(q) => pending_rz[q] += 1i64 << PHASE_BITS,       // π
                Gate::RZ(q, p) => pending_rz[q] += p as i64,
                Gate::H(q) => {
                    Self::flush(&mut optimized, &mut pending_rz, q);
                    optimized.push(Gate::H(q));
                }
                Gate::CZ(a, b) => {
                    Self::flush(&mut optimized, &mut pending_rz, a);
                    Self::flush(&mut optimized, &mut pending_rz, b);
                    optimized.push(Gate::CZ(a, b));
                }
                Gate::MCZ(ctrls) => {
                    for &q in &ctrls {
                        Self::flush(&mut optimized, &mut pending_rz, q);
                    }
                    optimized.push(Gate::MCZ(ctrls));
                }
            }
        }

        for q in 0..num_qubits {
            Self::flush(&mut optimized, &mut pending_rz, q);
        }

        optimized
    }

    #[inline(always)]
    fn flush(opt: &mut Vec<Gate>, pending_rz: &mut [i64], q: usize) {
        let count = pending_rz[q];
        if count == 0 {
            return;
        }
        pending_rz[q] = 0;

        // reduce mod 2^(PHASE_BITS+1)
        let norm = phase_from_u32(count as u32);
        if norm != 0 {
            opt.push(Gate::RZ(q, norm));
        }
    }
}

// --- 3. REDUCTION & SUMMATION ENGINE ---

#[derive(Clone, Debug)]
enum DicksonOp {
    Swap(usize, usize),
    Add(usize, usize),
}
struct DicksonPlan {
    ops: Vec<DicksonOp>,
    reduced_b: Vec<FixedBitSet>,
    rank: usize,
}

fn plan_dickson_z4(mut b: Vec<FixedBitSet>, m: usize) -> DicksonPlan {
    let mut ops = Vec::new();
    let (mut r, mut p) = (0, 0);
    while p + 1 < m {
        let mut pivot = None;
        'o: for i in p..m {
            for j in (i + 1)..m {
                if b[i].contains(j) {
                    pivot = Some((i, j));
                    break 'o;
                }
            }
        }
        let Some((i, j)) = pivot else { break; };
        if i != p {
            b.swap(p, i);
            ops.push(DicksonOp::Swap(p, i));
        }
        let j_act = if j == p { i } else { j };
        if j_act != p + 1 {
            b.swap(p + 1, j_act);
            ops.push(DicksonOp::Swap(p + 1, j_act));
        }
        let rp = b[p].clone();
        let rp1 = b[p + 1].clone();
        for k in (p + 2)..m {
            if b[k].contains(p) {
                b[k].symmetric_difference_with(&rp1);
                ops.push(DicksonOp::Add(k, p + 1));
            }
            if b[k].contains(p + 1) {
                b[k].symmetric_difference_with(&rp);
                ops.push(DicksonOp::Add(k, p));
            }
        }
        r += 2;
        p += 2;
    }
    DicksonPlan { ops, reduced_b: b, rank: r }
}

#[inline(always)]
fn apply_plan_z4(plan: &DicksonPlan, v: &mut [u8]) {
    for op in &plan.ops {
        match *op {
            DicksonOp::Swap(i, j) => v.swap(i, j),
            DicksonOp::Add(target, pivot) => v[pivot] = z4_add(v[pivot], v[target]),
        }
    }
}

fn eval_sum_canonical(plan: &DicksonPlan, v: &[u8]) -> Complex64 {
    let mut sum_val = Complex64::new(1.0, 0.0);
    let (b, r, m) = (&plan.reduced_b, plan.rank, v.len());
    let mut p = 0;
    while p < r {
        let mut pair_sum = Complex64::new(0.0, 0.0);
        let has_edge = b[p].contains(p + 1) as u8;
        for x1 in 0..=1u8 {
            for x2 in 0..=1u8 {
                let mut ph = 0u8;
                if has_edge != 0 && x1 != 0 && x2 != 0 {
                    ph = z4_add(ph, 2);
                }
                if x1 != 0 {
                    ph = z4_add(ph, v[p]);
                }
                if x2 != 0 {
                    ph = z4_add(ph, v[p + 1]);
                }
                pair_sum += match ph {
                    0 => Complex64::new(1.0, 0.0),
                    1 => Complex64::new(0.0, 1.0),
                    2 => Complex64::new(-1.0, 0.0),
                    3 => Complex64::new(0.0, -1.0),
                    _ => unreachable!(),
                };
            }
        }
        sum_val *= pair_sum;
        p += 2;
    }
    for k in r..m {
        match v[k] & 3 {
            0 => sum_val *= 2.0,
            1 => sum_val *= Complex64::new(1.0, 1.0),
            2 => return Complex64::new(0.0, 0.0),
            3 => sum_val *= Complex64::new(1.0, -1.0),
            _ => unreachable!(),
        }
    }
    sum_val
}

// --- 4. OPTIMIZED SIMULATION ENGINE ---

/// Gray-code incremental evaluator for remainder terms (poly.rem).
struct RemGrayTracker {
    terms_by_bit: Vec<Vec<usize>>,
    need: Vec<u16>,
    weight: Vec<Phase>,
    sat_count: Vec<u16>,
    rem_phase: Phase,
}

impl RemGrayTracker {
    fn new(poly: &CompiledPhasePoly, vpos_of_var: &[i32], nv: usize, x_full: &[u8]) -> Self {
        let mut terms_by_bit: Vec<Vec<usize>> = (0..nv).map(|_| Vec::new()).collect();
        let mut need: Vec<u16> = Vec::with_capacity(poly.rem.len());
        let mut weight: Vec<Phase> = Vec::with_capacity(poly.rem.len());
        let mut rem_phase: Phase = 0;

        for term in &poly.rem {
            let mut ok_fixed = true;
            let mut positions: Vec<usize> = Vec::new();

            for &var in &term.vars {
                let pos = vpos_of_var[var];
                if pos >= 0 {
                    positions.push(pos as usize);
                } else if x_full[var] == 0 {
                    ok_fixed = false;
                    break;
                }
            }
            if !ok_fixed {
                continue;
            }

            if positions.len() > 1 {
                positions.sort_unstable();
                positions.dedup();
            }

            let tid = need.len();
            let deg = positions.len() as u16;
            need.push(deg);
            weight.push(term.weight);

            if deg == 0 {
                rem_phase = phase_add(rem_phase, term.weight);
            } else {
                for &b in &positions {
                    terms_by_bit[b].push(tid);
                }
            }
        }

        let sat_count = vec![0u16; need.len()];
        Self { terms_by_bit, need, weight, sat_count, rem_phase }
    }

    #[inline(always)]
    fn current_phase(&self) -> Phase {
        self.rem_phase
    }

    #[inline(always)]
    fn flip_bit(&mut self, flip: usize, turning_on: bool) {
        let term_list = &self.terms_by_bit[flip];
        for &t in term_list {
            let was_active = self.sat_count[t] == self.need[t];
            if turning_on {
                self.sat_count[t] += 1;
            } else {
                self.sat_count[t] -= 1;
            }
            let now_active = self.sat_count[t] == self.need[t];

            if !was_active && now_active {
                self.rem_phase = phase_add(self.rem_phase, self.weight[t]);
            } else if was_active && !now_active {
                self.rem_phase = phase_sub(self.rem_phase, self.weight[t]);
            }
        }
    }
}

pub fn amplitude_clifford_t_accel(poly: &CompiledPhasePoly, input: &[u8], target_y: usize) -> Complex64 {
    let t = poly.num_vars;

    // fixed assignments for boundary conditions
    let mut fixed: Vec<Option<u8>> = vec![None; t];
    for i in 0..poly.num_qubits {
        fixed[i] = Some(input[i]);
    }
    for i in 0..poly.num_qubits {
        let ov = poly.output_vars[i];
        let bit = ((target_y >> i) & 1) as u8;
        match fixed[ov] {
            None => fixed[ov] = Some(bit),
            Some(v) if v != bit => return Complex64::new(0.0, 0.0),
            _ => {}
        }
    }

    // compute x_full for vars fixed by boundary conditions
    let mut x_full = vec![0u8; t];
    for i in 0..t {
        if let Some(v) = fixed[i] {
            x_full[i] = v;
        }
    }

    // Mark all unfixed vars that appear in remainder terms as vvars
    let mut is_vvar = vec![false; t];
    for term in &poly.rem {
        for &v in &term.vars {
            if fixed[v].is_none() {
                is_vvar[v] = true;
            }
        }
    }

    let mut vvars = Vec::new();
    let mut uvars = Vec::new();
    for i in 0..t {
        if fixed[i].is_none() {
            if is_vvar[i] {
                vvars.push(i);
            } else {
                uvars.push(i);
            }
        }
    }

    let (nv, nu) = (vvars.len(), uvars.len());

    // Fast map: original var idx -> vvars position or -1
    let mut var_to_vpos: Vec<i32> = vec![-1; t];
    for (pos, &var) in vvars.iter().enumerate() {
        var_to_vpos[var] = pos as i32;
    }

    // eps_base_z4 depends only on fixed vars
    let mut eps_base_z4 = poly.eps4 & 3;
    let mut f_list = Vec::new();
    for (idx, &v) in fixed.iter().enumerate() {
        if v == Some(1) {
            f_list.push(idx);
        }
    }
    for &f in &f_list {
        eps_base_z4 = z4_add(eps_base_z4, poly.v4[f] & 3);
        for &f2 in &f_list {
            if f < f2 && poly.b4[f].contains(f2) {
                eps_base_z4 = z4_add(eps_base_z4, 2);
            }
        }
    }

    // Build bu (quadratic adjacency among uvars), vu_base, and cross masks (uvar -> vvar indices)
    let mut bu: Vec<FixedBitSet> = (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();
    let mut vu_base = vec![0u8; nu];
    let mut cross_masks = vec![FixedBitSet::with_capacity(nv); nu];

    for (ui, &orig_u) in uvars.iter().enumerate() {
        vu_base[ui] = poly.v4[orig_u] & 3;

        // u-u quadratic
        for (uj, &orig_uj) in uvars.iter().enumerate() {
            if ui < uj && poly.b4[orig_u].contains(orig_uj) {
                bu[ui].insert(uj);
                bu[uj].insert(ui);
            }
        }

        // fixed-u interactions shift linear term by +2 (mod 4)
        for &f in &f_list {
            if poly.b4[orig_u].contains(f) {
                vu_base[ui] = z4_add(vu_base[ui], 2);
            }
        }

        // u-v cross interactions for fast updates when v flips
        for (vj, &orig_v) in vvars.iter().enumerate() {
            if poly.b4[orig_u].contains(orig_v) {
                cross_masks[ui].insert(vj);
            }
        }
    }

    let plan = plan_dickson_z4(bu, nu);

    let mut rem_tracker = RemGrayTracker::new(poly, &var_to_vpos, nv, &x_full);

    let mut amp = Complex64::new(0.0, 0.0);
    let mut cur_vu = vu_base.clone();
    let mut cur_eps = eps_base_z4;
    let mut cur_x_vvars = vec![0u8; nv];
    let mut vu_exec = vec![0u8; nu];

    let iters = 1usize << nv;
    for i in 0..iters {
        vu_exec.copy_from_slice(&cur_vu);
        apply_plan_z4(&plan, &mut vu_exec);

        let rem_phase: Phase = rem_tracker.current_phase();
        let eps_phase: Phase = phase_from_u32((cur_eps as u32) << (PHASE_BITS - 1));
        let total_phase = phase_add(rem_phase, eps_phase);

        let phase = Complex64::from_polar(1.0, phase_to_angle_rad(total_phase));
        amp += phase * eval_sum_canonical(&plan, &vu_exec);

        if i + 1 < iters {
            let flip = (i + 1).trailing_zeros() as usize;

            for ui in 0..nu {
                if cross_masks[ui].contains(flip) {
                    cur_vu[ui] = z4_add(cur_vu[ui], 2);
                }
            }

            let v_idx = vvars[flip];
            let turning_on = cur_x_vvars[flip] == 0;

            if turning_on {
                cur_eps = z4_add(cur_eps, poly.v4[v_idx] & 3);

                for &f in &f_list {
                    if poly.b4[f].contains(v_idx) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }
                for (j, &ov) in vvars.iter().enumerate() {
                    if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }

                cur_x_vvars[flip] = 1;
            } else {
                cur_eps = z4_sub(cur_eps, poly.v4[v_idx] & 3);

                for &f in &f_list {
                    if poly.b4[f].contains(v_idx) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }
                for (j, &ov) in vvars.iter().enumerate() {
                    if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }

                cur_x_vvars[flip] = 0;
            }

            rem_tracker.flip_bit(flip, turning_on);
        }
    }

    amp * (2f64).powf(-(poly.num_h as f64) / 2.0)
}

pub fn simulate_statevector(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    (0..(1usize << poly.num_qubits))
        .map(|y| amplitude_clifford_t_accel(poly, input, y))
        .collect()
}

#[inline(always)]
fn t_phase_unit() -> Phase {
    // π/4 = 2^(PHASE_BITS-2) units
    phase_from_u32(1u32 << (PHASE_BITS - 2))
}

#[inline(always)]
fn pi_phase_unit() -> Phase {
    // π = 2^PHASE_BITS units
    phase_from_u32(1u32 << PHASE_BITS)
}

pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    let mut wire: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut next_var = n;
    let mut num_h = 0usize;

    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];

    let mut rem: Vec<PhaseTerm> = Vec::new();

    let mut grow = |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
        while v4.len() < new_t {
            v4.push(0);
        }
        while b4.len() < new_t {
            b4.push(FixedBitSet::with_capacity(new_t));
        }
        for row in b4.iter_mut() {
            row.grow(new_t);
        }
    };

    for g in gates {
        match g {
            Gate::H(q) => {
                let prev = *wire[*q].last().unwrap();
                let cur = next_var;
                next_var += 1;
                num_h += 1;
                wire[*q].push(cur);
                grow(next_var, &mut b4, &mut v4);
                b4[prev].insert(cur);
                b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[*a].last().unwrap();
                let vb = *wire[*b].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                b4[va].insert(vb);
                b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[*q].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                v4[v] = z4_add(v4[v], 2);
            }
            Gate::S(q) => {
                let v = *wire[*q].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                v4[v] = z4_add(v4[v], 1);
            }
            Gate::T(q) => {
                let v = *wire[*q].last().unwrap();
                rem.push(PhaseTerm { weight: t_phase_unit(), vars: vec![v] });
            }
            Gate::RZ(q, phase) => {
                let v = *wire[*q].last().unwrap();
                if *phase != 0 {
                    rem.push(PhaseTerm { weight: *phase, vars: vec![v] });
                }
            }
            Gate::MCZ(ctrls) => {
                match ctrls.len() {
                    0 => {
                        // global -1
                        rem.push(PhaseTerm { weight: pi_phase_unit(), vars: vec![] });
                    }
                    1 => {
                        // Z
                        let v = *wire[ctrls[0]].last().unwrap();
                        grow(next_var, &mut b4, &mut v4);
                        v4[v] = z4_add(v4[v], 2);
                    }
                    2 => {
                        // CZ
                        let va = *wire[ctrls[0]].last().unwrap();
                        let vb = *wire[ctrls[1]].last().unwrap();
                        grow(next_var, &mut b4, &mut v4);
                        b4[va].insert(vb);
                        b4[vb].insert(va);
                    }
                    _ => {
                        // high-degree -1 as remainder AND-monomial with weight π
                        let mut vars = Vec::with_capacity(ctrls.len());
                        for &q in ctrls {
                            vars.push(*wire[q].last().unwrap());
                        }
                        rem.push(PhaseTerm { weight: pi_phase_unit(), vars });
                    }
                }
            }
        }
    }

    CompiledPhasePoly {
        num_qubits: n,
        num_vars: next_var,
        num_h,
        output_vars: wire.iter().map(|w| *w.last().unwrap()).collect(),
        b4,
        v4,
        eps4: 0,
        rem,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn approx_eq(a: Complex64, b: Complex64, eps: f64) {
        assert!(
            (a - b).norm() <= eps,
            "a={:?} b={:?} |a-b|={}",
            a,
            b,
            (a - b).norm()
        );
    }

    #[test]
    fn test_mcz_on_three_qubits_adds_pi_phase_only_on_all_ones() {
        let gates = vec![Gate::H(0), Gate::H(1), Gate::H(2), Gate::MCZ(vec![0, 1, 2])];
        let poly = compile_clifford_t(3, &gates);
        let input = vec![0u8, 0u8, 0u8];

        let mag = (1.0 / 8.0f64).sqrt();
        for y in 0..8usize {
            let a = amplitude_clifford_t_accel(&poly, &input, y);
            if y == 7 {
                approx_eq(a, c(-mag, 0.0), 1e-10);
            } else {
                approx_eq(a, c(mag, 0.0), 1e-10);
            }
        }
    }

    #[test]
    fn test_mcz_equiv_to_cz_for_two_controls() {
        let gates_cz = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let gates_mcz = vec![Gate::H(0), Gate::H(1), Gate::MCZ(vec![0, 1]), Gate::H(1)];

        let p1 = compile_clifford_t(2, &gates_cz);
        let p2 = compile_clifford_t(2, &gates_mcz);

        let input = vec![0u8, 0u8];
        for y in 0..4usize {
            let a1 = amplitude_clifford_t_accel(&p1, &input, y);
            let a2 = amplitude_clifford_t_accel(&p2, &input, y);
            approx_eq(a1, a2, 1e-10);
        }
    }

    #[test]
    fn test_rz_equiv_to_t() {
        let gates_t = vec![Gate::H(0), Gate::T(0), Gate::H(0)];
        let gates_rz = vec![Gate::H(0), Gate::RZ(0, phase_from_u32(1u32 << (PHASE_BITS - 2))), Gate::H(0)];

        let p1 = compile_clifford_t(1, &gates_t);
        let p2 = compile_clifford_t(1, &gates_rz);

        let input = vec![0u8];
        for y in 0..2usize {
            let a1 = amplitude_clifford_t_accel(&p1, &input, y);
            let a2 = amplitude_clifford_t_accel(&p2, &input, y);
            approx_eq(a1, a2, 1e-10);
        }
    }
}