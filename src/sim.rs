use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::f64::consts::PI;

pub use crate::qc::{read_qasm_file, write_qasm_file, write_qasm_string, QasmError};

/// Phase ring configuration for Clifford+T simulation.
///
/// Phases are represented in units of π / 2^PHASE_BITS.
/// A stored integer `p` corresponds to an angle θ = π * p / 2^PHASE_BITS.
/// This allows efficient encoding of dyadic rotations and Clifford+T gates.
///
/// # Theory
/// Clifford+T circuits can be simulated efficiently by representing quantum gates
/// as phase polynomials over Boolean variables. The Clifford group is handled
/// using Z4 arithmetic, while non-Clifford gates (T, RZ) are tracked as remainder terms
/// in a dyadic phase ring. This enables fast amplitude and statevector computation.
pub const PHASE_BITS: u32 = 16;
pub type Phase = u16;

/// Returns the modulus for the phase ring (2^PHASE_BITS).
#[inline(always)]
pub const fn phase_modulus() -> u32 {
    1u32 << PHASE_BITS
}

/// Returns the mask for the phase ring (2^PHASE_BITS - 1).
#[inline(always)]
pub const fn phase_mask_u32() -> u32 {
    (1u32 << PHASE_BITS) - 1
}

/// Modular addition in the phase ring.
#[inline(always)]
fn phase_add(a: Phase, b: Phase) -> Phase {
    a.wrapping_add(b)
}

/// Modular subtraction in the phase ring.
#[inline(always)]
fn phase_sub(a: Phase, b: Phase) -> Phase {
    a.wrapping_sub(b)
}

/// Converts a u32 to a phase value, masking to PHASE_BITS.
#[inline(always)]
fn phase_from_u32(x: u32) -> Phase {
    (x & phase_mask_u32()) as Phase
}

/// Converts a phase value to a real angle in radians.
#[inline(always)]
fn phase_to_angle_rad(p: Phase) -> f64 {
    PI * (p as f64) / (phase_modulus() as f64)
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
/// - RZ(q, phase): Arbitrary Z rotation by phase units
/// - CZ(a, b): Controlled-Z gate between qubits a and b
#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),
    Z(usize),
    S(usize),
    T(usize),
    RZ(usize, Phase),
    CZ(usize, usize),
}

/// Remainder term for non-Clifford phase tracking.
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
///
/// Supports builder methods for all Clifford+T gates.
pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    /// Create a new circuit with the given number of qubits.
    pub fn new(num_qubits: usize) -> Self {
        Self { num_qubits, gates: Vec::new() }
    }
    /// Add a Hadamard gate.
    pub fn h(&mut self, q: usize) {
        self.gates.push(Gate::H(q));
    }
    /// Add a Pauli-Z gate.
    pub fn z(&mut self, q: usize) {
        self.gates.push(Gate::Z(q));
    }
    /// Add a phase (S) gate.
    pub fn s(&mut self, q: usize) {
        self.gates.push(Gate::S(q));
    }
    /// Add a T gate.
    pub fn t(&mut self, q: usize) {
        self.gates.push(Gate::T(q));
    }
    /// Add an arbitrary RZ rotation.
    pub fn rz(&mut self, q: usize, phase: Phase) {
        self.gates.push(Gate::RZ(q, phase));
    }
    /// Add a controlled-Z gate.
    pub fn cz(&mut self, q1: usize, q2: usize) {
        self.gates.push(Gate::CZ(q1, q2));
    }

    /// Compile the circuit into a phase polynomial for simulation.
    pub fn compile(&self) -> CompiledPhasePoly {
        let optimized_gates = Transpiler::optimize(self.num_qubits, self.gates.clone());
        compile_clifford_t(self.num_qubits, &optimized_gates)
    }
}

// --- 2. MULTI-PASS TRANSPILER (PHASE TELEPORTATION) ---
//
// This pass coalesces Z/S/T/RZ phases on a wire until they hit an H or CZ.
// For speed, we avoid HashMap and use a fixed Vec<i32> sized by num_qubits.
// Note: this transpiler only tracks phases by qubit index for the *gate list*,
// not the internal phase polynomial vars (that's handled by compilation).

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
        let mut pending: Vec<i32> = vec![0; num_qubits]; // in units of π/4 (Z8-like) for fast normalization
        let mut pending_rz: Vec<i32> = vec![0; num_qubits]; // full precision units (π / 2^PHASE_BITS)

        // We can safely keep a single full-precision accumulator and flush to RZ.
        // Keep legacy Z/S/T folded as exact dyadic phases.
        for gate in gates {
            match gate {
                Gate::T(q) => {
                    pending_rz[q] += 1i32 << (PHASE_BITS - 2); // π/4
                }
                Gate::S(q) => {
                    pending_rz[q] += 1i32 << (PHASE_BITS - 1); // π/2
                }
                Gate::Z(q) => {
                    pending_rz[q] += 1i32 << PHASE_BITS; // π (mod 2π in this representation is 2^(PHASE_BITS+1), but RZ phases are mod 2π? See below.)
                    // We interpret RZ as diag(1, e^{iθ}); θ periodic mod 2π.
                    // Our phase ring here is mod 2^PHASE_BITS for π-units; that means periodic mod 2π corresponds to doubling.
                    // To keep consistent with the rest of the simulator (which uses ω = exp(iπ/2^PHASE_BITS)),
                    // Z is phase = 2^PHASE_BITS (i.e. θ=π) which reduces mod 2^(PHASE_BITS+1) in 2π periodicity.
                    // But since we only ever exponentiate exp(i * π * p / 2^PHASE_BITS), periodicity is p mod 2^(PHASE_BITS+1).
                    // For speed and compatibility with the original Z8 model (mod 2π handled by Complex exponential),
                    // we store phases mod 2^(PHASE_BITS+1) internally in pending_rz, and reduce when converting.
                }
                Gate::RZ(q, p) => {
                    pending_rz[q] += p as i32;
                }
                Gate::H(q) => {
                    Self::flush(&mut optimized, &mut pending_rz, q);
                    optimized.push(Gate::H(q));
                }
                Gate::CZ(a, b) => {
                    Self::flush(&mut optimized, &mut pending_rz, a);
                    Self::flush(&mut optimized, &mut pending_rz, b);
                    optimized.push(Gate::CZ(a, b));
                }
            }
        }

        for q in 0..num_qubits {
            Self::flush(&mut optimized, &mut pending_rz, q);
        }

        optimized
    }

    #[inline(always)]
    fn flush(opt: &mut Vec<Gate>, pending_rz: &mut [i32], q: usize) {
        let count = pending_rz[q];
        if count == 0 {
            return;
        }
        pending_rz[q] = 0;

        // Reduce mod 2^PHASE_BITS (we intentionally drop 2π periodicity nuances; exp handles it anyway).
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
            DicksonOp::Add(target, pivot) => {
                // v[pivot] = (v[pivot] + v[target]) % 4
                v[pivot] = z4_add(v[pivot], v[target]);
            }
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
        // 4 cases (x1,x2) in {0,1}^2
        for x1 in 0..=1u8 {
            for x2 in 0..=1u8 {
                // ph = (2*x1*x2*edge + v[p]*x1 + v[p+1]*x2) mod 4
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

#[derive(Clone, Debug)]
struct RemTermMask {
    weight: Phase,
    /// mask over vvars (only valid when nv <= 64)
    mask: u64,
    /// vars that are NOT in vvars (fixed/full vars), must be checked via x_full
    other_vars: Vec<usize>,
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

    // Identify which vars appear in non-Clifford remainder terms.
    // For maximum speed, we only enumerate over those vars + any unfixed vars that remainder depends on.
    let mut in_rem = vec![false; t];
    for term in &poly.rem {
        for &v in &term.vars {
            in_rem[v] = true;
        }
    }

    // Variable compression split:
    // - vvars: appear in remainder and are unfixed
    // - uvars: unfixed but not in remainder
    let mut vvars = Vec::new();
    let mut uvars = Vec::new();
    for i in 0..t {
        if fixed[i].is_none() {
            if in_rem[i] {
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

    // Precompute Dickson plan once (major optimization)
    let plan = plan_dickson_z4(bu, nu);

    // Prepare remainder terms:
    // If nv <= 64, encode vvar-dependence as masks for very fast activation checks.
    let use_mask_terms = nv <= 64;
    let mut rem_masks: Vec<RemTermMask> = Vec::new();

    if use_mask_terms {
        rem_masks.reserve(poly.rem.len());
        for term in &poly.rem {
            let mut mask: u64 = 0;
            let mut other_vars = Vec::new();
            for &v in &term.vars {
                let pos = var_to_vpos[v];
                if pos >= 0 {
                    mask |= 1u64 << (pos as u64);
                } else {
                    other_vars.push(v);
                }
            }
            rem_masks.push(RemTermMask { weight: term.weight, mask, other_vars });
        }
    }

    // State for Gray code enumeration over vvars
    let mut amp = Complex64::new(0.0, 0.0);

    // Z4 linear vector updated incrementally
    let mut cur_vu = vu_base.clone();
    let mut cur_eps = eps_base_z4;

    // vvar assignment mask in Gray code (only valid if nv<=64) and also as bytes for fallback paths
    let mut cur_mask_u64: u64 = 0;
    let mut cur_x_vvars = vec![0u8; nv];

    // Base remainder phase contributed by fixed vars:
    // Since remainder terms are AND-monomials, a term contributes weight iff all its vars are 1.
    // If the term depends only on fixed vars, it is constant across enumeration.
    let mut rem_base: Phase = 0;
    if use_mask_terms {
        for rt in &rem_masks {
            // term active if:
            //  - its vvar mask is subset of cur mask (currently 0, so only mask==0 qualifies)
            //  - all other_vars are 1 in x_full
            if rt.mask == 0 {
                let mut ok = true;
                for &ov in &rt.other_vars {
                    if x_full[ov] == 0 {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    rem_base = phase_add(rem_base, rt.weight);
                }
            }
        }
    } else {
        // conservative baseline: nothing precomputed
    }

    // Working buffer to avoid cloning cur_vu every iteration
    let mut vu_exec = vec![0u8; nu];

    let iters = 1usize << nv;
    for i in 0..iters {
        // copy current vu into exec buffer and apply plan
        vu_exec.copy_from_slice(&cur_vu);
        apply_plan_z4(&plan, &mut vu_exec);

        // remainder evaluation
        let mut rem_phase: Phase = rem_base;

        if use_mask_terms {
            // add only terms with non-empty vvar masks, active in current mask
            for rt in &rem_masks {
                if rt.mask == 0 {
                    continue; // already included in base
                }
                if (rt.mask & cur_mask_u64) != rt.mask {
                    continue;
                }
                // check other fixed vars
                let mut ok = true;
                for &ov in &rt.other_vars {
                    if x_full[ov] == 0 {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    rem_phase = phase_add(rem_phase, rt.weight);
                }
            }
        } else {
            // fallback: scan terms and check vars against cur_x_vvars or x_full
            for term in &poly.rem {
                let mut ok = true;
                for &v in &term.vars {
                    let pos = var_to_vpos[v];
                    if pos >= 0 {
                        if cur_x_vvars[pos as usize] == 0 {
                            ok = false;
                            break;
                        }
                    } else if x_full[v] == 0 {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    rem_phase = phase_add(rem_phase, term.weight);
                }
            }
        }

        // Combine remainder (Z_{2^PHASE_BITS}) with Clifford Z4 phase:
        // Original code used exp(i * π * ((r_z8 + 2*eps4) mod 8) / 4).
        // Here we generalize:
        //
        // - rem_phase is in units π/2^PHASE_BITS.
        // - eps4 is in Z4, representing multiples of π/2 (because original "4*eps4 mod 8" -> π/2 steps).
        //   In angle terms, eps4 contributes angle = π * (2*eps4)/4 = π*eps4/2.
        //   Convert that to our units: π*eps4/2 = π * (eps4 * 2^(PHASE_BITS-1)) / 2^PHASE_BITS.
        let eps_phase: Phase = phase_from_u32((cur_eps as u32) << (PHASE_BITS - 1));

        let total_phase = phase_add(rem_phase, eps_phase);
        let phase = Complex64::from_polar(1.0, phase_to_angle_rad(total_phase));
        amp += phase * eval_sum_canonical(&plan, &vu_exec);

        // Gray-code update
        if i + 1 < iters {
            let flip = (i + 1).trailing_zeros() as usize;

            // toggle mask bit
            let bit = 1u64 << flip;
            let bit_set = (cur_mask_u64 & bit) != 0;

            // update Z4 linear vector for uvars impacted by this vvar flip
            for ui in 0..nu {
                if cross_masks[ui].contains(flip) {
                    // flipping a v var toggles contribution by +2 (mod 4) to vu[ui]
                    cur_vu[ui] = z4_add(cur_vu[ui], 2);
                }
            }

            // update eps4 due to v-v and fixed-v interactions (Z4 part)
            let v_idx = vvars[flip];

            if !bit_set {
                // turning v bit 0->1
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

                cur_mask_u64 |= bit;
                cur_x_vvars[flip] = 1;
            } else {
                // turning v bit 1->0
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

                cur_mask_u64 &= !bit;
                cur_x_vvars[flip] = 0;
            }
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
    // π/4 = π * 2^(PHASE_BITS-2) / 2^PHASE_BITS
    phase_from_u32(1u32 << (PHASE_BITS - 2))
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
        match *g {
            Gate::H(q) => {
                let prev = *wire[q].last().unwrap();
                let cur = next_var;
                next_var += 1;
                num_h += 1;
                wire[q].push(cur);
                grow(next_var, &mut b4, &mut v4);
                b4[prev].insert(cur);
                b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[a].last().unwrap();
                let vb = *wire[b].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                b4[va].insert(vb);
                b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[q].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                v4[v] = z4_add(v4[v], 2);
            }
            Gate::S(q) => {
                let v = *wire[q].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                v4[v] = z4_add(v4[v], 1);
            }
            Gate::T(q) => {
                let v = *wire[q].last().unwrap();
                rem.push(PhaseTerm { weight: t_phase_unit(), vars: vec![v] });
            }
            Gate::RZ(q, phase) => {
                let v = *wire[q].last().unwrap();
                if phase != 0 {
                    rem.push(PhaseTerm { weight: phase, vars: vec![v] });
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
    fn test_simulate_statevector() {
        let gates = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly = compile_clifford_t(2, &gates);
        let input = vec![0u8, 0u8];
        let statevec = simulate_statevector(&poly, &input);
        let expected_output = vec![
            Complex64::new(0.5f64.sqrt(), 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.5f64.sqrt(), 0.0),
        ];
        for (amp, expected) in statevec.iter().zip(expected_output.iter()) {
            assert!(
                (amp - expected).norm() < 1e-10,
                "amp={:?} expected={:?}",
                amp,
                expected
            );
        }
    }

    #[test]
    fn test_bell_state_clifford_amplitudes() {
        let gates = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly = compile_clifford_t(2, &gates);
        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);
        let a01 = amplitude_clifford_t_accel(&poly, &input, 1);
        let a10 = amplitude_clifford_t_accel(&poly, &input, 2);

        let norm = (0.5f64).sqrt();
        approx_eq(a00, c(norm, 0.0), 1e-10);
        approx_eq(a11, c(norm, 0.0), 1e-10);
        approx_eq(a01, c(0.0, 0.0), 1e-10);
        approx_eq(a10, c(0.0, 0.0), 1e-10);
    }

    #[test]
    fn test_clifford_plus_s_phase() {
        let gates = vec![Gate::H(0), Gate::S(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly = compile_clifford_t(2, &gates);
        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);

        let norm = (0.5f64).sqrt();
        approx_eq(a00, c(norm, 0.0), 1e-10);
        approx_eq(a11, c(0.0, norm), 1e-10);
    }

    #[test]
    fn test_clifford_plus_t_phase() {
        let gates = vec![Gate::H(0), Gate::T(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly = compile_clifford_t(2, &gates);
        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);

        approx_eq(a00, c((0.5f64).sqrt(), 0.0), 1e-10);
        approx_eq(a11, c(0.5, 0.5), 1e-10);
    }

    #[test]
    fn test_three_qubit_mixed_norm_sanity() {
        let gates = vec![
            Gate::H(0),
            Gate::CZ(0, 1),
            Gate::T(1),
            Gate::H(2),
            Gate::CZ(1, 2),
            Gate::S(0),
            Gate::Z(2),
        ];
        let poly = compile_clifford_t(3, &gates);
        let input = vec![0u8, 0u8, 0u8];

        let mut norm_sq = 0.0f64;
        for y in 0..(1usize << 3) {
            let a = amplitude_clifford_t_accel(&poly, &input, y);
            norm_sq += a.norm_sqr();
        }
        assert!((norm_sq - 1.0).abs() < 1e-8, "norm_sq={}", norm_sq);
    }

    #[test]
    fn test_rz_equiv_to_t() {
        // RZ(π/4) should behave like T in this remainder-only model
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