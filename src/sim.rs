use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::collections::HashMap;
use std::f64::consts::PI;

pub use crate::qc::{read_qasm_file, write_qasm_file, write_qasm_string, QasmError};

// --- 1. CORE DATA STRUCTURES ---

#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize), Z(usize), S(usize), T(usize), CZ(usize, usize),
}

#[derive(Clone, Debug)]
pub struct Z8Term {
    pub weight: u8,
    pub vars: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct CompiledPhasePoly {
    pub num_qubits: usize,
    pub num_vars: usize,
    pub num_h: usize,
    pub output_vars: Vec<usize>,
    pub b4: Vec<FixedBitSet>,
    pub v4: Vec<u8>,
    pub eps4: u8,
    pub rem: Vec<Z8Term>,
}

pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self { Self { num_qubits, gates: Vec::new() } }
    pub fn h(&mut self, q: usize) { self.gates.push(Gate::H(q)); }
    pub fn z(&mut self, q: usize) { self.gates.push(Gate::Z(q)); }
    pub fn s(&mut self, q: usize) { self.gates.push(Gate::S(q)); }
    pub fn t(&mut self, q: usize) { self.gates.push(Gate::T(q)); }
    pub fn cz(&mut self, q1: usize, q2: usize) { self.gates.push(Gate::CZ(q1, q2)); }

    pub fn compile(&self) -> CompiledPhasePoly {
        let optimized_gates = Transpiler::optimize(self.gates.clone());
        compile_clifford_t(self.num_qubits, &optimized_gates)
    }
}

// --- 2. MULTI-PASS TRANSPILER (PHASE TELEPORTATION) ---

pub struct Transpiler;

impl Transpiler {
    pub fn optimize(gates: Vec<Gate>) -> Vec<Gate> {
        let mut current_gates = gates;
        loop {
            let start_count = current_gates.len();
            current_gates = Self::phase_teleportation_pass(current_gates);
            if current_gates.len() >= start_count { break; }
        }
        current_gates
    }

    fn phase_teleportation_pass(gates: Vec<Gate>) -> Vec<Gate> {
        let mut optimized = Vec::with_capacity(gates.len());
        let mut pending: HashMap<usize, i8> = HashMap::new();

        for gate in gates {
            match gate {
                Gate::T(q) => { *pending.entry(q).or_insert(0) += 1; }
                Gate::S(q) => { *pending.entry(q).or_insert(0) += 2; }
                Gate::Z(q) => { *pending.entry(q).or_insert(0) += 4; }
                Gate::H(q) | Gate::CZ(q, _) => {
                    Self::flush(&mut optimized, &mut pending, q);
                    if let Gate::CZ(_, q2) = gate { Self::flush(&mut optimized, &mut pending, q2); }
                    optimized.push(gate);
                }
            }
        }
        let keys: Vec<usize> = pending.keys().cloned().collect();
        for q in keys { Self::flush(&mut optimized, &mut pending, q); }
        optimized
    }

    fn flush(opt: &mut Vec<Gate>, phases: &mut HashMap<usize, i8>, q: usize) {
        if let Some(count) = phases.remove(&q) {
            let norm = (count % 8 + 8) % 8;
            if norm % 2 != 0 { opt.push(Gate::T(q)); }
            if (norm / 2) % 2 != 0 { opt.push(Gate::S(q)); }
            if (norm / 4) % 2 != 0 { opt.push(Gate::Z(q)); }
        }
    }
}

// --- 3. REDUCTION & SUMMATION ENGINE ---

#[derive(Clone, Debug)]
enum DicksonOp { Swap(usize, usize), Add(usize, usize) }
struct DicksonPlan { ops: Vec<DicksonOp>, reduced_b: Vec<FixedBitSet>, rank: usize }

fn plan_dickson_z4(mut b: Vec<FixedBitSet>, m: usize) -> DicksonPlan {
    let mut ops = Vec::new();
    let (mut r, mut p) = (0, 0);
    while p + 1 < m {
        let mut pivot = None;
        'o: for i in p..m { for j in (i+1)..m { if b[i].contains(j) { pivot = Some((i, j)); break 'o; } } }
        let Some((i, j)) = pivot else { break; };
        if i != p { b.swap(p, i); ops.push(DicksonOp::Swap(p, i)); }
        let j_act = if j == p { i } else { j };
        if j_act != p + 1 { b.swap(p + 1, j_act); ops.push(DicksonOp::Swap(p + 1, j_act)); }
        let rp = b[p].clone(); let rp1 = b[p + 1].clone();
        for k in (p + 2)..m {
            if b[k].contains(p) { b[k].symmetric_difference_with(&rp1); ops.push(DicksonOp::Add(k, p+1)); }
            if b[k].contains(p + 1) { b[k].symmetric_difference_with(&rp); ops.push(DicksonOp::Add(k, p)); }
        }
        r += 2; p += 2;
    }
    DicksonPlan { ops, reduced_b: b, rank: r }
}

#[inline(always)]
fn apply_plan_z4(plan: &DicksonPlan, v: &mut [u8]) {
    for op in &plan.ops {
        match *op {
            DicksonOp::Swap(i, j) => v.swap(i, j),
            DicksonOp::Add(target, pivot) => { v[pivot] = (v[pivot] + v[target]) % 4; }
        }
    }
}

fn eval_sum_canonical(plan: &DicksonPlan, v: &[u8]) -> Complex64 {
    let mut sum_val = Complex64::new(1.0, 0.0);
    let (b, r, m) = (&plan.reduced_b, plan.rank, v.len());
    let mut p = 0;
    while p < r {
        let mut pair_sum = Complex64::new(0.0, 0.0);
        for x1 in 0..=1 {
            for x2 in 0..=1 {
                let ph = (2 * x1 * x2 * (b[p].contains(p + 1) as u8) + v[p] * x1 + v[p+1] * x2) % 4;
                pair_sum += match ph { 0 => Complex64::new(1.0, 0.0), 1 => Complex64::new(0.0, 1.0), 2 => Complex64::new(-1.0, 0.0), 3 => Complex64::new(0.0, -1.0), _ => unreachable!() };
            }
        }
        sum_val *= pair_sum; p += 2;
    }
    for k in r..m {
        match v[k] % 4 {
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

pub fn amplitude_clifford_t_accel(poly: &CompiledPhasePoly, input: &[u8], target_y: usize) -> Complex64 {
    let t = poly.num_vars;
    let mut fixed = vec![None; t];
    for i in 0..poly.num_qubits { fixed[i] = Some(input[i]); }
    for i in 0..poly.num_qubits {
        let ov = poly.output_vars[i];
        let bit = ((target_y >> i) & 1) as u8;
        match fixed[ov] {
            None => fixed[ov] = Some(bit),
            Some(v) if v != bit => return Complex64::new(0.0, 0.0),
            _ => {}
        }
    }

    let mut in_rem = vec![false; t];
    for term in &poly.rem { for &v in &term.vars { in_rem[v] = true; } }
    let (mut vvars, mut uvars) = (Vec::new(), Vec::new());
    // Variable Compression: Only compress variables NOT in rem.
    for i in 0..t { if fixed[i].is_none() { if in_rem[i] { vvars.push(i); } else { uvars.push(i); } } }

    let (nv, nu) = (vvars.len(), uvars.len());
    let mut x_full = vec![0u8; t];
    for i in 0..t { if let Some(v) = fixed[i] { x_full[i] = v; } }

    let mut eps_base_z4 = poly.eps4 % 4;
    let f_list: Vec<usize> = fixed.iter().enumerate().filter_map(|(idx, &v)| if v == Some(1) { Some(idx) } else { None }).collect();
    for &f in &f_list {
        eps_base_z4 = (eps_base_z4 + poly.v4[f]) % 4;
        for &f2 in &f_list { if f < f2 && poly.b4[f].contains(f2) { eps_base_z4 = (eps_base_z4 + 2) % 4; } }
    }

    let mut bu: Vec<FixedBitSet> = (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();
    let mut vu_base = vec![0u8; nu];
    let mut cross_masks = vec![FixedBitSet::with_capacity(nv); nu];
    for (ui, &orig_u) in uvars.iter().enumerate() {
        vu_base[ui] = poly.v4[orig_u];
        for (uj, &orig_uj) in uvars.iter().enumerate() { if ui < uj && poly.b4[orig_u].contains(orig_uj) { bu[ui].insert(uj); bu[uj].insert(ui); } }
        for &f in &f_list { if poly.b4[orig_u].contains(f) { vu_base[ui] = (vu_base[ui] + 2) % 4; } }
        for (vj, &orig_v) in vvars.iter().enumerate() { if poly.b4[orig_u].contains(orig_v) { cross_masks[ui].insert(vj); } }
    }

    let plan = plan_dickson_z4(bu, nu); // PRE-DICKSON
    let mut amp = Complex64::new(0.0, 0.0);
    let (mut cur_vu, mut cur_eps) = (vu_base.clone(), eps_base_z4);
    let mut cur_x_vvars = vec![0u8; nv];

    for i in 0..(1usize << nv) {
        let mut vu_exec = cur_vu.clone();
        apply_plan_z4(&plan, &mut vu_exec);
        
        // Remainder calculation using current vvar state.
        let mut r_z8 = 0u8;
        for term in &poly.rem {
            if term.vars.iter().all(|&v| {
                if let Some(p) = vvars.iter().position(|&x| x == v) { cur_x_vvars[p] == 1 } else { x_full[v] == 1 }
            }) { r_z8 = (r_z8 + term.weight) & 7; }
        }

        let phase = Complex64::from_polar(1.0, PI * ((r_z8 + (cur_eps * 2)) % 8) as f64 / 4.0);
        amp += phase * eval_sum_canonical(&plan, &vu_exec);

        if i + 1 < (1usize << nv) {
            let flip = (i + 1).trailing_zeros() as usize; // GRAY CODE
            let bit_set = cur_x_vvars[flip] == 1;
            let v_idx = vvars[flip];
            for ui in 0..nu { if cross_masks[ui].contains(flip) { cur_vu[ui] = (cur_vu[ui] + 2) % 4; } }
            if !bit_set {
                cur_eps = (cur_eps + poly.v4[v_idx]) % 4;
                for &f in &f_list { if poly.b4[f].contains(v_idx) { cur_eps = (cur_eps + 2) % 4; } }
                for (j, &ov) in vvars.iter().enumerate() { if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) { cur_eps = (cur_eps + 2) % 4; } }
                cur_x_vvars[flip] = 1;
            } else {
                cur_eps = (cur_eps + 4 - poly.v4[v_idx]) % 4;
                for &f in &f_list { if poly.b4[f].contains(v_idx) { cur_eps = (cur_eps + 2) % 4; } }
                for (j, &ov) in vvars.iter().enumerate() { if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) { cur_eps = (cur_eps + 2) % 4; } }
                cur_x_vvars[flip] = 0;
            }
        }
    }
    amp * (2f64).powf(-(poly.num_h as f64) / 2.0)
}

pub fn simulate_statevector(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    (0..(1usize << poly.num_qubits)).map(|y| amplitude_clifford_t_accel(poly, input, y)).collect()
}

pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    let mut wire: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut next_var = n;
    let mut num_h = 0usize;
    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];
    let mut rem: Vec<Z8Term> = Vec::new();

    let mut grow = |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
        while v4.len() < new_t { v4.push(0); }
        while b4.len() < new_t { b4.push(FixedBitSet::with_capacity(new_t)); }
        for row in b4.iter_mut() { row.grow(new_t); }
    };

    for g in gates {
        match *g {
            Gate::H(q) => {
                let prev = *wire[q].last().unwrap();
                let cur = next_var; next_var += 1; num_h += 1;
                wire[q].push(cur); grow(next_var, &mut b4, &mut v4);
                b4[prev].insert(cur); b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[a].last().unwrap();
                let vb = *wire[b].last().unwrap();
                grow(next_var, &mut b4, &mut v4);
                b4[va].insert(vb); b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[q].last().unwrap();
                grow(next_var, &mut b4, &mut v4); v4[v] = (v4[v] + 2) % 4;
            }
            Gate::S(q) => {
                let v = *wire[q].last().unwrap();
                grow(next_var, &mut b4, &mut v4); v4[v] = (v4[v] + 1) % 4;
            }
            Gate::T(q) => {
                let v = *wire[q].last().unwrap();
                rem.push(Z8Term { weight: 1, vars: vec![v] });
            }
        }
    }
    CompiledPhasePoly {
        num_qubits: n, num_vars: next_var, num_h, 
        output_vars: wire.iter().map(|w| *w.last().unwrap()).collect(),
        b4, v4, eps4: 0, rem,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

    fn approx_eq(a: Complex64, b: Complex64, eps: f64) {
        assert!((a - b).norm() <= eps, "a={:?} b={:?} |a-b|={}", a, b, (a - b).norm());
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
            assert!((amp - expected).norm() < 1e-10, "amp={:?} expected={:?}", amp, expected);
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
        let gates = vec![Gate::H(0), Gate::CZ(0, 1), Gate::T(1), Gate::H(2), Gate::CZ(1, 2), Gate::S(0), Gate::Z(2)];
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
    fn ry() {
        let gates = vec![Gate::H(0), Gate::T(0), Gate::H(0)];
        let poly = compile_clifford_t(1, &gates);
        let input = vec![0u8];
        let a0 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a1 = amplitude_clifford_t_accel(&poly, &input, 1);

        approx_eq(a0, c(0.8535533905932737, 0.35355339059327373), 1e-3);
        approx_eq(a1, c(0.1464466094067262, -0.35355339059327373), 1e-3);
    }

    #[test]
    fn benchmark_like_clifford_vs_nonclifford_amplitude() {
        use std::time::Instant;
        let n = 12;

        let mut g1 = Vec::new();
        for q in 0..n { g1.push(Gate::H(q)); }
        for q in 0..(n - 1) { g1.push(Gate::CZ(q, q + 1)); }
        for q in (0..n).rev() { g1.push(Gate::H(q)); }
        let p1 = compile_clifford_t(n, &g1);

        let start = Instant::now();
        let a1 = amplitude_clifford_t_accel(&p1, &vec![0u8; n], 0);
        let t1 = start.elapsed();
        println!("[TIMING] clifford-ish amp = {:?} in {:?}", a1, t1);

        let mut g2 = g1.clone();
        for q in 0..n {
            if q % 2 == 0 { g2.push(Gate::T(q)); }
            if q % 3 == 0 { g2.push(Gate::S(q)); }
        }
        let p2 = compile_clifford_t(n, &g2);

        let start2 = Instant::now();
        let a2 = amplitude_clifford_t_accel(&p2, &vec![0u8; n], 0);
        let t2 = start2.elapsed();
        println!("[TIMING] clifford+T/S amp = {:?} in {:?}", a2, t2);

        assert!(a1.re.is_finite() && a1.im.is_finite());
        assert!(a2.re.is_finite() && a2.im.is_finite());
    }
}