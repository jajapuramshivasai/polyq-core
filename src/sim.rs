use fixedbitset::FixedBitSet;
use num_complex::Complex64;

pub use crate::qc;
pub use crate::qc::Circuit;
pub use crate::qc::{
    read_qasm_file,
    write_qasm_file,
    write_qasm_string,
    QasmError,
};

/// Represents the quantum gates supported by the simulator.
#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),         // Hadamard: Introduces a new variable and a weight-2 quadratic term in Z4.
    Z(usize),         // Pauli Z: Adds a weight-2 linear term in Z4.
    S(usize),         // Phase (S): Adds a weight-1 linear term in Z4 (NOW EVALUATED IN POLYNOMIAL TIME).
    T(usize),         // T gate: Adds a weight-1 linear term in Z8 (Pushed to remainder).
    CZ(usize, usize), // Controlled-Z: Adds a weight-2 quadratic term in Z4.
}

/// Represents a non-Clifford phase term (from T gates) in the polynomial.
#[derive(Clone, Debug)]
pub struct Z8Term {
    pub weight: u8,       
    pub vars: Vec<usize>, 
}

/// The compiled representation of the quantum circuit as a Phase Polynomial.
#[derive(Clone, Debug)]
pub struct CompiledPhasePoly {
    pub num_qubits: usize,       
    pub num_vars: usize,         
    pub num_h: usize,            
    pub output_vars: Vec<usize>, 
    pub b4: Vec<FixedBitSet>, // Symplectic adjacency matrix for Z4 quadratic terms (weight 2).
    pub v4: Vec<u8>,          // Z4 linear vector (values 0, 1, 2, 3).
    pub eps4: u8,             // Global Z4 phase offset.
    pub rem: Vec<Z8Term>,     // Remainder terms (Z8) exclusively for T gates.
}

/// Compiles a sequence of quantum gates into a `CompiledPhasePoly`.
pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    let mut wire: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut next_var = n;
    let mut num_h = 0usize;

    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];
    let eps4 = 0u8;

    let mut rem: Vec<Z8Term> = Vec::new();

    let mut grow_to = |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
        while v4.len() < new_t { v4.push(0); }
        while b4.len() < new_t { b4.push(FixedBitSet::with_capacity(new_t)); }
        for row in b4.iter_mut() { row.grow(new_t); }
    };

    for g in gates {
        match *g {
            Gate::H(q) => {
                let prev = *wire[q].last().unwrap();
                let cur = next_var;
                next_var += 1;
                num_h += 1;

                wire[q].push(cur);
                grow_to(next_var, &mut b4, &mut v4);

                b4[prev].insert(cur);
                b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[a].last().unwrap();
                let vb = *wire[b].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                
                b4[va].insert(vb);
                b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[q].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                // Z is a weight-2 linear term in Z4.
                v4[v] = (v4[v] + 2) % 4;
            }
            Gate::S(q) => {
                let v = *wire[q].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                // S is a weight-1 linear term in Z4. 
                // Absorbed natively! No longer pushed to `rem`.
                v4[v] = (v4[v] + 1) % 4;
            }
            Gate::T(q) => {
                let v = *wire[q].last().unwrap();
                // T remains a Z8 weight-1 term.
                rem.push(Z8Term { weight: 1, vars: vec![v] });
            }
        }
    }

    let output_vars = wire.iter().map(|w| *w.last().unwrap()).collect::<Vec<_>>();

    CompiledPhasePoly {
        num_qubits: n, num_vars: next_var, num_h, output_vars, b4, v4, eps4, rem,
    }
}

#[inline(always)]
fn eval_rem_mod8(rem: &[Z8Term], mask: usize, vvars: &[usize], x_full: &[u8]) -> u8 {
    let mut acc = 0u8;
    for t in rem {
        let mut active = true;
        for &var in &t.vars {
            let val = if let Some(p) = vvars.iter().position(|&x| x == var) {
                ((mask >> p) & 1) as u8
            } else {
                x_full[var]
            };
            if val == 0 {
                active = false;
                break;
            }
        }
        if active { acc = (acc + t.weight) & 7; }
    }
    acc
}

/// Applies a Z4-aware Dickson's Theorem to reduce the symplectic matrix.
fn dickson_reduce_z4(mut b: Vec<FixedBitSet>, mut v: Vec<u8>) -> (usize, Vec<u8>, Vec<FixedBitSet>) {
    let m = b.len();
    let mut r = 0usize;
    let mut p = 0usize;

    while p + 1 < m {
        let mut pivot: Option<(usize, usize)> = None;
        'outer: for i in p..m {
            for j in (i + 1)..m {
                if b[i].contains(j) {
                    pivot = Some((i, j));
                    break 'outer;
                }
            }
        }

        let Some((i, j)) = pivot else { break; };

        b.swap(p, i);
        v.swap(p, i);

        let mut j2 = j;
        if j2 == p { j2 = i; }
        b.swap(p + 1, j2);
        v.swap(p + 1, j2);

        let rp = b[p].clone();
        let rp1 = b[p + 1].clone();

        for k in (p + 2)..m {
            // Substitution: x_k -> x_k ^ x_{p+1}
            if b[k].contains(p) {
                b[k].symmetric_difference_with(&rp1);
                v[p + 1] = (v[p + 1] + v[k]) % 4; // Z4 Linear Update
                if v[k] % 2 != 0 {                // Z4 Quadratic Update
                    b[k].toggle(p + 1);
                    b[p + 1].toggle(k);
                }
            }
            // Substitution: x_k -> x_k ^ x_p
            if b[k].contains(p + 1) {
                b[k].symmetric_difference_with(&rp);
                v[p] = (v[p] + v[k]) % 4;         // Z4 Linear Update
                if v[k] % 2 != 0 {                // Z4 Quadratic Update
                    b[k].toggle(p);
                    b[p].toggle(k);
                }
            }
        }
        r += 2;
        p += 2;
    }
    (r, v, b)
}

/// Evaluates the exponential sum of a Z4 quadratic form.
fn z4_quadratic_exponential_sum(b: Vec<FixedBitSet>, v: Vec<u8>) -> Complex64 {
    let m = b.len();
    let (r, v2, b2) = dickson_reduce_z4(b, v);
    let mut sum_val = Complex64::new(1.0, 0.0);

    // Evaluate the isolated 2x2 hyperbolic blocks
    let mut p = 0;
    while p < r {
        let mut pair_sum = Complex64::new(0.0, 0.0);
        for x_p in 0..=1 {
            for x_p1 in 0..=1 {
                let quad = if b2[p].contains(p + 1) { 2 * x_p * x_p1 } else { 0 };
                let lin = v2[p] * x_p + v2[p + 1] * x_p1;
                let phase_z4 = (quad + lin) % 4;
                
                let phase_c = match phase_z4 {
                    0 => Complex64::new(1.0, 0.0),
                    1 => Complex64::new(0.0, 1.0),
                    2 => Complex64::new(-1.0, 0.0),
                    3 => Complex64::new(0.0, -1.0),
                    _ => unreachable!(),
                };
                pair_sum += phase_c;
            }
        }
        sum_val *= pair_sum;
        p += 2;
    }

    // Evaluate the isolated kernel variables
    for k in r..m {
        let term = match v2[k] % 4 {
            0 => Complex64::new(2.0, 0.0),
            1 => Complex64::new(1.0, 1.0),
            2 => Complex64::new(0.0, 0.0), // Destructive interference
            3 => Complex64::new(1.0, -1.0),
            _ => unreachable!(),
        };
        if term == Complex64::new(0.0, 0.0) {
            return Complex64::new(0.0, 0.0);
        }
        sum_val *= term;
    }

    sum_val
}

pub fn amplitude_clifford_t_accel(
    poly: &CompiledPhasePoly,
    input: &[u8],
    target_y: usize,
) -> Complex64 {
    let n = poly.num_qubits;
    let t = poly.num_vars;

    let mut fixed = vec![None; t];
    for i in 0..n { fixed[i] = Some(input[i]); }
    for i in 0..n {
        let ov = poly.output_vars[i];
        let bit = ((target_y >> i) & 1) as u8;
        match fixed[ov] {
            None => fixed[ov] = Some(bit),
            Some(v) if v != bit => return Complex64::new(0.0, 0.0),
            _ => {}
        }
    }

    let mut in_rem = vec![false; t];
    for term in &poly.rem {
        for &v in &term.vars { in_rem[v] = true; }
    }

    let mut internal = Vec::new();
    for i in 0..t {
        if fixed[i].is_none() { internal.push(i); }
    }

    let mut vvars = Vec::new();
    let mut uvars = Vec::new();
    for &v in &internal {
        if in_rem[v] { vvars.push(v); } else { uvars.push(v); }
    }

    let nv = vvars.len();
    let nu = uvars.len();
    let mut x_full = vec![0u8; t];
    for i in 0..t {
        if let Some(v) = fixed[i] { x_full[i] = v; }
    }

    let mut eps_base_z4 = poly.eps4 % 4;
    for i in 0..t {
        if let Some(1) = fixed[i] {
            eps_base_z4 = (eps_base_z4 + poly.v4[i]) % 4;
        }
    }

    let fixed_list: Vec<usize> = fixed.iter().enumerate()
        .filter_map(|(idx, &v)| if v == Some(1) { Some(idx) } else { None })
        .collect();
    for a in 0..fixed_list.len() {
        for b in (a + 1)..fixed_list.len() {
            let i = fixed_list[a];
            let j = fixed_list[b];
            if poly.b4[i].contains(j) {
                eps_base_z4 = (eps_base_z4 + 2) % 4;
            }
        }
    }

    let mut bu: Vec<FixedBitSet> = (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();
    let mut vu_base = vec![0u8; nu];

    for (ui, &orig_u) in uvars.iter().enumerate() { vu_base[ui] = poly.v4[orig_u]; }

    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (uj, &orig_uj) in uvars.iter().enumerate() {
            if ui < uj && poly.b4[orig_u].contains(orig_uj) {
                bu[ui].insert(uj);
                bu[uj].insert(ui);
            }
        }
        for f in &fixed_list {
            if poly.b4[orig_u].contains(*f) {
                vu_base[ui] = (vu_base[ui] + 2) % 4;
            }
        }
    }

    let mut cross: Vec<Vec<usize>> = vec![Vec::new(); nu];
    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (vj, &orig_v) in vvars.iter().enumerate() {
            if poly.b4[orig_u].contains(orig_v) { cross[ui].push(vj); }
        }
    }

    let total_v = 1usize << nv;
    let mut amp = Complex64::new(0.0, 0.0);
    let mut vu = vec![0u8; nu];

    for mask in 0..total_v {
        for (j, &orig_v) in vvars.iter().enumerate() {
            x_full[orig_v] = ((mask >> j) & 1) as u8;
        }

        let r_z8 = eval_rem_mod8(&poly.rem, mask, &vvars, &x_full);

        vu.copy_from_slice(&vu_base);
        for ui in 0..nu {
            let mut togg = 0u8;
            for &vj in &cross[ui] { togg ^= ((mask >> vj) & 1) as u8; }
            if togg == 1 { vu[ui] = (vu[ui] + 2) % 4; }
        }

        let mut eps_z4 = eps_base_z4;
        for (j, &orig_v) in vvars.iter().enumerate() {
            if ((mask >> j) & 1) == 1 {
                eps_z4 = (eps_z4 + poly.v4[orig_v]) % 4;
            }
        }
        for &f in &fixed_list {
            for (j, &orig_v) in vvars.iter().enumerate() {
                if ((mask >> j) & 1) == 1 && poly.b4[f].contains(orig_v) {
                    eps_z4 = (eps_z4 + 2) % 4;
                }
            }
        }
        for a in 0..nv {
            if ((mask >> a) & 1) == 0 { continue; }
            for b in (a + 1)..nv {
                if ((mask >> b) & 1) == 0 { continue; }
                if poly.b4[vvars[a]].contains(vvars[b]) {
                    eps_z4 = (eps_z4 + 2) % 4;
                }
            }
        }

        // Combine the Z8 remainder phase with the Z4 base phase mapping
        let total_phase_z8 = (r_z8 + (eps_z4 * 2)) % 8;
        let phase_c = Complex64::from_polar(1.0, std::f64::consts::PI * (total_phase_z8 as f64) / 4.0);

        let inner = z4_quadratic_exponential_sum(bu.clone(), vu.clone());
        amp += phase_c * inner;
    }

    let norm = (2f64).powf(-(poly.num_h as f64) / 2.0);
    amp * norm
}

pub fn simulate_statevector(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    let size = 1usize << poly.num_qubits;
    let mut outvec = vec![Complex64::new(0.0, 0.0); size];
    for i in 0..(1usize << poly.num_qubits) {
        outvec[i] = amplitude_clifford_t_accel(poly, input, i);
    }
    outvec
}

pub fn amplitude_clifford_t_Cached(_poly: &CompiledPhasePoly, _input: &[u8], _target_y: usize) -> Complex64 {
    unimplemented!()
}

pub fn amplitude_clifford_t_multithreaded(_poly: &CompiledPhasePoly, _input: &[u8], _target_y: usize) -> Complex64 {
    unimplemented!()
}

pub fn simulate_statevector_clifford_t_multiprocess(_poly: &CompiledPhasePoly, _input: &[u8]) -> Vec<Complex64> {
    unimplemented!()
}

pub fn optimised_transpile(gates: &[Gate], t_count: usize) -> CompiledPhasePoly {
    let num_qubits = gates.iter().fold(0usize, |acc, g| match *g {
        Gate::H(q) | Gate::Z(q) | Gate::S(q) | Gate::T(q) => acc.max(q + 1),
        Gate::CZ(a, b) => acc.max(a + 1).max(b + 1),
    });

    let opt_gates = gates.to_vec();
    if opt_gates.len() > (t_count as usize) {}

    let poly = compile_clifford_t(num_qubits, &opt_gates);

    if poly.num_vars > (t_count as usize) {}

    unimplemented!()
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