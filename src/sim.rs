use fixedbitset::FixedBitSet;
use num_complex::Complex64;

// the small QASM helper module lives at crate root; we re-export it here
// so that consumers of `sim` can still access `sim::Circuit` etc.
pub use crate::qc;

// convenience re-exports so callers can say `PolyQ::sim::Circuit` or import
// common helpers directly from this module without touching `qc`.
pub use crate::qc::Circuit;
pub use crate::qc::{
    read_qasm_file,
    write_qasm_file,
    write_qasm_string,
    QasmError,
};


/// Represents the quantum gates supported by the simulator.
/// In Montanaro's extended formalism, these are mapped to weights in Z8.
#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),         // Hadamard: Introduces a new variable and a weight-4 quadratic term.
    Z(usize),         // Pauli Z: Adds a weight-4 linear term.
    S(usize),         // Phase (S): Adds a weight-2 linear term.
    T(usize),         // T gate: Adds a weight-1 linear term.
    CZ(usize, usize), // Controlled-Z: Adds a weight-4 quadratic term.
}

/// Represents a non-Z2 phase term (from S or T gates) in the polynomial.
/// These terms cannot be solved directly via Dickson's theorem and are evaluated separately.
#[derive(Clone, Debug)]
pub struct Z8Term {
    pub weight: u8,       // Weight in Z8 (1 for T, 2 for S).
    pub vars: Vec<usize>, // The variables involved in this term.
}

/// The compiled representation of the quantum circuit as a Phase Polynomial.
/// The overall polynomial maps F2^t -> Z8.
#[derive(Clone, Debug)]
pub struct CompiledPhasePoly {
    pub num_qubits: usize,       // Initial number of variables (n).
    pub num_vars: usize,         // Total variables t = n + h (one new variable per H gate).
    pub num_h: usize,            // Total number of Hadamard gates.
    pub output_vars: Vec<usize>, // The final variable index residing on each qubit wire.
    pub b4: Vec<FixedBitSet>, // Symplectic adjacency matrix for weight-4 quadratic terms (CZ, H).
    pub v4: Vec<u8>,          // Vector for weight-4 linear terms (Z).
    pub eps4: u8,             // Global weight-4 phase/constant.
    pub rem: Vec<Z8Term>,     // Remainder terms (weights 1, 2) from T and S gates.
}

/// Compiles a sequence of quantum gates into a `CompiledPhasePoly`.
///
/// This iterates over the gates and applies rules to build the characteristic
/// boolean polynomial representing the circuit's action.
pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    // Tracks the current variable index on each qubit wire.
    let mut wire: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut next_var = n;
    let mut num_h = 0usize;

    // Initialize Clifford quadratic and linear structures.
    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];
    let eps4 = 0u8;

    let mut rem: Vec<Z8Term> = Vec::new();

    // Helper closure to dynamically expand the symplectic matrix and linear vector
    // whenever a Hadamard gate introduces a new variable.
    let mut grow_to = |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
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

    // Iterate through the circuit and apply polynomial updating rules.
    for g in gates {
        match *g {
            Gate::H(q) => {
                let prev = *wire[q].last().unwrap();
                let cur = next_var;
                next_var += 1;
                num_h += 1;

                wire[q].push(cur);
                grow_to(next_var, &mut b4, &mut v4);

                // Hadamard introduces a quadratic term x_prev * x_cur.
                b4[prev].insert(cur);
                b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[a].last().unwrap();
                let vb = *wire[b].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                // CZ introduces a quadratic term x_a * x_b.
                b4[va].insert(vb);
                b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[q].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                // Z introduces a linear term. Weight 4 in Z8 is equivalent to toggling bit in Z2.
                v4[v] ^= 1;
            }
            Gate::S(q) => {
                let v = *wire[q].last().unwrap();
                // S gate adds a linear term with weight 2 in Z8.
                rem.push(Z8Term {
                    weight: 2,
                    vars: vec![v],
                });
            }
            Gate::T(q) => {
                let v = *wire[q].last().unwrap();
                // T gate adds a linear term with weight 1 in Z8.
                rem.push(Z8Term {
                    weight: 1,
                    vars: vec![v],
                });
            }
        }
    }

    // Collect the final variable indices for each wire to use as outputs.
    let output_vars = wire.iter().map(|w| *w.last().unwrap()).collect::<Vec<_>>();

    CompiledPhasePoly {
        num_qubits: n,
        num_vars: next_var,
        num_h,
        output_vars,
        b4,
        v4,
        eps4,
        rem,
    }
}

/// Evaluates the non-Clifford (remainder) phase terms for a specific variable assignment.
/// Computes the sum of weights modulo 8.
#[inline(always)]
fn eval_rem_mod8(rem: &[Z8Term], mask: usize, vvars: &[usize], x_full: &[u8]) -> u8 {
    let mut acc = 0u8;

    for t in rem {
        let mut active = true;

        for &var in &t.vars {
            // determine the value of this variable: if it's a vvar use the mask,
            // otherwise look at the full assignment (covers fixed or uvars).
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

        if active {
            acc = (acc + t.weight) & 7;
        }
    }

    acc
}

/// Applies Dickson's Theorem to reduce a quadratic form to a canonical block-diagonal form.
/// This allows evaluating the exponential sum of Clifford circuits in O(m^3) time instead of O(2^m).
/// Returns the rank `r` of the symplectic matrix and the transformed linear vector.
fn dickson_reduce(mut b: Vec<FixedBitSet>, mut v: Vec<u8>) -> (usize, Vec<u8>) {
    let m = b.len();
    let mut r = 0usize;
    let mut p = 0usize;

    while p + 1 < m {
        // Find a pivot (a non-zero entry indicating a hyperbolic pair).
        let mut pivot: Option<(usize, usize)> = None;
        'outer: for i in p..m {
            for j in (i + 1)..m {
                if b[i].contains(j) {
                    pivot = Some((i, j));
                    break 'outer;
                }
            }
        }

        let Some((i, j)) = pivot else {
            break;
        };

        // Swap rows/cols to bring the pivot to the leading 2x2 block.
        b.swap(p, i);
        v.swap(p, i);

        let mut j2 = j;
        if j2 == p {
            j2 = i;
        }
        b.swap(p + 1, j2);
        v.swap(p + 1, j2);

        let rp = b[p].clone();
        let rp1 = b[p + 1].clone();

        // Perform congruent transformations to zero out the rest of the row/column.
        // This preserves the symplectic nature and the Hamming weight of the function.
        for k in (p + 2)..m {
            if b[k].contains(p) {
                b[k].symmetric_difference_with(&rp1);
                v[k] ^= v[p + 1];
            }
            if b[k].contains(p + 1) {
                b[k].symmetric_difference_with(&rp);
                v[k] ^= v[p];
            }
        }

        r += 2;
        p += 2;
    }

    (r, v)
}

/// Evaluates the exponential sum of a Z2 quadratic form using its Dickson-reduced rank.
fn z2_quadratic_exponential_sum(b: Vec<FixedBitSet>, v: Vec<u8>, eps: u8) -> f64 {
    let m = v.len();
    let (r, v2) = dickson_reduce(b, v);

    // If a kernel variable (outside the rank bounds) has a non-zero linear coefficient,
    // the function is balanced, and the exponential sum is exactly 0.
    for i in r..m {
        if v2[i] == 1 {
            return 0.0;
        }
    }

    // The magnitude depends solely on the number of variables and the matrix rank.
    let mag = (1u128 << (m - r / 2)) as f64;
    // Apply the global phase (epsilon).
    if eps == 1 { -mag } else { mag }
}

// previous helper removed: optimized amplitude implementation no longer needs a hashmap-based fixed set.

/// Calculates the transition amplitude <y | U | input> for the compiled circuit.
/// This function partitions variables into Clifford and non-Clifford sets to accelerate computation.
pub fn amplitude_clifford_t_accel(
    poly: &CompiledPhasePoly,
    input: &[u8],
    target_y: usize,
) -> Complex64 {
    let n = poly.num_qubits;
    let t = poly.num_vars;

    let mut fixed = vec![None; t];

    for i in 0..n {
        fixed[i] = Some(input[i]);
    }

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
        for &v in &term.vars {
            in_rem[v] = true;
        }
    }

    let mut internal = Vec::new();
    for i in 0..t {
        if fixed[i].is_none() {
            internal.push(i);
        }
    }

    let mut vvars = Vec::new();
    let mut uvars = Vec::new();

    for &v in &internal {
        if in_rem[v] {
            vvars.push(v);
        } else {
            uvars.push(v);
        }
    }

    let nv = vvars.len();
    let nu = uvars.len();

    let mut x_full = vec![0u8; t];

    for i in 0..t {
        if let Some(v) = fixed[i] {
            x_full[i] = v;
        }
    }

    let mut eps_base = poly.eps4 & 1;

    for i in 0..t {
        if let Some(v) = fixed[i] {
            if v == 1 && (poly.v4[i] & 1) == 1 {
                eps_base ^= 1;
            }
        }
    }

    // include quadratic contributions from pairs of fixed variables
    let fixed_list: Vec<usize> = fixed
        .iter()
        .enumerate()
        .filter_map(|(idx, &v)| if v == Some(1) { Some(idx) } else { None })
        .collect();
    for a in 0..fixed_list.len() {
        for b in (a + 1)..fixed_list.len() {
            let i = fixed_list[a];
            let j = fixed_list[b];
            if poly.b4[i].contains(j) {
                eps_base ^= 1;
            }
        }
    }

    let mut bu: Vec<FixedBitSet> =
        (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();

    let mut vu_base = vec![0u8; nu];

    for (ui, &orig_u) in uvars.iter().enumerate() {
        vu_base[ui] = poly.v4[orig_u] & 1;
    }

    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (uj, &orig_uj) in uvars.iter().enumerate() {
            if ui < uj && poly.b4[orig_u].contains(orig_uj) {
                bu[ui].insert(uj);
                bu[uj].insert(ui);
            }
        }

        for f in 0..t {
            if let Some(val) = fixed[f] {
                if val == 1 && poly.b4[orig_u].contains(f) {
                    vu_base[ui] ^= 1;
                }
            }
        }
    }

    let mut cross: Vec<Vec<usize>> = vec![Vec::new(); nu];

    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (vj, &orig_v) in vvars.iter().enumerate() {
            if poly.b4[orig_u].contains(orig_v) {
                cross[ui].push(vj);
            }
        }
    }

    let total_v = 1usize << nv;

    let mut amp = Complex64::new(0.0, 0.0);

    let mut vu = vec![0u8; nu];

    for mask in 0..total_v {
        for (j, &orig_v) in vvars.iter().enumerate() {
            x_full[orig_v] = ((mask >> j) & 1) as u8;
        }

        let r = eval_rem_mod8(&poly.rem, mask, &vvars, &x_full);

        let theta = std::f64::consts::PI * (r as f64) / 4.0;

        let phase = Complex64::from_polar(1.0, theta);

        vu.copy_from_slice(&vu_base);

        for ui in 0..nu {
            let mut togg = 0u8;

            for &vj in &cross[ui] {
                togg ^= ((mask >> vj) & 1) as u8;
            }

            vu[ui] ^= togg;
        }

        // compute eps including contributions from mask bits and fixed interactions
        let mut eps = eps_base;
        // contributions from vvars linear terms
        for (j, &orig_v) in vvars.iter().enumerate() {
            if ((mask >> j) & 1) == 1 && (poly.v4[orig_v] & 1) == 1 {
                eps ^= 1;
            }
        }
        // fixed-vvar quadratic interactions
        for i in 0..t {
            if let Some(fval) = fixed[i] {
                if fval == 1 {
                    for (j, &orig_v) in vvars.iter().enumerate() {
                        if ((mask >> j) & 1) == 1 && poly.b4[i].contains(orig_v) {
                            eps ^= 1;
                        }
                    }
                }
            }
        }
        // vvar-vvar quadratic interactions
        for a in 0..nv {
            if ((mask >> a) & 1) == 0 {
                continue;
            }
            for b in (a + 1)..nv {
                if ((mask >> b) & 1) == 0 {
                    continue;
                }
                if poly.b4[vvars[a]].contains(vvars[b]) {
                    eps ^= 1;
                }
            }
        }

        let inner = z2_quadratic_exponential_sum(bu.clone(), vu.clone(), eps);

        amp += phase * inner;
    }

    let norm = (2f64).powf(-(poly.num_h as f64) / 2.0);

    amp * norm
}

/*
Experimental
Start
*/

pub fn simulate_statevector(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    // This is a placeholder for a function that would simulate the entire state vector of the circuit.
    // It would iterate over all 2^n possible output states and compute their amplitudes using the accelerated method.
    // unimplemented!()
    // let mut outvec = Vec::new();
    let size = 1usize << poly.num_qubits;
    let mut outvec = vec![Complex64::new(0.0, 0.0); size];

    for i in 0..(1usize << poly.num_qubits) {
        let amp = amplitude_clifford_t_accel(poly, input, i);
        // println!("y={:0width$b} amp={:?}", i, amp, width=poly.num_qubits); //ddebug
        outvec[i] = amp;
    }
    // unimplemented!()
    outvec
}

#[test]
fn test_simulate_statevector() {
    // Test case for the state vector simulation function
    let gates = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
    let poly = compile_clifford_t(2, &gates);
    let input = vec![0u8, 0u8];
    let statevec = simulate_statevector(&poly, &input);
    let expected_output = vec![
        Complex64::new(0.5f64.sqrt(), 0.0), // |00>
        Complex64::new(0.0, 0.0),           // |01>
        Complex64::new(0.0, 0.0),           // |10>
        Complex64::new(0.5f64.sqrt(), 0.0), // |11>
    ];

    // println!("State vector: {:?}", statevec);
    for (amp, expected) in statevec.iter().zip(expected_output.iter()) {
        assert!(
            (amp - expected).norm() < 1e-10,
            "amp={:?} expected={:?}",
            amp,
            expected
        );
    }
}

//Hybrid parallelization

pub fn amplitude_clifford_t_multithreaded(
    poly: &CompiledPhasePoly,
    input: &[u8],
    target_y: usize,
) -> Complex64 {
    // This is a placeholder for the non-accelerated version of the amplitude calculation.
    // It would simply iterate over all 2^t variable assignments and evaluate the phase polynomial directly.
    // This is exponentially slow in t and is only intended for testing correctness against the accelerated version.
    unimplemented!()
}

pub fn simulate_statevector_clifford_t_multiprocess(
    poly: &CompiledPhasePoly,
    input: &[u8],
) -> Vec<Complex64> {
    // This is a placeholder for a multithreaded version of the state vector simulation.
    // It would use a parallel iterator (e.g., from the Rayon crate) to compute amplitudes for all output states in parallel.
    unimplemented!()
}

/*
Experimental
End
*/

/*
WIP
Start
*/

pub fn optimised_transpile(gates: &[Gate], t_count: usize) -> CompiledPhasePoly {
    // determine how many qubits the circuit acts on
    let num_qubits = gates.iter().fold(0usize, |acc, g| match *g {
        Gate::H(q) | Gate::Z(q) | Gate::S(q) | Gate::T(q) => acc.max(q + 1),
        Gate::CZ(a, b) => acc.max(a + 1).max(b + 1),
    });

    // ------------------------------------------------------------------
    // phase‑1: circuit‑level optimisation
    // ------------------------------------------------------------------
    // currently a no‑op; a real implementation would e.g. commute/cancel
    // Paulis, resynthesise Clifford subcircuits, reduce the T‑count, …
    let mut opt_gates = gates.to_vec();

    // “use” t_count so we don’t get an unused‑variable warning
    if opt_gates.len() > (t_count as usize) {
        // placeholder branch – in the future we would try to rewrite the
        // circuit so that the number of non‑Clifford gates ≤ t_count.
    }

    // ------------------------------------------------------------------
    // phase‑2: compile to a phase polynomial
    // ------------------------------------------------------------------
    let mut poly = compile_clifford_t(num_qubits, &opt_gates);

    // ------------------------------------------------------------------
    // phase‑3: polynomial‑level optimisation / approximate compilation
    // ------------------------------------------------------------------
    // again, this is just a stub.  A true optimiser might eliminate
    // zero rows/columns of `b4`, merge linear terms, fold fixed
    // variables, or even return an “approximate” polynomial if the
    // variable budget is exceeded.
    if poly.num_vars > (t_count as usize) {
        // record the fact that the budget was exceeded; leave `poly`
        // unchanged for now.
    }

    unimplemented!()
}

/*
WIP
End
*/

#[cfg(test)]
mod tests {
    // rayon::vec import removed; not used in tests

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
        let gates = vec![
            Gate::H(0),
            Gate::S(0),
            Gate::H(1),
            Gate::CZ(0, 1),
            Gate::H(1),
        ];
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
        let gates = vec![
            Gate::H(0),
            Gate::T(0),
            Gate::H(1),
            Gate::CZ(0, 1),
            Gate::H(1),
        ];
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
    fn ry() {
        let gates = vec![Gate::H(0), Gate::T(0), Gate::H(0)];
        let poly = compile_clifford_t(1, &gates);
        let input = vec![0u8];
        let a0 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a1 = amplitude_clifford_t_accel(&poly, &input, 1);

        //composere result [ 0.854+0.354j, 0.146-0.354j ]
        println!("a0={:?} a1={:?}", a0, a1);
        approx_eq(a0, c(0.8535533905932737, 0.35355339059327373), 1e-3);
        approx_eq(a1, c(0.1464466094067262, -0.35355339059327373), 1e-3);
    }

    #[test]
    fn benchmark_like_clifford_vs_nonclifford_amplitude() {
        use std::time::Instant;

        let n = 12;

        let mut g1 = Vec::new();
        for q in 0..n {
            g1.push(Gate::H(q));
        }
        for q in 0..(n - 1) {
            g1.push(Gate::CZ(q, q + 1));
        }
        for q in (0..n).rev() {
            g1.push(Gate::H(q));
        }
        let p1 = compile_clifford_t(n, &g1);

        let start = Instant::now();
        let a1 = amplitude_clifford_t_accel(&p1, &vec![0u8; n], 0);
        let t1 = start.elapsed();
        println!("[TIMING] clifford-ish amp = {:?} in {:?}", a1, t1);

        let mut g2 = g1.clone();
        for q in 0..n {
            if q % 2 == 0 {
                g2.push(Gate::T(q));
            }
            if q % 3 == 0 {
                g2.push(Gate::S(q));
            }
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
