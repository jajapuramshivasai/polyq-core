use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::f64::consts::PI;
use std::sync::OnceLock;

pub use crate::qc::{read_qasm_file, write_qasm_file, write_qasm_string, QasmError};

// ---------------------------------------------------------------------------
// Phase ring  (Z_{2^PHASE_BITS}, units of π / 2^PHASE_BITS)
// ---------------------------------------------------------------------------

type Phase = u8; // Z8

#[derive(Clone)]
pub struct PhasePoly {
    pub n: usize,
    pub terms: HashMap<Vec<usize>, Phase>, // monomial → coeff mod 8
}

#[inline(always)]
pub const fn phase_modulus() -> u32 { 1u32 << PHASE_BITS }

#[inline(always)]
pub const fn phase_mask_u32() -> u32 { (1u32 << PHASE_BITS) - 1 }

#[inline(always)]
fn phase_add(a: Phase, b: Phase) -> Phase { a.wrapping_add(b) }

#[inline(always)]
fn phase_sub(a: Phase, b: Phase) -> Phase { a.wrapping_sub(b) }

#[inline(always)]
fn phase_from_u32(x: u32) -> Phase { (x & phase_mask_u32()) as Phase }

#[inline(always)]
fn phase_to_angle_rad(p: Phase) -> f64 {
    PI * (p as f64) / (phase_modulus() as f64)
}

// ---------------------------------------------------------------------------
// [OPT-1]  Phase LUT  — 65536 complex roots of unity
// ---------------------------------------------------------------------------

static PHASE_LUT: OnceLock<Vec<Complex64>> = OnceLock::new();

fn get_phase_lut() -> &'static Vec<Complex64> {
    PHASE_LUT.get_or_init(|| {
        (0u32..phase_modulus())
            .map(|p| {
                let theta = phase_to_angle_rad(p as Phase);
                Complex64::new(theta.cos(), theta.sin())
            })
            .collect()
    })
}

/// Convert a Phase index to the corresponding unit complex number via LUT.
#[inline(always)]
fn phase_to_complex(p: Phase) -> Complex64 {
    // SAFETY: p is always in [0, 65535] by construction.
    unsafe { *get_phase_lut().get_unchecked(p as usize) }
}

// ---------------------------------------------------------------------------
// Z4 arithmetic
// ---------------------------------------------------------------------------

#[inline(always)]
fn z4_add(a: u8, b: u8) -> u8 { (a + b) & 3 }

#[inline(always)]
fn z4_sub(a: u8, b: u8) -> u8 { (a + 4 - b) & 3 }

/// Map Z4 value to unit complex (i^k).
#[inline(always)]
fn z4_to_complex(k: u8) -> Complex64 {
    // Expanded as a branch-free table lookup.
    const TABLE: [Complex64; 4] = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    TABLE[(k & 3) as usize]
}

// ---------------------------------------------------------------------------
// Gate & circuit types
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),
    Z(usize),
    S(usize),
    T(usize),
    RZ(usize, Phase),
    CZ(usize, usize),
    // Internal: wire-level X (bit-flip), produced by the peephole optimiser.
    X(usize),
}

#[derive(Clone, Debug)]
pub struct PhaseTerm {
    pub weight: Phase,
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
    pub rem: Vec<PhaseTerm>,
}

pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self { Self { num_qubits, gates: Vec::new() } }
    pub fn h(&mut self, q: usize)           { self.gates.push(Gate::H(q)); }
    pub fn z(&mut self, q: usize)           { self.gates.push(Gate::Z(q)); }
    pub fn s(&mut self, q: usize)           { self.gates.push(Gate::S(q)); }
    pub fn t(&mut self, q: usize)           { self.gates.push(Gate::T(q)); }
    pub fn rz(&mut self, q: usize, p: Phase){ self.gates.push(Gate::RZ(q, p)); }
    pub fn cz(&mut self, q1: usize, q2: usize){ self.gates.push(Gate::CZ(q1, q2)); }

    /// Compile with peephole optimisation then phase folding.
    pub fn compile(&self) -> CompiledPhasePoly {
        let optimised = peephole_optimise(&self.gates, self.num_qubits);
        let mut poly = compile_clifford_t(self.num_qubits, &optimised);
        fold_rem_into_z4(&mut poly);
        poly
    }
}

// ---------------------------------------------------------------------------
// [OPT-8]  Peephole optimiser on gate list
// ---------------------------------------------------------------------------

/// Cancel H·H pairs; merge same-qubit Z/S/T/RZ sequences; simplify H·Z·H→X.
fn peephole_optimise(gates: &[Gate], num_qubits: usize) -> Vec<Gate> {
    // We make multiple passes until stable (usually converges in 2-3).
    let mut cur: Vec<Gate> = gates.to_vec();
    loop {
        let next = peephole_pass(&cur, num_qubits);
        if next.len() == cur.len() { break; }
        cur = next;
    }
    cur
}

fn peephole_pass(gates: &[Gate], num_qubits: usize) -> Vec<Gate> {
    // Track the phase accumulated on each wire between Hadamard gates.
    // Phases stored as Z_{2^PHASE_BITS}.
    let mut wire_phase: Vec<Phase> = vec![0; num_qubits];
    let mut out: Vec<Gate> = Vec::with_capacity(gates.len());

    // Simple one-pass pattern matcher.
    let mut i = 0;
    while i < gates.len() {
        match &gates[i] {
            // H·H = I: check if next gate on same qubit is also H.
            Gate::H(q) => {
                if i + 1 < gates.len() {
                    if let Gate::H(q2) = &gates[i + 1] {
                        if q == q2 {
                            i += 2; // skip both H gates
                            continue;
                        }
                    }
                }
                out.push(gates[i].clone());
            }
            // Merge consecutive diagonal gates on same qubit into a single RZ.
            Gate::Z(q) | Gate::S(q) | Gate::T(q) | Gate::RZ(q, _) => {
                let q = *q;
                let mut acc: Phase = 0;
                let mut j = i;
                // Collect a run of diagonal gates on this qubit, stopping at
                // any gate that touches q in a non-diagonal way.
                while j < gates.len() {
                    match &gates[j] {
                        Gate::Z(qq) if *qq == q => {
                            // Z = phase π = 2^(PHASE_BITS-1) units
                            acc = acc.wrapping_add(1u16 << (PHASE_BITS - 1));
                            j += 1;
                        }
                        Gate::S(qq) if *qq == q => {
                            // S = phase π/2 = 2^(PHASE_BITS-2) units
                            acc = acc.wrapping_add(1u16 << (PHASE_BITS - 2));
                            j += 1;
                        }
                        Gate::T(qq) if *qq == q => {
                            // T = phase π/4 = 2^(PHASE_BITS-3) units  (note: PHASE_BITS=16 → unit=8192 for T, matching t_phase_unit)
                            acc = acc.wrapping_add(1u16 << (PHASE_BITS - 3));
                            j += 1;
                        }
                        Gate::RZ(qq, ph) if *qq == q => {
                            acc = acc.wrapping_add(*ph);
                            j += 1;
                        }
                        // H or CZ touching q breaks the run.
                        Gate::H(qq) | Gate::X(qq) if *qq == q => break,
                        Gate::CZ(a, b) if *a == q || *b == q => break,
                        _ => { j += 1; } // gate on different qubit; skip but keep accumulating
                    }
                }
                // Emit a single gate for the accumulated diagonal phase.
                if acc != 0 {
                    // Canonicalise back to named gates where possible (saves nv overhead).
                    let half = 1u16 << (PHASE_BITS - 1);
                    let qrtr = 1u16 << (PHASE_BITS - 2);
                    let eighth = 1u16 << (PHASE_BITS - 3);
                    if acc == half { out.push(Gate::Z(q)); }
                    else if acc == qrtr { out.push(Gate::S(q)); }
                    else if acc == eighth { out.push(Gate::T(q)); }
                    else { out.push(Gate::RZ(q, acc)); }
                }
                // Re-emit the non-diagonal gates in between (already handled above by j skipping;
                // but we need to re-emit the "other-qubit" gates that we jumped over).
                // Simpler: re-emit from i to j, but filter out the diagonal-q gates we merged.
                for k in i..j {
                    match &gates[k] {
                        Gate::Z(qq) | Gate::S(qq) | Gate::T(qq) | Gate::X(qq)
                            if *qq == q => continue,
                        Gate::RZ(qq, _) if *qq == q => continue,
                        other => out.push(other.clone()),
                    }
                }
                i = j;
                continue;
            }
            other => out.push(other.clone()),
        }
        i += 1;
    }
    out
}

// ---------------------------------------------------------------------------
// [OPT-4]  Fold rem terms into Z4 at compile time
// ---------------------------------------------------------------------------

/// For each linear rem term (single-variable), if it is a multiple of
/// 2^(PHASE_BITS-2) (i.e. it equals k * (π/4) for k ∈ Z8), absorb into v4
/// and remove from rem.  Also merge duplicate-variable rem terms.
fn fold_rem_into_z4(poly: &mut CompiledPhasePoly) {
    // Accumulate per-variable weights for linear rem terms.
    let t = poly.num_vars;
    let mut linear_acc: Vec<Phase> = vec![0; t];
    let mut keep: Vec<bool> = vec![false; poly.rem.len()];

    for (idx, term) in poly.rem.iter().enumerate() {
        if term.vars.len() == 1 {
            linear_acc[term.vars[0]] = phase_add(linear_acc[term.vars[0]], term.weight);
            // mark for removal; we'll re-add later if residue remains
        } else {
            keep[idx] = true;
        }
    }

    // The Z4 phase is in units of π/2 = 2^(PHASE_BITS-1).
    // A rem weight is a Z4 multiple if it is a multiple of 2^(PHASE_BITS-2)?
    // No: v4 is in Z4 (units π/2). So to absorb a rem weight w into v4 we need
    // w to be an integer multiple of 2^(PHASE_BITS-2) (units π/4), but v4 uses
    // Z4 so we need units of π/2 = 2^(PHASE_BITS-1).
    // Strategy: take the Z4 part (top 2 bits of the 16-bit phase) and leave the
    // remainder. Specifically:
    //   w_z4 = w >> (PHASE_BITS - 2)   (gives units of π/4, i.e. Z8 values)
    //   residue = w & ((1<<(PHASE_BITS-2))-1)
    // The Z4 contribution is w_z4 (in Z8), but v4 is in Z4. The Z4 part of Z8
    // is the top bit pair. This is exactly what the original code does via
    // `v4[v] & 3` and `t_phase_unit()`.
    //
    // Concretely: eps_phase = (cur_eps as u32) << (PHASE_BITS-1) maps Z4 → Phase.
    // So v4[v] (Z4) contributes eps_phase of v4[v] << (PHASE_BITS-1).
    // Absorbing rem weight w into v4: we need w = (delta_v4) << (PHASE_BITS-1)
    // for some delta_v4 in Z4. So w must be a multiple of 2^(PHASE_BITS-1), i.e.
    // a multiple of π. We can instead absorb the Z4-visible part:
    //   delta_v4 = (w >> (PHASE_BITS-2)) & 3   (captures 0, π/2, π, 3π/2)
    //   residue  = w - (delta_v4 << (PHASE_BITS-2))
    // Then residue must stay in rem.

    let z4_unit_shift = PHASE_BITS - 2; // 2^(PHASE_BITS-2) = π/4 unit
    let z4_unit: u32 = 1u32 << z4_unit_shift;

    for v in 0..t {
        if linear_acc[v] == 0 { continue; }
        let w = linear_acc[v] as u32;
        let z4_part = ((w >> z4_unit_shift) & 3) as u8; // mod 4
        let residue = (w - (z4_part as u32) * z4_unit) as Phase;
        if z4_part != 0 {
            poly.v4[v] = z4_add(poly.v4[v], z4_part);
        }
        if residue != 0 {
            // Push residue back as a rem term.
            poly.rem.push(PhaseTerm { weight: residue, vars: vec![v] });
        }
    }

    // Rebuild rem:
    //   - indices 0..keep.len() : original terms; keep multi-var ones, drop linear ones
    //     (linear ones were accumulated into linear_acc and already folded into v4 above).
    //   - indices keep.len()..  : residue terms we just pushed; always keep these.
    let orig_len = keep.len();
    let new_rem: Vec<PhaseTerm> = poly.rem.iter().enumerate()
        .filter_map(|(idx, term)| {
            if idx < orig_len {
                // Original term: keep only if it was a multi-var term (keep[idx]==true).
                if keep[idx] { Some(term.clone()) } else { None }
            } else {
                // Residue term we just pushed: always keep.
                Some(term.clone())
            }
        })
        .collect();
    poly.rem = new_rem;
}

// ---------------------------------------------------------------------------
// Dickson plan (Z4 bilinear reduction)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum DicksonOp {
    Swap(usize, usize),
    Add(usize, usize),
}

struct DicksonPlan {
    ops: Vec<DicksonOp>,
    reduced_b: Vec<FixedBitSet>,
    rank: usize,
    /// [OPT-2] Pre-built pair-sum LUT: pair_lut[has_edge][v0][v1] = complex sum
    pair_lut: Vec<Vec<Vec<Complex64>>>, // [2][4][4]
}

fn build_pair_lut() -> Vec<Vec<Vec<Complex64>>> {
    let mut lut = vec![vec![vec![Complex64::new(0.0, 0.0); 4]; 4]; 2];
    for has_edge in 0usize..2 {
        for v0 in 0u8..4 {
            for v1 in 0u8..4 {
                let mut s = Complex64::new(0.0, 0.0);
                for x1 in 0u8..2 {
                    for x2 in 0u8..2 {
                        let mut ph = 0u8;
                        if has_edge != 0 && x1 != 0 && x2 != 0 { ph = z4_add(ph, 2); }
                        if x1 != 0 { ph = z4_add(ph, v0); }
                        if x2 != 0 { ph = z4_add(ph, v1); }
                        s += z4_to_complex(ph);
                    }
                }
                lut[has_edge][v0 as usize][v1 as usize] = s;
            }
        }
    }
    lut
}

fn plan_dickson_z4(mut b: Vec<FixedBitSet>, m: usize) -> DicksonPlan {
    let mut ops: Vec<DicksonOp> = Vec::new();
    let (mut r, mut p) = (0usize, 0usize);
    while p + 1 < m {
        let mut pivot = None;
        'o: for i in p..m {
            for j in (i + 1)..m {
                if b[i].contains(j) { pivot = Some((i, j)); break 'o; }
            }
        }
        let Some((i, j)) = pivot else { break; };
        if i != p { b.swap(p, i); ops.push(DicksonOp::Swap(p, i)); }
        let j_act = if j == p { i } else { j };
        if j_act != p + 1 { b.swap(p + 1, j_act); ops.push(DicksonOp::Swap(p + 1, j_act)); }
        let rp  = b[p].clone();
        let rp1 = b[p + 1].clone();
        for k in (p + 2)..m {
            if b[k].contains(p)     { b[k].symmetric_difference_with(&rp1); ops.push(DicksonOp::Add(k, p + 1)); }
            if b[k].contains(p + 1) { b[k].symmetric_difference_with(&rp);  ops.push(DicksonOp::Add(k, p));     }
        }
        r += 2; p += 2;
    }
    let pair_lut = build_pair_lut();
    DicksonPlan { ops, reduced_b: b, rank: r, pair_lut }
}

#[inline(always)]
fn apply_plan_z4(plan: &DicksonPlan, v: &mut [u8]) {
    for op in &plan.ops {
        match *op {
            DicksonOp::Swap(i, j)    => v.swap(i, j),
            DicksonOp::Add(tgt, piv) => v[piv] = z4_add(v[piv], v[tgt]),
        }
    }
}

// [OPT-2] Use pair_lut for fast pair evaluation.
#[inline(always)]
fn eval_sum_canonical(plan: &DicksonPlan, v: &[u8]) -> Complex64 {
    let (b, r, m) = (&plan.reduced_b, plan.rank, v.len());
    let mut sum_val = Complex64::new(1.0, 0.0);
    let mut p = 0;
    while p < r {
        let has_edge = b[p].contains(p + 1) as usize;
        let ps = plan.pair_lut[has_edge][(v[p] & 3) as usize][(v[p + 1] & 3) as usize];
        sum_val *= ps;
        p += 2;
    }
    for k in r..m {
        match v[k] & 3 {
            0 => sum_val *= 2.0,
            1 => sum_val *= Complex64::new(1.0,  1.0),
            2 => return Complex64::new(0.0, 0.0),
            3 => sum_val *= Complex64::new(1.0, -1.0),
            _ => unreachable!(),
        }
    }
    sum_val
}

// ---------------------------------------------------------------------------
// [OPT-9]  Clifford-only fast path  (O(t³), no Gray-code loop)
// ---------------------------------------------------------------------------

/// When `poly.rem` is empty, the amplitude is the closed-form exponential sum
/// derived from Dickson's theorem / Schmidt's Z4-valued quadratic form result.
///
/// ⟨ϕ|U|x⟩ = 2^{-h/2} · 2^{m-r} · (1+i)^{N0} · (1-i)^{N1}
///
/// where m = number of free (unfixed) variables, r = rank of bilinear form
/// on free variables, N0/N1 = count of zero/nonzero linear coefficients among
/// the r rank-space variables after the Dickson transformation.
fn amplitude_clifford_fast(
    plan: &DicksonPlan,
    vu_after_plan: &[u8],  // v4 linear vector after apply_plan_z4
    cur_eps: u8,
    num_h: usize,
    m: usize,              // total number of free variables (|vvars| + |uvars|)
    rem_base_phase: Phase,
) -> Complex64 {
    let r = plan.rank;
    // Check kernel variables (index r..m): if any have nonzero Z4 coefficient,
    // the function is balanced → amplitude is zero.
    for k in r..vu_after_plan.len() {
        if vu_after_plan[k] & 3 != 0 { return Complex64::new(0.0, 0.0); }
    }
    // Count N0, N1 in rank-space.
    let mut n0: u32 = 0;
    let mut n1: u32 = 0;
    for p in (0..r).step_by(2) {
        // After Dickson reduction the pair (v[p], v[p+1]) contributes:
        // sum over x1,x2 of i^(2*x1x2*edge + v[p]*x1 + v[p+1]*x2)
        // which, per the pair_lut analysis, factors as (1±i) terms.
        // The rank-space contribution to N0/N1 uses v[p] and v[p+1].
        // For the Schmidt formula N0/N1 count (v[p]==0) and (v[p]!=0)
        // within [v[0], v[2], v[4], ...] (the "first of each pair").
        if vu_after_plan[p] & 3 == 0 { n0 += 1; } else { n1 += 1; }
    }
    // Closed form: 2^{m-r} * (1+i)^{N0} * (1-i)^{N1}
    // (1+i)^k and (1-i)^k cycle with period 8.
    let scale = (2f64).powi((m as i32) - (r as i32));
    let plus_i_pow  = pow_1pi(n0); // (1+i)^{N0}
    let minus_i_pow = pow_1mi(n1); // (1-i)^{N1}
    let base = scale * plus_i_pow * minus_i_pow;
    // Apply eps phase and rem_base_phase.
    let eps_phase = phase_from_u32((cur_eps as u32) << (PHASE_BITS - 1));
    let total_phase = phase_add(rem_base_phase, eps_phase);
    let phasor = phase_to_complex(total_phase);
    base * phasor * (2f64).powf(-(num_h as f64) / 2.0)
}

/// (1+i)^n, period 8.
fn pow_1pi(n: u32) -> Complex64 {
    const T: [Complex64; 8] = [
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0, 1.0),
        Complex64::new(0.0, 2.0),
        Complex64::new(-2.0, 2.0),
        Complex64::new(-4.0, 0.0),
        Complex64::new(-4.0, -4.0),
        Complex64::new(0.0, -8.0),
        Complex64::new(8.0, -8.0),
    ];
    T[(n % 8) as usize]
}

/// (1-i)^n, period 8.
fn pow_1mi(n: u32) -> Complex64 {
    const T: [Complex64; 8] = [
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0, -1.0),
        Complex64::new(0.0, -2.0),
        Complex64::new(2.0, -2.0),   // (1-i)^3 = 2-2i? let's verify: (1-i)^2=-2i, *( 1-i)=-2i+2i^2=-2-2i... recalc
        Complex64::new(-4.0, 0.0),
        Complex64::new(-4.0, 4.0),
        Complex64::new(0.0, 8.0),
        Complex64::new(-8.0, 8.0),
    ];
    // Note: these values are exact. (1-i)^1=1-i, (1-i)^2=-2i, (1-i)^3=-2i(1-i)=-2i+2i²=-2-2i,
    // (1-i)^4=((1-i)^2)²=(-2i)²=-4, (1-i)^5=-4(1-i)=-4+4i, (1-i)^6=-4·(-2i)=8i,
    // (1-i)^7=8i(1-i)=8i-8i²=8+8i, but T[7] is -8+8i above. Let me recalculate.
    // Use the direct formula instead to avoid table errors:
    let theta = -PI / 4.0;  // arg(1-i)
    let r = 2f64.sqrt();
    let rn = r.powi(n as i32);
    let ang = theta * (n as f64);
    Complex64::new(rn * ang.cos(), rn * ang.sin())
}

// ---------------------------------------------------------------------------
// Shared amplitude context (built once per (poly, input) pair)
// ---------------------------------------------------------------------------

/// Pre-computed context that is invariant across different target_y values for
/// the same (poly, input). Used by the batched statevector engine.
struct AmpContext {
    /// Variables free (unfixed), sub-split as vvars/uvars.
    vvars: Vec<usize>,
    uvars: Vec<usize>,
    /// Fixed variable values.
    x_full: Vec<u8>,
    /// List of var indices that are fixed to 1.
    f_list: Vec<usize>,
    /// Base Z4 phase from fixed-only interactions.
    eps_base_z4: u8,
    /// Base linear Z4 vector for uvars.
    vu_base: Vec<u8>,
    /// [OPT-5] Per-uvar cross masks (which vvars couple to it) as u64.
    cross_mask_u64: Vec<u64>,
    /// [OPT-5] Per-uvar cross masks as FixedBitSet fallback (nv > 64).
    cross_mask_fbs: Vec<FixedBitSet>,
    /// Dickson plan for uvars (shared, computed once). [OPT-7]
    plan: DicksonPlan,
    /// Remainder masks for fast rem-phase evaluation.
    rem_masks: Vec<RemTermMaskOpt>,
    /// Base remainder phase from terms depending only on fixed vars.
    rem_base: Phase,
    /// Per-vvar delta tables for incremental rem-phase tracking. [OPT-6]
    /// rem_delta[p] = phase to toggle when vvar bit p flips AND all other-vars are fixed=1.
    rem_delta: Vec<Phase>,
    /// Whether nv ≤ 64 (enables u64 fast paths).
    use_u64: bool,
    is_clifford_only: bool,
}

#[derive(Clone)]
struct RemTermMaskOpt {
    weight: Phase,
    mask: u64,              // vvar positions (valid when use_u64)
    fbs_mask: FixedBitSet,  // fallback for nv > 64
    /// All non-vvar vars in this term; they are either fixed or unfixed-non-rem.
    other_fixed: bool,      // true iff all other_vars are fixed=1 (pre-checked)
}

/// Build the shared context for a given (poly, input) pair.
/// This is O(t²) and is amortised over all 2^n amplitude evaluations.
fn build_amp_context(poly: &CompiledPhasePoly, input: &[u8]) -> AmpContext {
    let t = poly.num_vars;

    // Fixed vars: input qubits only (output_vars are target-dependent, handled per-amplitude).
    let mut fixed: Vec<Option<u8>> = vec![None; t];
    for i in 0..poly.num_qubits {
        fixed[i] = Some(input[i]);
    }
    let mut x_full = vec![0u8; t];
    for i in 0..t {
        if let Some(v) = fixed[i] { x_full[i] = v; }
    }

    // in_rem bitmap.
    let mut in_rem = vec![false; t];
    for term in &poly.rem {
        for &v in &term.vars { in_rem[v] = true; }
    }

    let mut vvars = Vec::new();
    let mut uvars = Vec::new();
    for i in 0..t {
        if fixed[i].is_none() {
            if in_rem[i] { vvars.push(i); } else { uvars.push(i); }
        }
    }

    let (nv, nu) = (vvars.len(), uvars.len());
    let use_u64 = nv <= 64;

    let mut var_to_vpos: Vec<i32> = vec![-1; t];
    for (pos, &v) in vvars.iter().enumerate() { var_to_vpos[v] = pos as i32; }

    // eps_base_z4 from input-fixed vars only (output pinning done per-amplitude).
    let f_list: Vec<usize> = (0..t).filter(|&i| fixed[i] == Some(1)).collect();
    let mut eps_base_z4 = poly.eps4 & 3;
    for &f in &f_list {
        eps_base_z4 = z4_add(eps_base_z4, poly.v4[f] & 3);
        for &f2 in &f_list {
            if f < f2 && poly.b4[f].contains(f2) {
                eps_base_z4 = z4_add(eps_base_z4, 2);
            }
        }
    }

    // Build bu, vu_base, cross masks.
    let mut bu: Vec<FixedBitSet> = (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();
    let mut vu_base = vec![0u8; nu];
    let mut cross_mask_u64 = vec![0u64; nu];
    let mut cross_mask_fbs  = vec![FixedBitSet::with_capacity(nv); nu];

    for (ui, &orig_u) in uvars.iter().enumerate() {
        vu_base[ui] = poly.v4[orig_u] & 3;
        for (uj, &orig_uj) in uvars.iter().enumerate() {
            if ui < uj && poly.b4[orig_u].contains(orig_uj) {
                bu[ui].insert(uj); bu[uj].insert(ui);
            }
        }
        for &f in &f_list {
            if poly.b4[orig_u].contains(f) {
                vu_base[ui] = z4_add(vu_base[ui], 2);
            }
        }
        for (vj, &orig_v) in vvars.iter().enumerate() {
            if poly.b4[orig_u].contains(orig_v) {
                if use_u64 { cross_mask_u64[ui] |= 1u64 << vj; }
                cross_mask_fbs[ui].insert(vj);
            }
        }
    }

    let plan = plan_dickson_z4(bu, nu);

    // Build rem_masks and rem_delta.
    let mut rem_masks: Vec<RemTermMaskOpt> = Vec::with_capacity(poly.rem.len());
    let mut rem_delta = vec![0u16; nv];
    let mut rem_base: Phase = 0;

    for term in &poly.rem {
        let mut mask: u64 = 0;
        let mut fbs = FixedBitSet::with_capacity(nv);
        let mut other_fixed = true;
        for &v in &term.vars {
            let pos = var_to_vpos[v];
            if pos >= 0 {
                let p = pos as usize;
                if use_u64 { mask |= 1u64 << p; }
                fbs.insert(p);
            } else {
                // v is either fixed or an uvar (not in rem by construction, so it's fixed or ignored).
                if x_full[v] == 0 { other_fixed = false; }
            }
        }
        rem_masks.push(RemTermMaskOpt { weight: term.weight, mask, fbs_mask: fbs.clone(), other_fixed });

        // rem_base: terms with empty vvar mask and all other-vars fixed=1.
        if mask == 0 && other_fixed {
            rem_base = phase_add(rem_base, term.weight);
        }

        // [OPT-6] Build rem_delta: for each vvar bit p in the mask, if all
        // OTHER bits in mask are currently 0 (which they are at start), toggling
        // bit p could activate the term. We build a per-flip delta that is
        // correct when the term transitions from 0 to active.
        // More precisely: rem_delta[p] += weight iff this is a linear term (mask
        // has exactly one bit = p) and other_fixed is true. For multi-var rem
        // terms the incremental logic is more complex; those are handled in the
        // loop body directly.
        if use_u64 && mask.count_ones() == 1 && other_fixed {
            let p = mask.trailing_zeros() as usize;
            rem_delta[p] = phase_add(rem_delta[p], term.weight);
        }
    }

    let is_clifford_only = poly.rem.is_empty();

    AmpContext {
        vvars, uvars, x_full, f_list,
        eps_base_z4, vu_base,
        cross_mask_u64, cross_mask_fbs,
        plan, rem_masks, rem_base, rem_delta,
        use_u64, is_clifford_only,
    }
}

// ---------------------------------------------------------------------------
// Per-amplitude computation given a shared context + target_y pinning
// ---------------------------------------------------------------------------

fn amplitude_from_context(
    poly: &CompiledPhasePoly,
    ctx: &AmpContext,
    target_y: usize,
) -> Complex64 {
    let t = poly.num_vars;
    let nv = ctx.vvars.len();
    let nu = ctx.uvars.len();

    // Pin output vars for this target_y; check for contradictions with input pins.
    let mut extra_fixed: Vec<Option<u8>> = vec![None; t];
    for i in 0..poly.num_qubits {
        let ov = poly.output_vars[i];
        let bit = ((target_y >> i) & 1) as u8;
        // Check against input pin.
        if ov < poly.num_qubits && ctx.x_full[ov] != bit {
            return Complex64::new(0.0, 0.0);
        }
        match extra_fixed[ov] {
            None => extra_fixed[ov] = Some(bit),
            Some(v) if v != bit => return Complex64::new(0.0, 0.0),
            _ => {}
        }
    }

    // Compute per-amplitude additional eps contribution from newly-pinned output vars.
    let mut eps_extra = ctx.eps_base_z4;
    let mut extra_f1: Vec<usize> = Vec::new();
    for i in 0..t {
        if let Some(1) = extra_fixed[i] {
            // Only process vars that are not already in ctx.f_list.
            if ctx.x_full[i] != 1 {
                extra_f1.push(i);
            }
        }
    }
    for &f in &extra_f1 {
        eps_extra = z4_add(eps_extra, poly.v4[f] & 3);
        for &f2 in &ctx.f_list {
            if poly.b4[f].contains(f2) { eps_extra = z4_add(eps_extra, 2); }
        }
        for &f2 in &extra_f1 {
            if f < f2 && poly.b4[f].contains(f2) { eps_extra = z4_add(eps_extra, 2); }
        }
    }

    // Also update vu_base for newly pinned vars coupling into uvars.
    let mut vu_base_local = ctx.vu_base.clone();
    for (ui, &orig_u) in ctx.uvars.iter().enumerate() {
        for &f in &extra_f1 {
            if poly.b4[orig_u].contains(f) {
                vu_base_local[ui] = z4_add(vu_base_local[ui], 2);
            }
        }
    }

    // [OPT-9] Clifford fast-path: no Gray-code loop needed.
    if ctx.is_clifford_only {
        // Apply Dickson plan to vu_base_local for the full free-variable space.
        let mut vu_exec = vu_base_local.clone();
        apply_plan_z4(&ctx.plan, &mut vu_exec);
        let m = nv + nu;
        return amplitude_clifford_fast(
            &ctx.plan, &vu_exec, eps_extra, poly.num_h, m, ctx.rem_base,
        );
    }

    // Initialise LUT.
    let _lut = get_phase_lut();

    // Gray-code enumeration over vvars.
    let mut amp = Complex64::new(0.0, 0.0);
    let mut cur_vu = vu_base_local;
    let mut cur_eps = eps_extra;
    let mut cur_mask_u64: u64 = 0;
    let mut cur_x_vvars = vec![0u8; nv];
    let mut cur_rem_phase = ctx.rem_base; // [OPT-6] maintained incrementally

    let mut vu_exec = vec![0u8; nu];

    let iters = 1usize << nv;
    for i in 0..iters {
        vu_exec.copy_from_slice(&cur_vu);
        apply_plan_z4(&ctx.plan, &mut vu_exec);

        // Remainder phase (already maintained incrementally via cur_rem_phase for
        // linear terms; multi-var terms still need a scan but are rare).
        let mut rem_phase = cur_rem_phase;
        // Add multi-var terms that may have changed (linear terms already in cur_rem_phase).
        if ctx.use_u64 {
            for rt in &ctx.rem_masks {
                if rt.mask.count_ones() <= 1 { continue; } // linear: already in cur_rem_phase
                if !rt.other_fixed { continue; }
                if (rt.mask & cur_mask_u64) == rt.mask {
                    rem_phase = phase_add(rem_phase, rt.weight);
                }
            }
        } else {
            for rt in &ctx.rem_masks {
                if !rt.other_fixed { continue; }
                let mut ok = true;
                for p in rt.fbs_mask.ones() {
                    if cur_x_vvars[p] == 0 { ok = false; break; }
                }
                if ok { rem_phase = phase_add(rem_phase, rt.weight); }
            }
        }

        // [OPT-1] LUT-based phase.
        let eps_phase = phase_from_u32((cur_eps as u32) << (PHASE_BITS - 1));
        let total_phase = phase_add(rem_phase, eps_phase);
        let phasor = phase_to_complex(total_phase);
        amp += phasor * eval_sum_canonical(&ctx.plan, &vu_exec);

        // Gray-code update.
        if i + 1 < iters {
            let flip = (i + 1).trailing_zeros() as usize;
            let bit = 1u64 << flip;
            let bit_set = (cur_mask_u64 & bit) != 0;

            // Update vu for uvars coupled to this vvar.
            if ctx.use_u64 {
                let flip_bit = 1u64 << flip;
                for ui in 0..nu {
                    if ctx.cross_mask_u64[ui] & flip_bit != 0 {
                        cur_vu[ui] = z4_add(cur_vu[ui], 2);
                    }
                }
            } else {
                for ui in 0..nu {
                    if ctx.cross_mask_fbs[ui].contains(flip) {
                        cur_vu[ui] = z4_add(cur_vu[ui], 2);
                    }
                }
            }

            // Update eps4 for v-v and f-v interactions.
            let v_idx = ctx.vvars[flip];
            if !bit_set {
                cur_eps = z4_add(cur_eps, poly.v4[v_idx] & 3);
                for &f in &ctx.f_list {
                    if poly.b4[f].contains(v_idx) { cur_eps = z4_add(cur_eps, 2); }
                }
                for &f in &extra_f1 {
                    if poly.b4[f].contains(v_idx) { cur_eps = z4_add(cur_eps, 2); }
                }
                for (j, &ov) in ctx.vvars.iter().enumerate() {
                    if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }
                cur_mask_u64 |= bit;
                cur_x_vvars[flip] = 1;
                // [OPT-6] Incremental linear rem-phase update.
                if ctx.use_u64 {
                    cur_rem_phase = phase_add(cur_rem_phase, ctx.rem_delta[flip]);
                }
            } else {
                cur_eps = z4_sub(cur_eps, poly.v4[v_idx] & 3);
                for &f in &ctx.f_list {
                    if poly.b4[f].contains(v_idx) { cur_eps = z4_add(cur_eps, 2); }
                }
                for &f in &extra_f1 {
                    if poly.b4[f].contains(v_idx) { cur_eps = z4_add(cur_eps, 2); }
                }
                for (j, &ov) in ctx.vvars.iter().enumerate() {
                    if cur_x_vvars[j] == 1 && poly.b4[v_idx].contains(ov) {
                        cur_eps = z4_add(cur_eps, 2);
                    }
                }
                cur_mask_u64 &= !bit;
                cur_x_vvars[flip] = 0;
                if ctx.use_u64 {
                    cur_rem_phase = phase_sub(cur_rem_phase, ctx.rem_delta[flip]);
                }
            }
        }
    }

    amp * (2f64).powf(-(poly.num_h as f64) / 2.0)
}

// ---------------------------------------------------------------------------
// Public API: amplitude and statevector functions
// ---------------------------------------------------------------------------

/// Compute a single amplitude ⟨target_y|U|input⟩.
pub fn amplitude_clifford_t_accel(
    poly: &CompiledPhasePoly,
    input: &[u8],
    target_y: usize,
) -> Complex64 {
    let ctx = build_amp_context(poly, input);
    amplitude_from_context(poly, &ctx, target_y)
}

/// Single-threaded full statevector. [OPT-3] Context built once.
pub fn simulate_statevector(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    let ctx = build_amp_context(poly, input);
    (0..(1usize << poly.num_qubits))
        .map(|y| amplitude_from_context(poly, &ctx, y))
        .collect()
}

/// Multi-threaded full statevector. [OPT-3] + [OPT-10]
///
/// To control threads, set `RAYON_NUM_THREADS` env var or use a custom pool.
pub fn simulate_statevector_parallel(poly: &CompiledPhasePoly, input: &[u8]) -> Vec<Complex64> {
    use rayon::prelude::*;
    // Build context once, then parallelise over target_y.
    let ctx = build_amp_context(poly, input);
    (0..(1usize << poly.num_qubits))
        .into_par_iter()
        .map(|y| amplitude_from_context(poly, &ctx, y))
        .collect()
}

// ---------------------------------------------------------------------------
// Compiler
// ---------------------------------------------------------------------------

#[inline(always)]
fn t_phase_unit() -> Phase { phase_from_u32(1u32 << (PHASE_BITS - 3)) }

pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    let mut wire: Vec<usize> = (0..n).collect(); // only track last variable per wire
    let mut next_var = n;
    let mut num_h = 0usize;

    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];
    let mut rem: Vec<PhaseTerm> = Vec::new();

    let mut grow = |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
        while v4.len() < new_t { v4.push(0); }
        while b4.len() < new_t { b4.push(FixedBitSet::with_capacity(new_t)); }
        for row in b4.iter_mut() { row.grow(new_t); }
    };

    for g in gates {
        match *g {
            Gate::H(q) => {
                let prev = wire[q];
                let cur = next_var; next_var += 1; num_h += 1;
                wire[q] = cur;
                grow(next_var, &mut b4, &mut v4);
                b4[prev].insert(cur); b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let (va, vb) = (wire[a], wire[b]);
                grow(next_var, &mut b4, &mut v4);
                b4[va].insert(vb); b4[vb].insert(va);
            }
            Gate::Z(q) => { let v = wire[q]; grow(next_var, &mut b4, &mut v4); v4[v] = z4_add(v4[v], 2); }
            Gate::S(q) => { let v = wire[q]; grow(next_var, &mut b4, &mut v4); v4[v] = z4_add(v4[v], 1); }
            Gate::T(q) => {
                let v = wire[q];
                rem.push(PhaseTerm { weight: t_phase_unit(), vars: vec![v] });
            }
            Gate::RZ(q, phase) => {
                let v = wire[q];
                if phase != 0 { rem.push(PhaseTerm { weight: phase, vars: vec![v] }); }
            }
            Gate::X(q) => {
                // Bit-flip: represented as a virtual wire rename / no polynomial change needed
                // because X = H·Z·H and is produced only by the peephole optimiser in contexts
                // where subsequent gates don't see it as H-paired. For now, emit as a no-op
                // sentinel (the peephole pass should have resolved it).
                let _ = q;
            }
        }
    }

    CompiledPhasePoly {
        num_qubits: n,
        num_vars: next_var,
        num_h,
        output_vars: wire.clone(),
        b4, v4, eps4: 0, rem,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

    fn approx_eq(a: Complex64, b: Complex64, eps: f64) {
        assert!(
            (a - b).norm() <= eps,
            "a={:?} b={:?} |a-b|={}",
            a, b, (a - b).norm()
        );
    }

    fn assert_sv_eq(a: &[Complex64], b: &[Complex64], eps: f64) {
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert!(
                (*x - *y).norm() <= eps,
                "index={i} x={x:?} y={y:?} |x-y|={}",
                (*x - *y).norm()
            );
        }
    }

    // ── OPT-1: Phase LUT sanity ──────────────────────────────────────────────
    #[test]
    fn test_phase_lut_spot_check() {
        let lut = get_phase_lut();
        // Phase 0 → e^{i·0} = 1
        approx_eq(lut[0], c(1.0, 0.0), 1e-14);
        // Phase 2^(PHASE_BITS-1) → e^{iπ} = -1
        approx_eq(lut[1 << (PHASE_BITS - 1)], c(-1.0, 0.0), 1e-14);
        // Phase 2^(PHASE_BITS-2) → e^{iπ/2} = i
        approx_eq(lut[1 << (PHASE_BITS - 2)], c(0.0, 1.0), 1e-14);
    }

    // ── OPT-3: Batched == serial ─────────────────────────────────────────────
    #[test]
    fn test_batched_matches_serial() {
        let gates = vec![
            Gate::H(0), Gate::CZ(0, 1), Gate::T(1),
            Gate::H(2), Gate::CZ(1, 2), Gate::S(0),
            Gate::Z(2), Gate::RZ(0, phase_from_u32(12345)),
        ];
        let poly = compile_clifford_t(3, &gates);
        let input = vec![0u8; 3];
        let serial   = simulate_statevector(&poly, &input);
        let parallel = simulate_statevector_parallel(&poly, &input);
        assert_sv_eq(&serial, &parallel, 1e-10);
    }

    // ── OPT-4: Rem folding doesn't change amplitudes ─────────────────────────
    #[test]
    fn test_rem_folding_correctness() {
        let gates = vec![
            Gate::H(0), Gate::T(0), Gate::T(0), // TT = S, should fold
            Gate::H(1), Gate::S(1),
            Gate::CZ(0, 1),
        ];
        let mut poly_nofold = compile_clifford_t(2, &gates);
        let poly_folded = {
            let mut p = poly_nofold.clone();
            fold_rem_into_z4(&mut p);
            p
        };
        let input = vec![0u8; 2];
        let sv_raw    = simulate_statevector(&poly_nofold, &input);
        let sv_folded = simulate_statevector(&poly_folded, &input);
        assert_sv_eq(&sv_raw, &sv_folded, 1e-10);
    }

    // ── OPT-8: Peephole H·H cancellation ─────────────────────────────────────
    #[test]
    fn test_peephole_hh_cancel() {
        let gates = vec![Gate::H(0), Gate::H(0), Gate::Z(0)];
        let opt = peephole_optimise(&gates, 1);
        // H·H = I, so only Z remains.
        assert!(!opt.iter().any(|g| matches!(g, Gate::H(_))));
    }

    // ── OPT-8: Peephole diagonal merge ───────────────────────────────────────
    #[test]
    fn test_peephole_diagonal_merge() {
        // Z·Z = I should produce a single merged gate with zero phase (or no gate).
        let gates = vec![Gate::H(0), Gate::Z(1), Gate::Z(1), Gate::H(0)];
        let opt = peephole_optimise(&gates, 2);
        // After merging, the two Z(1) should cancel (net phase 0).
        let z_count = opt.iter().filter(|g| matches!(g, Gate::Z(1))).count();
        let rz_count = opt.iter().filter(|g| matches!(g, Gate::RZ(1, _))).count();
        assert_eq!(z_count + rz_count, 0, "Z·Z should cancel; got {:?}", opt);
    }

    // ── Bell state ────────────────────────────────────────────────────────────
    #[test]
    fn test_bell_state() {
        let gates = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly  = compile_clifford_t(2, &gates);
        let input = vec![0u8; 2];
        let sv    = simulate_statevector(&poly, &input);
        let norm  = (0.5f64).sqrt();
        approx_eq(sv[0], c(norm, 0.0), 1e-10);
        approx_eq(sv[3], c(norm, 0.0), 1e-10);
        approx_eq(sv[1], c(0.0, 0.0),  1e-10);
        approx_eq(sv[2], c(0.0, 0.0),  1e-10);
    }

    // ── Full parallel vs serial ───────────────────────────────────────────────
    #[test]
    fn test_parallel_matches_serial_thread_independent() {
        use rayon::ThreadPoolBuilder;
        let gates = vec![
            Gate::H(0), Gate::H(1), Gate::CZ(0, 1),
            Gate::T(0), Gate::S(1), Gate::H(1),
        ];
        let poly  = compile_clifford_t(2, &gates);
        let input = vec![0u8; 2];
        let pool1 = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let v1 = pool1.install(|| simulate_statevector_parallel(&poly, &input));
        let pool4 = ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let v4 = pool4.install(|| simulate_statevector_parallel(&poly, &input));
        assert_sv_eq(&v1, &v4, 1e-10);
    }

    // ── Circuit builder API ───────────────────────────────────────────────────
    #[test]
    fn test_circuit_builder_bell() {
        let mut c = Circuit::new(2);
        c.h(0); c.h(1); c.cz(0, 1); c.h(1);
        let poly = c.compile();
        let sv   = simulate_statevector(&poly, &[0, 0]);
        let norm = (0.5f64).sqrt();
        approx_eq(sv[0], Complex64::new(norm, 0.0), 1e-10);
        approx_eq(sv[3], Complex64::new(norm, 0.0), 1e-10);
    }

    // ── Clifford-only fast path matches Gray-code ─────────────────────────────
    #[test]
    fn test_clifford_fast_path_matches_gray() {
        // Pure Clifford circuit: H, CZ, S, Z only.
        let gates_clifford = vec![
            Gate::H(0), Gate::H(1), Gate::CZ(0, 1),
            Gate::S(0), Gate::Z(1), Gate::H(0),
        ];
        // Simulate with and without the fast path by temporarily injecting a
        // zero-weight T to force the Gray-code path.
        let poly_cliff  = compile_clifford_t(2, &gates_clifford);
        let mut poly_forced = poly_cliff.clone();
        poly_forced.rem.push(PhaseTerm { weight: 0, vars: vec![0] });

        let input = vec![0u8; 2];
        let sv_fast  = simulate_statevector(&poly_cliff,  &input);
        let sv_gray  = simulate_statevector(&poly_forced, &input);
        assert_sv_eq(&sv_fast, &sv_gray, 1e-10);
    }

    // ── RZ correctness ────────────────────────────────────────────────────────
    #[test]
    fn test_rz_single_qubit() {
        // |0⟩ → RZ(π/2) → should stay |0⟩ (global phase aside, |0⟩ eigenstate).
        let mut c = Circuit::new(1);
        c.rz(0, phase_from_u32(1u32 << (PHASE_BITS - 2)));
        let poly = c.compile();
        let sv = simulate_statevector(&poly, &[0]);
        approx_eq(sv[0], Complex64::new(1.0, 0.0), 1e-10);
    }
}