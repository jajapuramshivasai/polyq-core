//! Crate root for PolyQ.  The bulk of the implementation lives in
//! [`sim`] and the small QASM helper module [`qc`].

pub mod qc;
pub mod sim;

// re-export a convenient public API at the crate root so callers can write
// `PolyQ::Gate` rather than `PolyQ::sim::Gate`.
pub use sim::{
    amplitude_clifford_t_accel, compile_clifford_t, simulate_statevector, Circuit, CompiledPhasePoly,
    Gate, Phase, PhaseTerm, PHASE_BITS, read_qasm_file, write_qasm_file, write_qasm_string, QasmError,
};