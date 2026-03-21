//! Crate root for PolyQ.  The bulk of the implementation lives in
//! [`sim`] and the small QASM helper module [`qc`].

pub mod qc;
pub mod sim;

// re-export a convenient public API at the crate root so callers can write
// `PolyQ::Gate` rather than `PolyQ::sim::Gate`.
pub use sim::{
    Circuit, CompiledPhasePoly, Gate, PHASE_BITS, Phase, PhaseTerm, QasmError,
    amplitude_clifford_t_accel, compile_clifford_t, read_qasm_file, simulate_statevector,
    write_qasm_file, write_qasm_string,
};
