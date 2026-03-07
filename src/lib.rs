//! Crate root for PolyQ.  The bulk of the implementation lives in
//! [`sim`] and the small QASM helper module [`qc`].

pub mod qc;
pub mod sim;

// re-export a convenient public API at the crate root so callers can write
// `PolyQ::Gate` rather than `PolyQ::sim::Gate`.
pub use sim::{
    Gate,
    Z8Term,
    CompiledPhasePoly,
    compile_clifford_t,
    amplitude_clifford_t_accel,
    simulate_statevector,
    simulate_statevector_clifford_t_multiprocess,
    amplitude_clifford_t_multithreaded,
    Circuit,
    read_qasm_file,
    write_qasm_file,
    write_qasm_string,
    QasmError,
};

