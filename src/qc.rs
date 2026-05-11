use crate::{compile_clifford_t, Gate, Phase, PHASE_BITS, CompiledPhasePoly};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// A very small error type for the trivial QASM parser we ship with the
/// library.  Only the `h`, `z`, `s`, `t`, `rz` and `cz` instructions are
/// recognised; everything else (including comments) is silently ignored.
#[derive(Debug)]
pub enum QasmError {
    Io(io::Error),
    Parse(String),
}

impl From<io::Error> for QasmError {
    fn from(err: io::Error) -> QasmError {
        QasmError::Io(err)
    }
}

/// Convenient in‑memory representation of a circuit.
#[derive(Clone, Debug, Default)]
pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    /// Start with a clean circuit acting on `num_qubits` wires.
    pub fn new(num_qubits: usize) -> Self {
        Self { num_qubits, gates: Vec::new() }
    }

    pub fn from_gates(num_qubits: usize, gates: Vec<Gate>) -> Self {
        Self { num_qubits, gates }
    }

    // builder helpers --------------------------------------------------
    pub fn h(&mut self, q: usize) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::H(q));
    }
    pub fn z(&mut self, q: usize) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::Z(q));
    }
    pub fn s(&mut self, q: usize) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::S(q));
    }
    pub fn t(&mut self, q: usize) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::T(q));
    }
    /// Fixed-precision dyadic RZ: `phase` is in units of π/2^PHASE_BITS.
    pub fn rz(&mut self, q: usize, phase: Phase) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::RZ(q, phase));
    }
    pub fn cz(&mut self, a: usize, b: usize) {
        self.num_qubits = self.num_qubits.max(a + 1).max(b + 1);
        self.gates.push(Gate::CZ(a, b));
    }

    /// Compile this circuit into the usual phase‑polynomial form used by the simulator.
    pub fn compile(&self) -> CompiledPhasePoly {
        compile_clifford_t(self.num_qubits, &self.gates)
    }
}

// QASM <-> Circuit conversions ----------------------------------------

fn parse_qubit(token: &str) -> Result<usize, QasmError> {
    if let Some(start) = token.find('[') {
        if let Some(end) = token.find(']') {
            let num = &token[start + 1..end];
            return num
                .parse()
                .map_err(|e| QasmError::Parse(format!("invalid qubit index: {}", e)));
        }
    }
    Err(QasmError::Parse(format!("bad qubit token `{}`", token)))
}

/// Parse an RZ argument in a speed-focused, dyadic-friendly way.
///
/// Supported syntaxes:
/// - `rz k q[0];`            where k is an integer in [0,2^PHASE_BITS)
/// - `rz k/2^PHASE_BITS q[0];`  (we accept `2^16` literally)
///
/// We interpret the argument as `phase` in units of π/2^PHASE_BITS.
/// So `rz 16384 q[0];` is π/4 when PHASE_BITS=16.
fn parse_rz_phase(token: &str) -> Result<Phase, QasmError> {
    let tok = token.trim();

    // form: "k"
    if let Ok(v) = tok.parse::<u32>() {
        return Ok((v & ((1u32 << PHASE_BITS) - 1)) as Phase);
    }

    // form: "k/2^16"
    if let Some((num, den)) = tok.split_once('/') {
        let num_v: u32 = num.trim().parse().map_err(|e| QasmError::Parse(format!("bad rz numerator: {}", e)))?;
        let den = den.trim();
        if den.starts_with("2^") {
            let pow_str = &den[2..];
            let pow: u32 = pow_str.parse().map_err(|e| QasmError::Parse(format!("bad rz denominator power: {}", e)))?;
            if pow != PHASE_BITS {
                return Err(QasmError::Parse(format!(
                    "rz denominator must be 2^{} for this build (got 2^{})",
                    PHASE_BITS, pow
                )));
            }
            return Ok((num_v & ((1u32 << PHASE_BITS) - 1)) as Phase);
        }
    }

    Err(QasmError::Parse(format!("bad rz argument `{}`", token)))
}

/// Read a text file containing a (restricted) OpenQASM-like circuit and return a [`Circuit`].
pub fn read_qasm_file(path: &str) -> Result<Circuit, QasmError> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let mut circ = Circuit::default();

    for line in reader.lines() {
        let line = line?;
        for stmt in line.split(';') {
            let stmt = stmt.trim();
            if stmt.is_empty() {
                continue;
            }
            let stmt = stmt.replace(',', " ");
            let parts: Vec<&str> = stmt.split_whitespace().filter(|s| !s.is_empty()).collect();
            if parts.is_empty() || parts[0].starts_with('%') || parts[0].starts_with("//") {
                continue;
            }

            match parts.as_slice() {
                ["qubit", q] => {
                    if let Ok(n) = parse_qubit(q) {
                        circ.num_qubits = circ.num_qubits.max(n + 1);
                    }
                }
                ["h", q] => circ.h(parse_qubit(q)?),
                ["z", q] => circ.z(parse_qubit(q)?),
                ["s", q] => circ.s(parse_qubit(q)?),
                ["t", q] => circ.t(parse_qubit(q)?),
                ["rz", arg, q] => {
                    let qq = parse_qubit(q)?;
                    let ph = parse_rz_phase(arg)?;
                    circ.rz(qq, ph);
                }
                ["cz", a, b] => {
                    let aa = parse_qubit(a)?;
                    let bb = parse_qubit(b)?;
                    circ.cz(aa, bb);
                }
                _ => {
                    // ignore unknown statements
                }
            }
        }
    }

    Ok(circ)
}

/// Return a QASM string describing the circuit.
///
/// For speed and exactness, we emit `rz <phase>` as an integer phase in units π/2^PHASE_BITS.
/// Example (PHASE_BITS=16):
///   rz 16384 q[0];
pub fn write_qasm_string(circ: &Circuit) -> String {
    let mut out = String::new();
    for g in &circ.gates {
        match *g {
            Gate::H(q) => out.push_str(&format!("h q[{}];\n", q)),
            Gate::Z(q) => out.push_str(&format!("z q[{}];\n", q)),
            Gate::S(q) => out.push_str(&format!("s q[{}];\n", q)),
            Gate::T(q) => out.push_str(&format!("t q[{}];\n", q)),
            Gate::RZ(q, ph) => out.push_str(&format!("rz {} q[{}];\n", ph, q)),
            Gate::CZ(a, b) => out.push_str(&format!("cz q[{}],q[{}];\n", a, b)),
        }
    }
    out
}

pub fn write_qasm_file(path: &str, circ: &Circuit) -> io::Result<()> {
    let mut file = File::create(path)?;
    file.write_all(write_qasm_string(circ).as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    // helper to parse directly from string (not exported)
    fn read_qasm_string(src: &str) -> Result<Circuit, QasmError> {
        let mut circ = Circuit::default();
        for line in src.lines() {
            for stmt in line.split(';') {
                let stmt = stmt.trim();
                if stmt.is_empty() {
                    continue;
                }
                let stmt = stmt.replace(',', " ");
                let parts: Vec<&str> = stmt.split_whitespace().filter(|s| !s.is_empty()).collect();
                if parts.is_empty() {
                    continue;
                }
                match parts.as_slice() {
                    ["h", q] => circ.h(parse_qubit(q)?),
                    ["z", q] => circ.z(parse_qubit(q)?),
                    ["s", q] => circ.s(parse_qubit(q)?),
                    ["t", q] => circ.t(parse_qubit(q)?),
                    ["rz", arg, q] => circ.rz(parse_qubit(q)?, parse_rz_phase(arg)?),
                    ["cz", a, b] => circ.cz(parse_qubit(a)?, parse_qubit(b)?),
                    _ => (),
                }
            }
        }
        Ok(circ)
    }

    #[test]
    fn round_trip_qasm() {
        let src = "h q[0]; cz q[0],q[1]; rz 16384 q[1]; h q[1];";
        let circ = read_qasm_string(src).unwrap();
        assert_eq!(circ.num_qubits, 2);
        assert_eq!(
            circ.gates,
            vec![Gate::H(0), Gate::CZ(0, 1), Gate::RZ(1, 16384u16), Gate::H(1)]
        );

        let back = write_qasm_string(&circ);
        assert!(back.contains("h q[0];"));
        assert!(back.contains("cz q[0],q[1];"));
        assert!(back.contains("rz 16384 q[1];"));

        let tmp = std::env::temp_dir().join("polyq_roundtrip.qasm");
        write_qasm_file(tmp.to_str().unwrap(), &circ).unwrap();
        let circ2 = read_qasm_file(tmp.to_str().unwrap()).unwrap();
        assert_eq!(circ2.num_qubits, circ.num_qubits);
        assert_eq!(circ2.gates, circ.gates);
    }
}