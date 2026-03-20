use crate::{compile_clifford_t, CompiledPhasePoly, Gate, Phase, PHASE_BITS};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// A very small error type for the trivial QASM parser we ship with the
/// library.
///
/// Supported instructions:
/// - `h q[i];`
/// - `z q[i];`
/// - `s q[i];`
/// - `t q[i];`
/// - `rz <k> q[i];` where k is an integer phase in units of π/2^PHASE_BITS (now modulo 2π)
/// - `cz q[i],q[j];`
/// - `mcz q[a],q[b],q[c],...;`
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
    pub fn new(num_qubits: usize) -> Self {
        Self { num_qubits, gates: Vec::new() }
    }

    pub fn from_gates(num_qubits: usize, gates: Vec<Gate>) -> Self {
        Self { num_qubits, gates }
    }

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
    pub fn rz(&mut self, q: usize, phase: Phase) {
        self.num_qubits = self.num_qubits.max(q + 1);
        self.gates.push(Gate::RZ(q, phase));
    }
    pub fn cz(&mut self, a: usize, b: usize) {
        self.num_qubits = self.num_qubits.max(a + 1).max(b + 1);
        self.gates.push(Gate::CZ(a, b));
    }
    pub fn mcz(&mut self, ctrls: &[usize]) {
        for &q in ctrls {
            self.num_qubits = self.num_qubits.max(q + 1);
        }
        self.gates.push(Gate::MCZ(ctrls.to_vec()));
    }

    pub fn compile(&self) -> CompiledPhasePoly {
        compile_clifford_t(self.num_qubits, &self.gates)
    }
}

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

fn parse_qubits(tokens: &[&str]) -> Result<Vec<usize>, QasmError> {
    let mut out = Vec::with_capacity(tokens.len());
    for &t in tokens {
        out.push(parse_qubit(t)?);
    }
    Ok(out)
}

/// Parse an RZ argument in a dyadic-friendly way.
///
/// Supported syntaxes:
/// - `rz k q[0];`
/// - `rz k/2^PHASE_BITS q[0];`
///
/// With 2π-periodic phases, the valid integer range is [0, 2^(PHASE_BITS+1)).
fn parse_rz_phase(token: &str) -> Result<Phase, QasmError> {
    let tok = token.trim();

    // form: "k"
    if let Ok(v) = tok.parse::<u32>() {
        // mask to 2^(PHASE_BITS+1)
        let mask = (1u32 << (PHASE_BITS + 1)) - 1;
        return Ok((v & mask) as Phase);
    }

    // form: "k/2^PHASE_BITS"
    if let Some((num, den)) = tok.split_once('/') {
        let num_v: u32 = num
            .trim()
            .parse()
            .map_err(|e| QasmError::Parse(format!("bad rz numerator: {}", e)))?;
        let den = den.trim();
        if den.starts_with("2^") {
            let pow_str = &den[2..];
            let pow: u32 = pow_str
                .parse()
                .map_err(|e| QasmError::Parse(format!("bad rz denominator power: {}", e)))?;

            if pow != PHASE_BITS {
                return Err(QasmError::Parse(format!(
                    "rz denominator must be 2^{} for this build (got 2^{})",
                    PHASE_BITS, pow
                )));
            }

            // Here, "k/2^PHASE_BITS" means k * (π/2^PHASE_BITS), so integer is k.
            let mask = (1u32 << (PHASE_BITS + 1)) - 1;
            return Ok((num_v & mask) as Phase);
        }
    }

    Err(QasmError::Parse(format!("bad rz argument `{}`", token)))
}

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
                ["mcz", rest @ ..] => {
                    let qs = parse_qubits(rest)?;
                    circ.mcz(&qs);
                }
                _ => {}
            }
        }
    }

    Ok(circ)
}

pub fn write_qasm_string(circ: &Circuit) -> String {
    let mut out = String::new();
    for g in &circ.gates {
        match g {
            Gate::H(q) => out.push_str(&format!("h q[{}];\n", q)),
            Gate::Z(q) => out.push_str(&format!("z q[{}];\n", q)),
            Gate::S(q) => out.push_str(&format!("s q[{}];\n", q)),
            Gate::T(q) => out.push_str(&format!("t q[{}];\n", q)),
            Gate::RZ(q, ph) => out.push_str(&format!("rz {} q[{}];\n", ph, q)),
            Gate::CZ(a, b) => out.push_str(&format!("cz q[{}],q[{}];\n", a, b)),
            Gate::MCZ(ctrls) => {
                out.push_str("mcz ");
                for (i, q) in ctrls.iter().enumerate() {
                    if i != 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("q[{}]", q));
                }
                out.push_str(";\n");
            }
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
                    ["mcz", rest @ ..] => {
                        let qs = parse_qubits(rest)?;
                        circ.mcz(&qs);
                    }
                    _ => (),
                }
            }
        }
        Ok(circ)
    }

    #[test]
    fn round_trip_qasm_including_mcz() {
        let src = "h q[0]; mcz q[0],q[1],q[2]; rz 65536 q[1];";
        let circ = read_qasm_string(src).unwrap();
        assert_eq!(circ.num_qubits, 3);
        assert_eq!(circ.gates, vec![Gate::H(0), Gate::MCZ(vec![0, 1, 2]), Gate::RZ(1, 65536u32)]);

        let back = write_qasm_string(&circ);
        assert!(back.contains("mcz q[0],q[1],q[2];"));
        assert!(back.contains("rz 65536 q[1];"));
    }
}