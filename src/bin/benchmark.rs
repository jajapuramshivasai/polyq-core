use PolyQ::sim::{Circuit, simulate_statevector};
use PolyQ::qc::read_qasm_file;
use std::time::Instant;
use std::fs::File;
use std::io::Read;
use std::path::Path;


//clifford circuits
fn build_bv(n: usize, secret: &[u8]) -> Circuit {

// Build a Bernstein–Vazirani circuit for a secret bitstring.
    let mut circ = Circuit::new(n);

    // initial Hadamards
    for q in 0..n {
        circ.h(q);
    }

    // oracle: apply Z on wires where secret bit = 1
    for (q, &b) in secret.iter().enumerate() {
        if b == 1 {
            circ.z(q);
        }
    }

    // final Hadamards
    for q in 0..n {
        circ.h(q);
    }

    circ
}


//non‑clifford circuits
fn build_qft(n: usize) -> Circuit {
    let mut circ = Circuit::new(n);
    // Apply Hadamard to all qubits
    for i in 0..n {
        circ.h(i);
        for j in (i + 1)..n {
            let k = j - i + 1;
            // Approximate controlled phase exp(i*pi/2^k) using Clifford+T gates
            if k == 1 {
                circ.cz(j, i); // CX as CZ for Clifford simulator
                circ.s(i);
                circ.cz(j, i);
            } else if k == 2 {
                circ.cz(j, i);
                circ.t(i);
                circ.cz(j, i);
            } else {
                circ.cz(j, i);
                for _ in 0..(1 << (k - 2)) {
                    circ.t(i);
                }
                circ.cz(j, i);
            }
        }
    }
    // Swap qubits
    for i in 0..(n / 2) {
        let a = i;
        let b = n - i - 1;
        circ.cz(a, b); // No native swap, use CZ as placeholder
        circ.cz(b, a);
        circ.cz(a, b);
    }
    circ
}

fn benchmark_bv(n: usize) {
    // predetermined secret: alternating 1,0 pattern
    let secret: Vec<u8> = (0..n).map(|i| if i % 2 == 0 { 1 } else { 0 }).collect();
    println!("BV benchmark ({} qubits)", n);
    println!("secret = {:?}", secret);
    let circ = build_bv(n, &secret);
    let poly = circ.compile();
    let input = vec![0u8; n];
    let start = Instant::now();
    let state = simulate_statevector(&poly, &input);
    let elapsed = start.elapsed();
    let expected = secret
        .iter()
        .enumerate()
        .fold(0usize, |acc, (i, &b)| acc | ((b as usize) << i));
    println!("simulation time: {:?}", elapsed);
    println!("amplitude at expected index {} = {:?}", expected, state[expected]);
    let norm_sq: f64 = state.iter().map(|c| c.norm_sqr()).sum();
    println!("state norm squared = {}", norm_sq);
}

fn simulate_qasm_statevector(qasm_path: &str) {
    println!("\nSimulating QASM file: {}", qasm_path);
    // Use PolyQ's QASM parser
    let circ = match read_qasm_file(qasm_path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Failed to parse QASM: {:?}", e);
            return;
        }
    };
    let n = circ.num_qubits;
    let poly = circ.compile();
    let input = vec![0u8; n];
    let start = Instant::now();
    let state = simulate_statevector(&poly, &input);
    let elapsed = start.elapsed();
    println!("QASM simulation time: {:?}", elapsed);
    // Print norm squared for sanity
    let norm_sq: f64 = state.iter().map(|c| c.norm_sqr()).sum();
    println!("state norm squared = {}", norm_sq);
}

fn main() {
    // cargo run -p benchmark --release
    //cargo build --release --bin benchmark 
    // cargo flamegraph --release --bin benchmark

    simulate_qasm_statevector("dataset/random_circuit_q27_h10_t5_s10_z25_cz10/random_circuit_q27_h10_t5_s10_z25_cz10.qasm");
}