#![allow(warnings)]


use PolyQ::amplitude_clifford_t_accel;
use PolyQ::sim::{self, Circuit, simulate_statevector, simulate_statevector_parallel};
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
    // Apply Hadamard and controlled phase rotations
    for i in 0..n {
        circ.h(i);
        for j in (i + 1)..n {
            // Controlled RZ: phase = PI / 2^(j-i)
            let k = j - i + 1;
            let phase = 1u16 << (16 - k); // PHASE_BITS = 16
            circ.rz(i, phase);
        }
    }
    // Swap qubits to reverse order
    for i in 0..(n / 2) {
        let a = i;
        let b = n - i - 1;
        // If swap gate is not available, use CZ as placeholder
        circ.cz(a, b);
        circ.cz(b, a);
        circ.cz(a, b);
    }
    circ
}
fn benchmark_qft_statevector(n: usize) {
    println!("QFT benchmark ({} qubits)", n);
    let circ = build_qft(n);
    let poly = circ.compile();
    let input = vec![0u8; n];
    let start = Instant::now();
    let state = simulate_statevector(&poly, &input);
    let elapsed = start.elapsed();
    println!("QFT simulation time: {:?}", elapsed);
    // Print norm squared for sanity
    let norm_sq: f64 = state.iter().map(|c| c.norm_sqr()).sum();
    println!("state norm squared = {}", norm_sq);
}
fn benchmark_qft_amplitude(n: usize) {
    println!("QFT benchmark ({} qubits)", n);
    let circ = build_qft(n);
    let poly = circ.compile();
    let input = vec![0u8; n];
    let target_index = 0; // amplitude of |0...0⟩
    let start = Instant::now();
    let amp = amplitude_clifford_t_accel(&poly, &input, target_index);
    let elapsed = start.elapsed();
    println!("Amplitude at index {} = {:?}", target_index, amp);
    println!("QFT amplitude simulation time: {:?}", elapsed);
}
fn benchmark_bv_statevector(n: usize) {
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


fn benchmark_bv_amplitude(n: usize) {
    // predetermined secret: alternating 1,0 pattern
    let secret: Vec<u8> = (0..n).map(|i| if i % 2 == 0 { 1 } else { 0 }).collect();
    println!("BV benchmark ({} qubits)", n);
    println!("secret = {:?}", secret);
    let circ = build_bv(n, &secret);
    let poly = circ.compile();
    let input = vec![0u8; n];
    let expected = secret
        .iter()
        .enumerate()
        .fold(0usize, |acc, (i, &b)| acc | ((b as usize) << i));
    let start = Instant::now();
    let amp = amplitude_clifford_t_accel(&poly, &input, expected);
    let elapsed = start.elapsed();
    println!("amplitude at expected index {} = {:?}", expected, amp);
    println!("amplitude simulation time: {:?}", elapsed);
}

fn simulate_qasm_amplitude(qasm_path: &str) {
    //benchmark clifford+t amplitude simulation for a single amplitude
    println!("\nSimulating amplitude for QASM file: {}", qasm_path);
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
    // Choose a target index (e.g., all-zeros)
    let target_index = 0;
    let start = Instant::now();
    let amp = amplitude_clifford_t_accel(&poly, &input, target_index);
    let elapsed = start.elapsed();
    println!("Amplitude at index {} = {:?}", target_index, amp);
    println!("Amplitude simulation time: {:?}", elapsed);

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

fn simulate_qasm_statevector_parallel(qasm_path: &str) {
    println!("\nSimulating QASM file in parallel: {}", qasm_path);
    // Print number of Rayon threads
    let num_threads = rayon::current_num_threads();
    println!("Rayon num_threads = {}", num_threads);
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
    let state = simulate_statevector_parallel(&poly, &input);
    let elapsed = start.elapsed();
    println!("Parallel QASM simulation time: {:?}", elapsed);
    // Print norm squared for sanity
    let norm_sq: f64 = state.iter().map(|c| c.norm_sqr()).sum();
    println!("state norm squared = {}", norm_sq);
}


fn main() {

    //pkill -f cargo

    // cargo run --release --bin benchmark 
    // cargo flamegraph --release --bin benchmark
    // samply record cargo run --release
    // RUSTFLAGS="-C force-frame-pointers=yes" samply record cargo run --release

    // benchmark_bv_amplitude(20);
    // benchmark_bv_statevector(18);

    // benchmark_qft_amplitude(18);
    // benchmark_qft_statevector(20);


    simulate_qasm_amplitude("experiments/test/test.qasm");
    simulate_qasm_statevector("experiments/test/test.qasm");
    simulate_qasm_statevector_parallel("experiments/test/test.qasm");
   
}
