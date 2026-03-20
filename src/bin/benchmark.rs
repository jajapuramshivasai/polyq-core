#![allow(dead_code)]
#![allow(unused_mut)]

use PolyQ::amplitude_clifford_t_accel;
use PolyQ::qc::read_qasm_file;
use PolyQ::sim::{Circuit, simulate_statevector};
use std::time::Instant;

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
            let phase = 1u32 << (16 - k); // PHASE_BITS = 16
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
    println!(
        "amplitude at expected index {} = {:?}",
        expected, state[expected]
    );
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

fn main() {
    // cargo run --release --bin benchmark
    // cargo flamegraph --release --bin benchmark
    // samply record cargo run --release
    // RUSTFLAGS="-C force-frame-pointers=yes" samply record cargo run --release

    benchmark_bv_amplitude(20);
    benchmark_bv_statevector(18);

    // benchmark_qft_amplitude(18);
    // benchmark_qft_statevector(20);

    // simulate_qasm_statevector("dataset/random_circuit_q27_h10_t5_s10_z25_cz10/random_circuit_q27_h10_t5_s10_z25_cz10.qasm");
    //dataset/random_circuit_q27_h10_t0_s10_z25_cz10
    // simulate_qasm_amplitude("dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm");
    // simulate_qasm_statevector("dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm");

    // simulate_qasm_amplitude("experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm");
    // simulate_qasm_statevector("experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm");
}

/*

BV benchmark (18 qubits)
secret = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
simulation time: 655.003416ms
amplitude at expected index 87381 = Complex { re: 1.0, im: 0.0 }
state norm squared = 1

Simulating amplitude for QASM file: dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm
Amplitude at index 0 = Complex { re: 0.125, im: 0.0 }
Amplitude simulation time: 2.792µs

Simulating QASM file: dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm
QASM simulation time: 6.537444458s
state norm squared = 1

Simulating amplitude for QASM file: experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm
Amplitude at index 0 = Complex { re: 0.03125, im: 0.0 }
Amplitude simulation time: 13.125µs

Simulating QASM file: experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm
QASM simulation time: 7.208341875s
state norm squared = 1



BV benchmark (18 qubits)
secret = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
simulation time: 447.444083ms
amplitude at expected index 87381 = Complex { re: 1.0, im: 0.0 }
state norm squared = 1

Simulating amplitude for QASM file: dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm
Amplitude at index 0 = Complex { re: 0.125, im: 0.0 }
Amplitude simulation time: 2.584µs

Simulating QASM file: dataset/random_circuit_q27_h10_t0_s10_z25_cz10/random_circuit_q27_h10_t0_s10_z25_cz10.qasm
QASM simulation time: 5.915788s
state norm squared = 1

Simulating amplitude for QASM file: experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm
Amplitude at index 0 = Complex { re: 0.03125, im: 0.0 }
Amplitude simulation time: 28.25µs

Simulating QASM file: experiments/dataset/random_circuit_q27_h10_t70_s10_z25_cz10/random_circuit_q27_h10_t70_s10_z25_cz10.qasm
QASM simulation time: 6.535101875s
state norm squared = 1

*/
