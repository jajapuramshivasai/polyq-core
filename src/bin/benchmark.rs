use PolyQ::sim::{Circuit, simulate_statevector};
use std::time::Instant;

/// Build a Bernstein–Vazirani circuit for a secret bitstring.
fn build_bv(n: usize, secret: &[u8]) -> Circuit {
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

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.get(1).map(|s| s == "csv").unwrap_or(false) {
        // open file in project root
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create("bv_times.csv").expect("failed to create output file");
        writeln!(file, "qubits,state_secs,amp_secs,amp_re,amp_im,amp_norm").unwrap();
        for n in 4..=15 {
            let secret: Vec<u8> = (0..n).map(|i| if i % 2 == 0 { 1 } else { 0 }).collect();
            let circ = build_bv(n, &secret);
            let poly = circ.compile();
            let input = vec![0u8; n];
            // time full statevector
            let start_state = Instant::now();
            let state = simulate_statevector(&poly, &input);
            let elapsed_state = start_state.elapsed();
            let state_secs = elapsed_state.as_secs_f64();
            // compute expected index
            let expected = secret
                .iter()
                .enumerate()
                .fold(0usize, |acc, (i, &b)| acc | ((b as usize) << i));
            // time amplitude routine
            let start_amp = Instant::now();
            let amp = PolyQ::sim::amplitude_clifford_t_accel(&poly, &input, expected);
            let elapsed_amp = start_amp.elapsed();
            let amp_secs = elapsed_amp.as_secs_f64();
            let norm = amp.norm();
            writeln!(file, "{},{},{},{},{},{}", n, state_secs, amp_secs, amp.re, amp.im, norm).unwrap();
        }
        println!("wrote bv_times.csv to project root");
        return;
    }

    let n: usize = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

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

    // compute little‑endian index of secret
    let expected = secret
        .iter()
        .enumerate()
        .fold(0usize, |acc, (i, &b)| acc | ((b as usize) << i));

    println!("simulation time: {:?}", elapsed);
    println!("amplitude at expected index {} = {:?}", expected, state[expected]);
    let norm_sq: f64 = state.iter().map(|c| c.norm_sqr()).sum();
    println!("state norm squared = {}", norm_sq);
}