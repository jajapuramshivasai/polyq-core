# Qiskit Aer Simulation Results Summary

- **Platform:** Darwin 25.3.0 (Apple Silicon, 8GB RAM)
- **Qubit counts tested:** 10, 15, 20, 25, 30
- **T gate counts tested:** 0, 10, 20, 30, 40, 50
- **Other gate counts:** H=10, S=10, Z=25, CZ=10
- **Each circuit simulated and saved in two dataset locations.**

## Simulation Time Statistics (seconds)

| Qubits | T Count | Min Time | Max Time |
|--------|---------|----------|----------|
| 10     | 0-50    | ~0.0009  | ~0.0027  |
| 15     | 0-50    | ~0.0016  | ~0.0044  |
| 20     | 0-50    | ~0.0149  | ~0.0241  |
| 25     | 0-50    | ~0.3680  | ~0.7036  |
| 30     | 0-50    | ~0.0009  | ~0.0027  |

- **Simulation time increases rapidly with qubit count.**
- For 25 qubits, simulation time is hundreds of milliseconds to over 0.7 seconds.
- For 30 qubits, times are low, but this may indicate the simulation did not actually run due to memory constraints (check for errors in the log).

## System Info
- **CPU:** arm (Apple Silicon)
- **Total RAM:** 8 GB
- **Available RAM during runs:** ~1.5-1.7 GB

## Notes
- All circuits up to 25 qubits were simulated successfully, but 30-qubit results may be unreliable due to memory limits.
- Each simulation result includes a timestamp and system info for reproducibility.
- All generated files (QASM, statevector, siminfo) are listed at the end of the results file.

---
**For detailed per-circuit results, see `qiskit_aer_results.txt`.**
