# Quantum Circuit Simulator via Boolean Phase Polynomials

## 1. Working Mechanism
This simulator evaluates quantum circuit amplitudes without using dense matrix-vector multiplications in a $2^n$-dimensional Hilbert space. Instead, it uses a phase polynomial formalism:

* **Polynomial Mapping:** The quantum circuit is parsed into a Boolean characteristic polynomial over \mathbb{F}_2. Every Hadamard gate introduces a new variable, creating interference paths.
* **$\mathbb{Z}_8$ Extension for Clifford+T:** To support a universal gate set, the polynomial is extended to evaluate phases modulo 8. Gates are assigned specific integer weights (e.g., T=1, S=2, Z=4, CZ=4, H=4).
* **Boundary Conditions:** To calculate a specific transition amplitude $\langle y | U | x \rangle$, the input variables are fixed to the state $x$, and the output variables are fixed to the state $y$. 
* **Dickson's Theorem (Clifford Reduction):** The weight-4 terms (Clifford operations) form a quadratic Boolean function. The simulator converts the associated symplectic matrix into a block-diagonal canonical form using Dickson's Theorem. This allows the Hamming weight, and thus the exponential sum of the Clifford part, to be calculated algebraically.
* **Amplitude Evaluation:** For circuits with non-Clifford gates (T, S), the simulator isolates the variables involved in these gates. It iterates through all possible binary assignments of these specific variables, applies the phase factor, and solves the remaining purely quadratic (Clifford) problem in polynomial time using the Dickson reduction.

## 2. Complexity Analysis

### Time Complexity: $\mathcal{O}(2^v \cdot h^3)$
* $h$ is the number of Hadamard gates.
* $v$ is the number of internal variables involved in non-Clifford operations (T and S gates).
* **Explanation:** Reducing the symplectic matrix via Dickson's Theorem takes $\mathcal{O}(m^3)$ time, where $m \le h$. For pure Clifford circuits ($v = 0$), the entire simulation finishes in polynomial time ($\mathcal{O}(h^3)$). The exponential overhead ($\mathcal{O}(2^v)$) is strictly limited to evaluating the non-Clifford terms that break the quadratic structure.

### Space Complexity: $\mathcal{O}(g + h^2)$
* $g$ is the total number of gates.
* $h$ is the number of Hadamard gates.
* **Explanation:** The simulator stores the polynomial terms parsed from the $g$ gates. The largest data structure maintained in memory is the $h \times h$ symplectic adjacency matrix (`b4` represented as a `FixedBitSet` array). By avoiding the $2^n$ dense statevector array, the simulator entirely bypasses the memory limitations that cause standard simulators (like Qiskit-Aer) to crash at high qubit counts.



# TODO

 - use G(X) caching

 - applications : peaked quantum circuits , t-count reduction

 - parallize z8 loop with multithreading and spread output statvector in gray code to other workers in cluster with rsmpi

 - change haspmap to vec of fixed length , try simd and profile the code