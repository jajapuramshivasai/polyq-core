apply the following optimisations

Here is the complete summary of the four major classical simulation optimizations 

### 1. Dickson Block Reduction over $\mathbb{Z}_4$ (The Mathematical Core)
This foundational optimization replaces exponential-time statevector array multiplications with polynomial-time structural graph analysis.

* **Concept:** A quantum circuit (Clifford+T) is compiled into a characteristic Boolean polynomial acting over a $\mathbb{Z}_4$ phase ring. Because exact diagonalization is impossible for alternating matrices without $S$ gates, we use Dickson's theorem (a symplectic Gauss-Jordan equivalent) to decouple the variable interactions into isolated $2 \times 2$ blocks.
* **Procedure:** 1. Construct the alternating Boolean matrix $B$ to track quadratic interactions ($2x_i x_j$) and the linear vector $v$ to track phase weights ($1x_i, 2x_i$).
    2. Iteratively locate a coupled pivot pair, move it to the block diagonal, and use symmetric row/column XOR operations to clear all other edges in the graph.
    
    3. Apply these identical operations to update the linear phase vector $v$ using $\mathbb{Z}_4$ arithmetic.
    4. Calculate the total amplitude by iterating over the 4 states of each isolated block and checking kernel variables for destructive interference (which forces the sum to 0).
* **Complexity Reduction:** Collapses the evaluation of a single transition amplitude from $O(2^h)$ exponential time down to $O(m^3)$ polynomial time, where $m$ is the number of internal variables.

---

### 2. The Transfer Matrix (Algebraic Statevector Batching)
Generating the full statevector naively would mean re-running the $O(m^3)$ matrix reduction for all $2^n$ possible output bitstrings. This optimization algebraically bypasses that loop entirely.

* **Concept:** The adjacency matrix $B$ is entirely invariant to output boundary conditions; only the linear phase vector $v$ shifts when the target bitstring changes. Because the sequence of operations in the Dickson reduction is a purely linear transformation, its effects can be distributed and precomputed.
* **Procedure:**
    1. Run the $O(m^3)$ matrix reduction exactly once on the base state.
    2. Precompute the isolated transformation of the linear vector for every individual output bit. Assemble these resulting shift vectors into an $m \times n$ Transfer Matrix.
    3. Iterate through the $2^n$ output space using a Gray code, so only one target bit flips per step.
    
    4. In the hot loop, skip the matrix reduction completely. Simply update your current canonical vector by adding the $k$-th column of the Transfer Matrix via a single, fast array addition.
* **Complexity Reduction:** Plummets the statevector generation complexity from a devastating $O(2^n \cdot m^3)$ down to an initial $O(m^3)$ compile step followed by a highly optimized $O(2^n \cdot m)$ array addition loop.

---

### 3. Pure Integer Tracking (The Arithmetical Optimization)
Evaluating the decoupled $2 \times 2$ blocks millions of times using `Complex64` structs throttles the CPU and ruins hardware vectorization.

* **Concept:** The complex exponential sum of any $2 \times 2$ Dickson block over $\mathbb{Z}_4$ is heavily constrained. Every non-zero block definitively evaluates to a magnitude that is a power of $\sqrt{2}$ and a phase angle that is a multiple of $\pi/4$.
* **Procedure:** 1. Strip all floating-point math from the inner loops.
    2. Initialize two integers: `mag_half_powers` (to track magnitude exponents) and `phase_pi4` (to track rotational phase modulo 8).
    3. Read the coefficients $(v_1, v_2)$ of a decoupled block and use a fast bitwise match statement to instantly map the state to basic integer additions.
    4. Perform the floating-point conversion strictly once at the end of the calculation by pulling from a pre-allocated array of the 8 roots of unity.
* **Complexity Reduction:** Constant-factor overhead drops massively. The inner-loop block evaluations execute in roughly ~5 CPU cycles using pure integer arithmetic, perfectly enabling LLVM SIMD auto-vectorization.

---

### 4. Static vs. Dynamic Graph Partitioning (Constant Folding)
In circuits with localized entanglement, many internal variables never interact with the final output qubits, meaning they are re-evaluated redundantly in the hot loop.

* **Concept:** Treat the polynomial matrix as a dependency graph. Sub-graphs that do not share edges with the output boundary are mathematically invariant across the entire $2^n$ statevector generation space. 
* **Procedure:**
    1. During the pre-computation phase, scan matrix $B$ and partition the canonical variables into a "Static Set" (no output edges) and a "Dynamic Set" (shares edges with outputs).
    2. Run the integer summation logic on the Static Set exactly once to compute a baseline `base_mag` and `base_phase`.
    3. Inside the $2^n$ hot loop, load the baseline values and execute the integer summation strictly on the Dynamic Set.
* **Complexity Reduction:** Shrinks the innermost evaluation loop from length $O(m)$ down to $O(k)$, where $k$ is the number of dynamic variables. For heavily but locally entangled circuits, this skips the bulk of the circuit for every statevector step.

---