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


The goal is to reduce the complexity from roughly

```
O(2^(n + H))
```

to something closer to

```
O(2^T)
```

or even

```
O(2^(T/2))
```

for many circuits.

---

# PolyQ Simulator Speedup Techniques

```markdown
# PolyQ Simulator Speedup Techniques

## Baseline Pipeline

Circuit
→ Compile phase polynomial
→ Partition variables (vvars / uvars)
→ Enumerate 2^nv masks
→ Evaluate quadratic exponential sum
→ Accumulate amplitudes
```

---

# 1. Phase Polynomial Compilation (Already Implemented)

Represent the circuit as

```
f(x) = 4(xᵀBx + vᵀx + ε) + rem(x) mod 8
```

Where

```
B : quadratic matrix
v : linear vector
rem : T/S terms
```

### Complexity

```
O(|gates|)
```

### Data Structures

```
B : FixedBitSet adjacency matrix
v : bit vector
rem : list of Z8 terms
```

---

# 2. Clifford / Non-Clifford Variable Partitioning (Already Implemented)

Split variables:

```
internal variables
    ↓
vvars : appear in T/S terms
uvars : purely Clifford
```

Amplitude becomes:

```
A = Σ_v ω^{rem(v)} Σ_u (-1)^{uᵀBu + vᵀCu}
```

### Complexity

```
O(2^nv)
```

instead of

```
O(2^t)
```

where

```
t = n + H
```

---

# 3. Dickson Reduction (Already Implemented)

Evaluate

```
Σ_u (-1)^{uᵀBu + vᵀu}
```

using canonical form.

### Gauss Sum Formula

```
Σ_x (-1)^{xᵀAx} = ± 2^(n - r/2)
```

where

```
r = rank(A)
```

### Complexity

```
O(n^3)
```

---

# 4. Pre-Reduction of Quadratic Matrix (Major Optimization)

Currently Dickson reduction is recomputed every iteration.

Instead precompute canonical form.

### Idea

```
PᵀBP = canonical
```

Then only transform linear vector each iteration.

### Pseudocode

```
function preprocess(B):
    (P, B_canonical) = dickson_reduce_once(B)
    return (P, B_canonical)

function evaluate(v):
    v' = Pᵀ v
    return fast_gauss_sum(B_canonical, v')
```

### Complexity

```
Before: O(2^nv * nu^3)
After : O(nu^3 + 2^nv * nu)
```

---

# 5. Avoid Matrix Cloning

Current code clones matrices each iteration.

```
bu.clone()
vu.clone()
```

Instead reuse memory buffers.

### Pseudocode

```
preallocate vu

for mask in 0..2^nv:
    update vu in place
    evaluate exponential sum
```

### Speedup

```
5× – 20×
```

---

# 6. Cross-Interaction Bitmask Optimization

Precompute interactions between `vvars` and `uvars`.

### Current

```
for ui:
    for vj in cross[ui]:
        toggle
```

### Optimized

Store bitmasks.

```
cross_mask[ui]
```

### Pseudocode

```
toggle = parity(mask & cross_mask[ui])
vu[ui] ^= toggle
```

### Speedup

```
3× – 10×
```

---

# 7. Gray-Code Enumeration

Enumerate masks so only **one bit flips each step**.

### Current

```
mask = binary counter
```

### Optimized

```
mask = gray_code(counter)
```

### Pseudocode

```
prev_mask = 0

for i in 1..2^nv:
    mask = gray(i)
    changed = mask XOR prev_mask
    j = index_of(changed)

    update vu using cross[j]

    prev_mask = mask
```

### Complexity

```
O(nu) per iteration
```

instead of

```
O(nu * nv)
```

---

# 8. Bitmask Representation of T/S Terms

Current evaluation scans each term.

### Instead

Encode terms as masks.

### Pseudocode

```
for term in terms:
    if popcount(mask & term.mask) == term.size:
        phase += term.weight
```

### Speedup

```
5× – 20×
```

---

# 9. Variable Compression (Linear Substitution)

Reduce number of `vvars`.

### Idea

Use quadratic interactions to express variables as XOR combinations.

```
x = P y
```

### Algorithm

1. Build matrix

```
M[i,j] = interaction(u_i, v_j)
```

2. Compute rank.

3. Express dependent variables.

### Pseudocode

```
M = build_uv_matrix()

basis = gaussian_elimination(M)

for v in dependent_columns:
    substitute v = XOR(basis vars)
```

### Effect

```
nv → rank(M)
```

### Speedup

```
10× – 4000× depending on circuit
```

---

# 10. Fast Walsh-Hadamard Transform (FWHT)

Compute sums for all masks simultaneously.

### Transform

```
(a,b) → (a+b , a-b)
```

### Pseudocode

```
function FWHT(vec):
    n = len(vec)
    step = 1
    while step < n:
        for i in range(0,n,2*step):
            for j in range(step):
                u = vec[i+j]
                v = vec[i+j+step]

                vec[i+j] = u + v
                vec[i+j+step] = u - v
        step *= 2
```

### Complexity

```
O(nv * 2^nv)
```

---

# 11. Z₄ Quadratic Form Method

Rewrite phase polynomial:

```
f(x) = aᵀx + 2xᵀBx   (mod 4)
```

Use Z₄ Gaussian elimination.

### Complexity

```
O(n^2)
```

instead of

```
O(n^3)
```

### Speedup

```
5× – 15×
```

---

# 12. Fourier Transform Technique (FTT)

Use Fourier expansion:

```
ω^x = ½[(1+i) + (1-i)(-1)^x]
```

Each T gate becomes **two Clifford branches**.

### Complexity

```
O(2^(T/2))
```

instead of

```
O(2^T)
```

---

# 13. Tensor Network Fallback

For circuits with low interaction graph treewidth.

### Algorithm

```
build interaction graph
estimate treewidth

if treewidth < nv:
    use tensor contraction
```

### Complexity

```
O(2^treewidth)
```

---

# 14. Parallelization

Amplitude computations are independent.

### Strategy

```
parallel over output states
parallel over mask enumeration
```

### Pseudocode

```
parallel for y in 0..2^n:
    amp[y] = compute_amplitude(y)
```

Use:

```
Rayon / OpenMP / GPU kernels
```

---

# 15. Hybrid Simulation Strategy

Choose algorithm dynamically.

### Pseudocode

```
if nv < 25:
    phase polynomial simulator

elif treewidth small:
    tensor contraction

else:
    stabilizer decomposition
```

---

# Final Complexity Landscape

| Method               | Complexity     |
| -------------------- | -------------- |
| Naive statevector    | O(2^n)         |
| Phase polynomial     | O(2^T)         |
| Compressed variables | O(2^(rank))    |
| FTT                  | O(2^(T/2))     |
| Tensor contraction   | O(2^treewidth) |

---

# Expected Combined Speedups

| Optimization          | Speedup      |
| --------------------- | ------------ |
| Pre-Dickson reduction | 100× – 1000× |
| Variable compression  | 10× – 4000×  |
| FWHT                  | 5× – 20×     |
| Z4 method             | 5× – 15×     |
| Gray code + bitmasks  | 10× – 50×    |

Potential total improvement:

```
10^4 – 10^6
```

over naive simulation.

```

---

## Final remark

Your current simulator **already implements the most important ideas**:

- Phase polynomial representation  
- Clifford/non-Clifford split  
- Dickson reduction  

Adding the optimizations above would turn PolyQ into something **very close to a research-grade Clifford+T simulator**.

---

If you want, I can also give you **one final section that is extremely useful for your BTP thesis**:

```

"Which of these optimizations are publishable / novel vs standard?"

```

because a few of them could actually form a **small research contribution.**
```
