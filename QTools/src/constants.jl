"""
Pauli basis for 2×2 Hermitian operators.

Defines the set `{σ₀, σ₁, σ₂, σ₃}` where:

- `σ₀` = identity matrix,
- `σ₁` = Pauli-X,
- `σ₂` = Pauli-Y,
- `σ₃` = Pauli-Z.

These form an orthogonal (but not orthonormal) basis under the Hilbert-Schmidt inner product:
`⟨A, B⟩ = Tr(A†B)`.

To obtain an orthonormal basis, normalize each σₖ as `σₖ / √2`.
"""
const paulimatrices = [
    [1.0 0.0; 0.0 1.0],    # σ₀
    [0.0 1.0; 1.0 0.0],    # σ₁
    [0.0 -im; im 0.0],     # σ₂
    [1.0 0.0; 0.0 -1.0],   # σ₃
]


"""
Canonical operator basis for 2×2 Hermitian operators.

Defines the set of 2×2 matrices `T₀₀`, `T₀₁`, `T₁₀`, `T₁₁`, where each `Tᵢⱼ` corresponds to the operator
|i⟩⟨j| in Dirac notation.

These matrices:
- form a complete basis for 2×2 complex matrices,
- are orthonormal under the Hilbert-Schmidt inner product ⟨A, B⟩ = Tr(A†B).
"""
const canonical_op_basis = [
     [1.0 0.0; 0.0 0.0],
     [0.0 1.0; 0.0 0.0],
     [0.0 0.0; 1.0 0.0],
     [0.0 0.0; 0.0 1.0],
]


"""
Permutation matrix to reshuffle A-form into B-form matrix.

Reference:
V. Jagadish and F. Petruccione, “An invitation to quantum channels”.
"""
const A2Bmatrix = [
#! format: off
#    1  2  3  4   5  6  7  8   9 10 11 12  13 14 15 16
     1  0  0  0   0  0  0  0   0  0  0  0   0  0  0  0   # 1
     0  0  0  0   0  0  0  0   1  0  0  0   0  0  0  0   # 2
     0  0  1  0   0  0  0  0   0  0  0  0   0  0  0  0   # 3
     0  0  0  0   0  0  0  0   0  0  1  0   0  0  0  0   # 4
     0  0  0  0   1  0  0  0   0  0  0  0   0  0  0  0   # 5
     0  0  0  0   0  0  0  0   0  0  0  0   1  0  0  0   # 6
     0  0  0  0   0  0  1  0   0  0  0  0   0  0  0  0   # 7
     0  0  0  0   0  0  0  0   0  0  0  0   0  0  1  0   # 8
     0  1  0  0   0  0  0  0   0  0  0  0   0  0  0  0   # 9
     0  0  0  0   0  0  0  0   0  1  0  0   0  0  0  0   # 10
     0  0  0  1   0  0  0  0   0  0  0  0   0  0  0  0   # 11
     0  0  0  0   0  0  0  0   0  0  0  1   0  0  0  0   # 12
     0  0  0  0   0  1  0  0   0  0  0  0   0  0  0  0   # 13
     0  0  0  0   0  0  0  0   0  0  0  0   0  1  0  0   # 14
     0  0  0  0   0  0  0  1   0  0  0  0   0  0  0  0   # 15
     0  0  0  0   0  0  0  0   0  0  0  0   0  0  0  1   # 16
#! format: on
]

