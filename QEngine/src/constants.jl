"""
Pauli basis for 2×2 Hermitian operators.

Defines the set of Pauli matrices {σ₀, σ₁, σ₂, σ₃}, where:
- σ₀ is the 2×2 identity matrix,
- σ₁ is the Pauli-X matrix,
- σ₂ is the Pauli-Y matrix,
- σ₃ is the Pauli-Z matrix.

These matrices form an orthogonal basis for 2×2 Hermitian matrices.
They are not orthonormal with reference to the Hilbert-Schmidt inner product 
defined as <A,B> = Tr[A†B]. Their normalized version is obtained considering
σₖ → σₖ/√2.
"""
const σ0 = Matrix{ComplexF64}([1.0 0.0; 0.0 1.0])   # Identity matrix
const σ1 = Matrix{ComplexF64}([0.0 1.0; 1.0 0.0])   # Pauli-X
const σ2 = Matrix{ComplexF64}([0.0 -im; im 0.0])    # Pauli-Y
const σ3 = Matrix{ComplexF64}([1.0 0.0; 0.0 -1.0])  # Pauli-Z
const paulimatrices = [σ0, σ1, σ2, σ3]

"""
Canonical operator basis for a two-level quantum system.

Defines the set of 2×2 matrices {T₀₀, T₀₁, T₁₀, T₁₁}, where each Tᵢⱼ corresponds to the operator
|i⟩⟨j| in Dirac notation. These matrices form a complete basis for the space of 2×2 matrices.
These matrices provide an orthonormal basis with reference to the Hilbert-Schmidt inner product 
defined as <A,B> = Tr[A†B].

Explicitly:
T₀₀ = |0⟩⟨0| = [1 0; 0 0]  
T₀₁ = |0⟩⟨1| = [0 1; 0 0]  
T₁₀ = |1⟩⟨0| = [0 0; 1 0]  
T₁₁ = |1⟩⟨1| = [0 0; 0 1]

Grouped in the list `canonical_op_basis = [T₀₀, T₀₁, T₁₀, T₁₁]`.
"""
const T00 = Matrix{ComplexF64}([1.0 0.0; 0.0 0.0])
const T01 = Matrix{ComplexF64}([0.0 1.0; 0.0 0.0])
const T10 = Matrix{ComplexF64}([0.0 0.0; 1.0 0.0])
const T11 = Matrix{ComplexF64}([0.0 0.0; 0.0 1.0])
const canonical_op_basis = [T00, T01, T10, T11]

"""
Permutation matrix to reshuffle a A-form matrix into a B-form matrix.
For reference see V. Jagadish and F. Petruccione, “An invitation to quantum channels”.
"""
#! format: off
const A2Bmatrix = Matrix{ComplexF64}([
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
])
#! format: on
