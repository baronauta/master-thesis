# Pauli basis
const σ0 = Complex{Float64}[1.0 0.0; 0.0 1.0]   # Identity matrix
const σ1 = Complex{Float64}[0.0 1.0; 1.0 0.0]   # Pauli-X
const σ2 = Complex{Float64}[0.0 -im; im 0.0]    # Pauli-Y
const σ3 = Complex{Float64}[1.0 0.0; 0.0 -1.0]  # Pauli-Z
const pauli_basis = [σ0, σ1, σ2, σ3]

# Computational basis
const T00 = Matrix{ComplexF64}([1.0 0.0; 0.0 0.0])
const T01 = Matrix{ComplexF64}([0.0 1.0; 0.0 0.0])
const T10 = Matrix{ComplexF64}([0.0 0.0; 1.0 0.0])
const T11 = Matrix{ComplexF64}([0.0 0.0; 0.0 1.0])
const T_basis = [T00, T01, T10, T11]

# Permutation matrix to reshuffle a A-form matrix into a B-form matrix.
# For reference see V. Jagadish and F. Petruccione, “An invitation to quantum channels”.
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
