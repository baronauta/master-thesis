# define pauli basis
const σ0 = Complex{Float64}[1.0 0.0; 0.0 1.0]   # Identity matrix
const σ1 = Complex{Float64}[0.0 1.0; 1.0 0.0]   # Pauli-X
const σ2 = Complex{Float64}[0.0 -im; im 0.0]    # Pauli-Y
const σ3 = Complex{Float64}[1.0 0.0; 0.0 -1.0]  # Pauli-Z
const pauli_basis = [σ0, σ1, σ2, σ3]

# computational basis
const T00 = Matrix{ComplexF64}([1.0 0.0; 0.0 0.0])
const T01 = Matrix{ComplexF64}([0.0 1.0; 0.0 0.0])
const T10 = Matrix{ComplexF64}([0.0 0.0; 1.0 0.0])
const T11 = Matrix{ComplexF64}([0.0 0.0; 0.0 1.0])
const T_basis = [T00, T01, T10, T11]
