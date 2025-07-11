using LinearAlgebra

function diagonalize_tridiagonal(freqs::Vector{Float64}, coups::Vector{Float64})
    # Chain Hamiltonian for the environment: H = c† A c
    #   - c vector of chain operators
    #   - A three-diagonal matrix
    # Matrix A
    #   - major diagonal: frequencies ωn
    #   - sub-diagonal and super-diagonal: couplings κn
    #     (exclude κ0 that is the coupling between the qubit and the env)
    k = coups[2:end]
    A = Tridiagonal(k, freqs, k)
    # Any Hermitian (here also real) matrix can be diagonalized by a unitary matrix.
    #   A = P D Pᵀ,
    # with D being a diagonal matrix and P being unitary matrix, Pᵀ=P⁻¹=P†.
    # Eigensolver:
    # F.values contains the eigenvalues (diagonal entries of D);
    # F.vectors contains the eigenvectors (columns of P).
    F = eigen(A)
    D = Diagonal(F.values)
    P = F.vectors

    modes = diag(D)
    U = Array(transpose(P))    # I prefer A = U† D U, so define U = Pᵀ

    return modes, U
end

"""
Compute the occupation of environment modes.

# Arguments
- `U::Matrix{Float64}`: Mode transformation matrix (n_modes × n_sites).
- `cdagc::Matrix{ComplexF64}`: Correlation matrix ⟨cₖ† cⱼ⟩ in the site basis.
"""
function envmodes_occupation(U::Matrix{Float64}, cdagc::Matrix{ComplexF64})
    # Number of chain sites
    nchain = size(cdagc, 1)
    occupations = zeros(ComplexF64, nchain)  # modes
    for n = 1:nchain  # loop over modes (columns)
        # Compute occupation of mode 'n'
        # <bₙ†bₙ> = < (∑ₖ Uₙₖ cₖ)† (∑ⱼ Uₙⱼ cⱼ) > = ∑ₖ ∑ⱼ Uₙₖ* Uₙⱼ < cₖ† cⱼ >
        occupations[n] =
            sum(conj(U[n, k]) * U[n, j] * cdagc[k, j] for k = 1:nchain for j = 1:nchain)

    end
    return occupations
end