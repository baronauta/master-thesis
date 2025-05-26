function normalmodes(freqs::Vector{Float64}, coups::Vector{Float64})
    # Chain Hamiltonian for the env is Hₑ = c† A c, where A is a three-diagonal matrix:
    # major diagonal formed by frequencies ωn (n = 0...N-1),
    # sub-diagonal and super-diagonal formed by couplings κn (n = 1...N-1).
    # Note: κ0 (coupling tls-env) must be excluded because it doesn't belong with Hₑ.
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
    return D, P
end

function envmodes_occupation(freqs::Vector{Float64}, coups::Vector{Float64}, measN::Matrix{Float64})
    # Any Hermitian (here also real) matrix can be diagonalized by a unitary matrix.
    #   A = P D Pᵀ,
    # with D being a diagonal matrix and P being unitary matrix, Pᵀ=P⁻¹=P†.
    D, P = normalmodes(freqs, coups)
    modes = diag(D)
    # I prefer to see A = U† D U. Define U = Pᵀ. 
    U = transpose(P)
    # In fact, given Hₑ = c† A c with c is vector of annihilation operators,
    # normal modes decomposition of Hₑ is straightforwardly found as
    #   Hₑ = c† (U† D U) c =  t† D t,
    # with t:=Uc, t†:=c†U†. Thus, occupation number of the new modes:
    #   <tₙ† tₙ> = ∑ₖ U[n,k]^2 * <cₖ† cₖ>, with k = 0,...,N-1.

    # Precompute U² values (squared coefficients for each mode)
    U_squared = U .^ 2

    # measN[time, modes], for a fixed time t_idx, I have measN[t_idx,k] = <cₖ† cₖ>.
    # I fictitiously enlarge the chain, see Riva et al.: PRB 108, 195138 (2023), Appendix E.
    if size(measN, 2) < size(U_squared, 2)
        missing_cols = size(U_squared, 2) - size(measN, 2)
        measN = hcat(measN, zeros(size(measN, 1), missing_cols))
    end

    nchain = size(U_squared, 1)  # number of modes
    tt = size(measN, 1)          # number of time steps
    occupations = zeros(Float64, tt, nchain)  # time × mode
    for n = 1:nchain  # loop over modes (columns)
        for t = 1:tt  # loop over time (rows)
            # Compute occupation of mode `n` at time `t`
            # <tₙ†tₙ> = < (∑ₖ Uₙₖ cₖ)† (∑ⱼ Uₙⱼ cⱼ) > = ∑ₖ |Uₙₖ|² < cₖ† cₖ >
            occupations[t, n] = sum(U_squared[n, k] * measN[t, k] for k = 1:nchain)
        end
    end
    return modes, occupations
end


# function write_envmodes(dir::AbstractString, modes::Vector{Float64}, occupations::Matrix{Float64})
#     # Write to a file normal modes 
#     open(joinpath(dir, "envmodes_modes.dat"), "w") do io
#         writedlm(io, modes, ',')
#     end
#     # Write to a file normal modes occupation (time × modes)
#     open(joinpath(dir, "envmodes_occ.dat"), "w") do io
#         writedlm(io, occupations, ',')
#     end
# end

# function read_modes(dir::AbstractString)
#     return readdlm(joinpath(dir, "envmodes_modes.dat"), '\n', Float64)[:, 1]
# end

# function read_occupations(dir::AbstractString)
#     return readdlm(joinpath(dir, "envmodes_occ.dat"), ',', Float64, '\n')
# end
