# ─────────────────────────────────────────────────────────────
# Quantum dynamics
#
# Let S(H) be the set of physical state on Hilbert space H.
# - Dynamical map: 
#   Φ(t): S(H) → S(H); ρ(t) = Φ(t)[ρ(0)].
# - Generator:
#   L(t): S(H) → S(H); d/dt ρ(t) = L(t)[ρ(t)].
#
# ─────────────────────────────────────────────────────────────

function horrendus_filtering(qmap::Matrix{ComplexF64})
    eps = 1e-10

    A = [
        1 1 1 1
        1 1 0 0
        0 0 1 1
        0 0 1 1
    ]
    B = [
        0 0 0 0
        0 0 eps eps
        eps eps 0 0
        eps eps 0 0
    ]

    return A .* qmap .+ B
end

"""
    qmaptomography(dirdata::String)

Performs quantum process tomography from measurement data in `dirdata`. 
Returns a named tuple with the quantum map (`qmap`) and associated `time` vector.
"""
function qmaptomography(dirdata::String)
    # Load measurements and extract results
    meas_files = ["Up", "Dn", "+", "i"]
    meas_data = Dict(
        name => get_measurements(dirdata * "/measurements_$name.dat", "densitymatrix")
        for name in meas_files
    )
    Up = meas_data["Up"].result
    Down = meas_data["Dn"].result
    Plus = meas_data["+"].result
    Trans = meas_data["i"].result
    # Extract time (assumes all share the same timeline)
    time = meas_data["Up"].time

    # Compute evolved Pauli basis
    evolved_basis =
        [Up .+ Down, 2 .* Plus .- Up .- Down, 2 .* Trans .- Up .- Down, Up .- Down]

    # Compute the quantum map as a vector of matrices
    qmap = map(eachindex(Up)) do t
        [0.5 * tr(σi' * evolved_basis[j][t]) for σi in paulimatrices, j = 1:4]
    end
    qmap = [reshape(q, 4, 4) for q in qmap]

    # qmap = horrendus_filtering.(qmap)

    # Return results as named tuple
    return (qmap = qmap, time = time)
end

"""
    qmap_to_gen(qmap::Vector{Matrix{ComplexF64}}, time::Vector{Float64})

Computes the time-local generator `gen` from a time-dependent quantum map `qmap`, 
using finite difference methods. Returns a named tuple with `gen` and `time`.
"""
function qmap_to_gen(qmap::Vector{Matrix{ComplexF64}}, time::Vector{Float64})
    # The generator L(t) is related to the dynamical map Φ(t)
    # by means of L(t) = dΦ/dt Φ^(-1).
    N = length(time)
    gen = Vector{Matrix{ComplexF64}}(undef, N)
    # Integration step (assumes constant timestep)
    dt = time[2] - time[1]

    # First point: forward difference
    derivative = (qmap[2] - qmap[1]) / dt
    gen[1] = derivative * inv(qmap[1])

    # Interior points: centered difference
    for t = 2:N-1
        derivative = (qmap[t+1] - qmap[t-1]) / (2 * dt)
        gen[t] = derivative * inv(qmap[t])
    end

    # Last point: backward difference
    derivative = (qmap[N] - qmap[N-1]) / dt
    gen[N] = derivative * inv(qmap[N])

    return (gen = gen, time = time)
end

"""
    changebasis(Φ::Matrix{ComplexF64}, from_basis::Vector{Matrix{ComplexF64}}, to_basis::Vector{Matrix{ComplexF64}})

Change the representation of a matrix between orthonormal operator bases.

Transforms the matrix `Φ`, which is assumed to be expressed in the operator basis `from_basis`,
into its representation in the target basis `to_basis`. Both bases must be orthonormal with
respect to the Hilbert-Schmidt inner product ⟨A, B⟩ = Tr[A†B] (check it by yourself).

The transformation is performed via a change-of-basis matrix `Λ` with entries Λₖₗ = Tr[(fromₖ)† * toₗ].

!!! warning
    Numerical rounding errors may occur during the transformation, especially when converting
    between significantly different bases (e.g., Pauli ↔ canonical). Small components (|x| < 1e-12)
    are rounded to zero to reduce numerical noise.
"""
function changebasis(
    Φ::Matrix{ComplexF64},
    from_basis::Vector{Matrix{ComplexF64}},
    to_basis::Vector{Matrix{ComplexF64}},
)
    # Hilbert-Schmidt inner product: <A,B> = Tr[A†B]
    Λ = [tr(A' * B) for A in from_basis, B in to_basis]
    Φnew = Λ' * Φ * Λ
    # I noticed that when transforming between Pauli basis to canonical basis, this computation
    # suffers from rounding errors.
    return Φnew = map(x -> abs(real(x)) < 1e-12 && abs(imag(x)) < 1e-12 ? 0.0 : x, Φnew)
end

"""
    changebasis(Φ::Matrix{ComplexF64})

Change the representation of a matrix between orthonormal operator bases:
- from the Pauli basis {σ₀/√2, σ₁/√2, σ₂/√2, σ₃/√2};
- to the canonical basis {|0⟩⟨0|, |0⟩⟨1|, |1⟩⟨0|, |1⟩⟨1|}.

Both bases are orthonormal with respect to the Hilbert-Schmidt inner product ⟨A, B⟩ = Tr[A†B].
"""
function changebasis(Φ::Matrix{ComplexF64})
    # Transformation matrix is Λ = 1/√2 ⋅ M, where M is defined as follows.
    # Given Φ in Pauli basis its representation in the canonical basis is
    # Λ† Φ Λ = 1/2 M† Φ M.
    # I do this to avoid numerical roundings of √2.
    M = Matrix{ComplexF64}([1 0 0 1; 0 1 1 0; 0 im -im 0; 1 0 0 -1])
    return 1 / 2 * M' * Φ * M
end

"""
    krausdecomposition(M::Matrix{ComplexF64}

Computes the Kraus decomposition of a Hermitian superoperator.  
Returns a tuple `(eigvals, krausoperators)`, where `eigvals` are the eigenvalues 
and `krausoperators` are the corresponding Kraus operators.
"""
function krausdecomposition(M::Matrix{ComplexF64})
    # M is Hermitian ⇒ admits spectral decomposition:
    # i.e. there exists a unitary Q and real eigenvalues λₖ such that
    #       M = Q Λ Q† = ∑ₖ λₖ q(k) q(k)†,
    # where the k-th column of Q is the eigenvector q(k) of M,
    # and λₖ are the corresponding eigevalues.
    if !(M ≈ M')
        throw(ArgumentError("Input matrix is not Hermitian"))
    end
    E = eigen(M)
    # E.values contains the eigenvalues;
    # E.vectors contains the eigenvectors (stored as columns).
    eigvals = E.values
    Q = E.vectors

    # Kraus decomposition of a superoperator S: S(H) → S(H):
    #       S(ρ) = ∑ₖ λₖ Eₖ ρ Eₖ†, 
    # where
    # - λₖ eigenvalues of matrix;
    # - Eₖ = ∑ₐ qₐ(k) τₐ, where q(k)=Q[:,k] is k-th column of Q and {τₐ} are the canonical operator basis. 
    krausoperators = Vector{Matrix{ComplexF64}}()
    for k = 1:length(eigvals)
        val = sum(Q[:, k] .* canonical_op_basis)
        push!(krausoperators, val)
    end
    return eigvals, krausoperators
end

"""
    dynamics_krausform(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}}, ρ::Matrix{ComplexF64}) -> Matrix{ComplexF64}

Applies a quantum operation to a density matrix ρ using its Kraus decomposition.

# Arguments
- `eigvals`: Eigenvalues `λₖ` from the decomposition.
- `E`: Kraus operators `Eₖ`.
- `ρ`: Density matrix to which the operation is applied.

# Returns
The transformed density matrix under the Kraus form:

    Φ(ρ) = ∑ₖ λₖ Eₖ ρ Eₖ†
"""
function dynamics_krausform(
    eigvals::Vector{Float64},
    E::Vector{Matrix{ComplexF64}},
    ρ::Matrix{ComplexF64},
)
    sum([eigvals[k] * E[k] * ρ * E[k]' for k in axes(eigvals, 1)])
end

"""
    computeKs(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})

Computes the effective Hamiltonian part Ks(t) of a time-local quantum generator at time t, 
given its Kraus decomposition.

# Arguments
- `eigvals`: eigenvalues λₖ(t) from the Kraus decomposition of the generator.
- `E`: matrix whose columns are the Kraus operators Eₖ(t).

# Returns
A `Matrix{ComplexF64}` representing the effective Hamiltonian `Ks(t)`, computed as:

    Ks(t) = - (i/4) * ∑ₖ λₖ(t) [Tr(Eₖ) * Eₖ† − Tr(Eₖ†) * Eₖ].
"""
function computeKs(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})
    total =
        sum([eigvals[k] * (tr(E[k]) * E[k]' - tr(E[k]') * E[k]) for k = 1:length(eigvals)])
    return -im / 4.0 * total
end

"""
    computeKs_check(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})

Alternative method to compute the effective Hamiltonian part Ks(t) of a time-local quantum generator 
using its action on a canonical operator basis.

# Arguments
- `eigvals`: Eigenvalues `λₖ` from the Kraus decomposition.
- `E`: Kraus operators `Eₖ`.

# Returns
A `Matrix{ComplexF64}` representing the effective Hamiltonian `Ks(t)`, computed via:

    Ks(t) = - (i/4) * ∑ₐ [ τₐ† * S(τₐ) − S(τₐ) * τₐ† ]

where `S(·)` is the superoperator defined by the Kraus decomposition, and `{τₐ}` is the canonical operator basis.
"""
function computeKs_check(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})
    total = sum([
        canonical_op_basis[a]' * dynamics_krausform(E, eigvals, canonical_op_basis[a]) -
        dynamics_krausform(E, eigvals, canonical_op_basis[a]) * canonical_op_basis[a]'
        for a in axes(canonical_op_basis, 1)
    ])
    return -im / 4.0 * total
end

"""
    computeKs(dirdata::String)

Compute the effective hamiltonian Ks(t) of a time-local quantum generator from measurement data in `dirdata`.  
Returns a named tuple with time-dependent matrices `Ks` and associated `time` vector.
"""
function computeKs(dirdata::String; check = false)
    # Map tomography: from data obtain the dynamical map Φ for every t in the Pauli basis.
    qmap = qmaptomography(dirdata)
    # Change dynamics description from dynamical map Φ to generator L = dΦ/dt Φ^(-1).
    L = qmap_to_gen(qmap.qmap, qmap.time)
    time = L.time
    gen = L.gen

    # Ks(t)
    Ks = Vector{Matrix{ComplexF64}}()
    for t = 1:length(time)
        # gen[t] is in the Pauli basis, rewrite in computational basis
        gen_canonical = changebasis(gen[t])
        # gen_canonical is an A-represention matrix of the quantum generator L,
        # kraus decomposition should be applied to B-represention matrix (Choi matrix).
        choimatrix = reshape(A2Bmatrix * vec(gen_canonical), (4, 4))
        # Kraus decomposition of the generator at time t
        eigvals, krausoperators = krausdecomposition(choimatrix)
        if check
            push!(Ks, computeKs_check(eigvals, krausoperators))
        else
            push!(Ks, computeKs(eigvals, krausoperators))
        end
    end
    return (Ks = Ks, time = time)
end

"""
    computeUs(dirdata::String)

Compute the real part of system-environment interaction energy for the "Up" initial state.

# Arguments
- `dirdata::String`: path to directory containing measurement and Hamiltonian data.

# Returns
- `Us`: vector of real interaction energies tr(Ks * ρs).
- `time`: corresponding time points.
"""
function computeUs(dirdata::String)
    effective_hamiltonian = computeKs(dirdata)
    # Load ρs(t) for system prepared in the "Up" state
    meas = get_measurements(dirdata * "/measurements_Up.dat", "densitymatrix")
    Up = meas.result
    # I may have a discontinuous behaviour in t=0 since Ks evaluate effective system while
    # Up[1] is the TLS prepared in a system that has not already interact with the env.
    Us = real(tr.(effective_hamiltonian.Ks .* Up[1:end-1]))
    return (Us = Us, time = ffective_hamiltonian.time)
end

"""
    chain_occupation(dirdata::String)

Extract occupation numbers of each site in the environment chain.

# Arguments
- `dirdata::String`: path to directory containing occupation measurements.

# Returns
- `ns`: matrix of occupation numbers; each row is referred to a different time.
- `sites`: site indices.
- `time`: array containing time of measurements.
"""
function chain_occupation(dirdata::String)
    meas = get_measurements(dirdata * "/measurements_N.dat", "N")
    # MPS of sys⊗env: environment sites from {2} to {NN+1}
    # with NN number of chain sites.
    # Measurments performed for each sites for each time in meas.time:
    # collect result into `data`.
    NN = length(meas.result)
    data = [[meas.result["N{$i}_re"][t] for i = 2:NN+1] for t = 1:length(meas.time)]
    sites = 1:NN
    return (ns = data, sites = sites, time = meas.time)
end

"""
Diagonalize extended chain Hamiltonian to obtain normal-mode frequencies and transformation matrix.

# Arguments
- `filename::String`: path to JSON config file.
- `extdendN`: extend chain length for fictitious modes, default 500.

# Returns
- `D`: Diagonal matrix of mode frequencies.
- `P`: Unitary transformation matrix (eigenvectors).
"""
function _normalmodes(filename::String; extendN = 500)
    # Chain Hamiltonian for the env is H_E = c† A c, where A is a three-diagonal matrix:
    # major diagonal formed by frequencies ωn (n = 0...N-1),
    # sub-diagonal and super-diagonal formed by couplings κn (n = 1...N-1).
    params = JSON.parsefile(filename)
    # I fictitiously enlarge the chain, see Riva et al.: PRB 108, 195138 (2023), Appendix E.
    # chain_lenght < PolyChaos_nquad !
    params["chain_length"] = extendN
    if params["PolyChaos_nquad"] < extendN
        params["PolyChaos_nquad"] = extendN * 10
    end
    coeff = chain_coefficients(params)
    f = coeff.frequencies
    # Couplings κn with n = 0...N-1. But κ0, the one for the TLS-env interaction,
    # must be excluded because it doesn't belong with H_E.
    # Thus, I consider κn with n = 1...N-1.
    k = coeff.couplings[2:end]
    A = Tridiagonal(k, f, k)

    # Any Hermitian (here also real) matrix can be diagonalized by a unitary matrix.
    #   A = U D Uᵀ,
    # with D being a diagonal matrix and U  being unitary matrix, Uᵀ=U⁻¹=U†.
    # Eigensolver:
    # F.values contains the eigenvalues (diagonal entries of D);
    # F.vectors contains the eigenvectors (columns of U).
    # Note: matrix A can be reconstructed as F.vectors * Diagonal(F.values) * inv(F.vectors).
    # Note: since U=F.vectors unitary, I have inv(U) = transpose(U) = U'
    F = eigen(A)
    D = Diagonal(F.values)
    U = F.vectors
    return D, U
end


function envmodes_occupation(meastime, U_squared, c_site_occupations)
    mode_occupations = Vector{Vector{Float64}}()
    NN = size(U_squared,1)
    T = length(meastime)
    for n in 1:NN  # loop over modes
        mode_values = Float64[]
        for t in 1:T  # loop over time
            occ = sum(U_squared[k, n] * c_site_occupations[k][t] for k in 1:NN)
            push!(mode_values, occ)
        end
        push!(mode_occupations, mode_values)
    end
    return mode_occupations
end

"""
    envmodes_occupation(dirdata::String)

Compute occupation of environment normal modes and write to CSV.
- Writes `envmodes_data.dat` with time and occupations.
- Writes `envmodes_modes.dat` with mode frequencies.

# Arguments
- `dirdata::String`: path to directory for measurements and output.
"""
function envmodes_occupation(dirdata::String; extendN = 700)
    # Any Hermitian (here also real) matrix can be diagonalized by a unitary matrix.
    #   A = U D Uᵀ,
    # with D being a diagonal matrix and U  being unitary matrix, Uᵀ=U⁻¹=U†.
    # (Optional)
    # Extend fictitiously the chain with 0s entries; chain length: N ↦ NN.
    D, U = _normalmodes(dirdata * "/config.json"; extendN = extendN)
    modes = diag(D)
    NN = length(modes)
    # I prefer to see A = U† D U. In fact, given Hₑ = c† A c 
    # with c is vector of annihilation operators,
    # normal modes decomposition of Hₑ is straightforwardly found as
    #   Hₑ = t† D t, with t:=Uc, t†:=c†U†
    # Occupation number of the new modes: 
    #   <tₙ† tₙ> = ∑k U[k,n]^2 * <cₖ† cₖ†>, with k = 0,...,N-1.
    U = transpose(U)
    # Precompute U^2 values (squared coefficients for each mode)
    U_squared = U.^2

    # MPS of sys⊗env: environment sites from {2} to {N+1}
    # <cⱼ† cⱼ> = meas["N{i}_re"] with i from 2 to N+1,
    # <cⱼ† cⱼ> = 0 with i from N+2 to NN+1.
    # Load measurement data: <c† c> values over time
    meas = get_measurements(dirdata * "/measurements_N.dat", "N")
    result, meastime = meas.result, meas.time
    # Build list of per-site occupation time series
    c_site_occupations = Vector{Vector{Float64}}()
    for site_index in 2:NN+1
        key = "N{$(site_index)}_re"
        values = haskey(result, key) ? result[key] : zeros(lenth(meastime))
        push!(c_site_occupations, values)
    end

    # Compute occupation numbers of each normal mode over time
    mode_occupations = envmodes_occupation(time, U_squared, c_site_occupations)

    # Write mode occupations over sampled time indices to file
    open("$dirdata/envmodes_data.dat", "w") do io
        # Header
        println(io, join(["time"; ["mode_$(i)" for i = 1:NN]], ","))
        # Transpose mode_occupations so each row is time → [mode_1, mode_2, ..., mode_N]
        for idx in eachindex(meastime)
            row = [meastime[idx]]
            for mode in mode_occupations
                push!(row, mode[idx])
            end
            println(io, join(row, ","))
        end
    end

    # Write normal modes ti file
    open("$dirdata/envmodes_modes.dat", "w") do io
        println(io, "frequency")
        for f in modes
            println(io, f)
        end
    end
end

function read_envmodes_occupation(dirdata::String)
    (rawdata, header) = readdlm("$dirdata/envmodes_modes.dat", ',', Float64, header = true)
    modes = rawdata[:] # flatten from Matrix to Vector
    (rawdata, header) = readdlm("$dirdata/envmodes_data.dat", ',', Float64, header = true)
    meastime, ns = rawdata[:, 1], rawdata[:, 2:end]
    return modes, ns, meastime
end

"""
    effective_freqs(dirdata::String, time::Vector{Float64})

Compute frequency transitions ωₛ′ = (E₊ - E₋)/ħ from the effective Hamiltonian Kₛ at specified time points.
"""
function effective_freqs(dirdata::String, time::Vector{Float64})
    # Compute Ks and align its values corresponding to the
    # specified time points, assuming all times in `time` 
    # exist in `Ks.time`.
    effective_hamiltonian = computeKs(dirdata)
    inds = [
        findfirst(x -> isapprox(x, t; atol = 1e-8), effective_hamiltonian.time) for
        t in time
    ]
    Ks = effective_hamiltonian.Ks[inds]
    # Print error message if sliced `ks` length doesn't match
    # `time` lenght.
    if length(Ks) !== length(time)
        println("Error: time points mismatch!")
        return 0
    end
    # Frequency transition of Ks (effective Hamiltonian):
    # ωs' = (E+ - E-)/ħ with E± eigvals of Ks.
    E = eigen.(Ks)
    Eminus = [eig.values[1] for eig in E]
    Eplus = [eig.values[2] for eig in E]
    freqs = Eplus - Eminus
    return freqs
end
