# ─────────────────────────────────────────────────────────────
# Quantum dynamics
#
# Let B(H) be the space of bounded linear operators on Hilbert space H.
# Dynamical map Φ(t): B(H) → B(H);
# dynamics of a state is described by ρ(t) = Φ(t)[ρ(0)].
# Generator L(t): B(H) → B(H);
# dynamics of a state is described by d/dt ρ(t) = L(t)[ρ(t)].
#
# ─────────────────────────────────────────────────────────────

function qmaptomography(dirdata)
    # Load measurements and extract results
    meas_files = ["Up", "Dn", "+", "i"]
    meas_data = Dict(
        name => get_measurements(dirdata * "/measurements_$name.dat", "densitymatrix")
        for name in meas_files
    )
    Up, Down, Plus, Trans = meas_data["Up"].result,
    meas_data["Dn"].result,
    meas_data["+"].result,
    meas_data["i"].result
    # Extract time (assumes all share the same timeline)
    time = meas_data["Up"].time

    # Compute evolved Pauli basis
    evolved_basis =
        [Up .+ Down, 2 .* Plus .- Up .- Down, 2 .* Trans .- Up .- Down, Up .- Down]

    # Compute the quantum map as a vector of matrices
    qmap = map(eachindex(Up)) do t
        [0.5 * tr(σi' * evolved_basis[j][t]) for σi in pauli_basis, j = 1:4]
    end
    qmap = [reshape(q, 4, 4) for q in qmap]

    # Return results as named tuple
    return (qmap = qmap, time = time)
end

function qmap_to_gen(qmap::Vector{Matrix{ComplexF64}}, time::Vector{Float64})
    # The generator L(t) is related to the dynamical map Φ(t)
    # by means of L = dΦ/dt Φ^(-1).
    gen = Vector{Matrix{ComplexF64}}()

    # Integration step (assumes constant timestep)
    dt = time[2] - time[1]
    # Derivative computed as forward difference, i.e. df(t)/dt = {f(t+1)-f(t)} / dt.
    # (df(t)/dt)_{t=t_{n+1}} = {f(t_{n+1})-f(t_n)}/dt is not defined 
    # because t_{n+1} is not defined.
    for t = 1:length(time)-1
        derivative = (qmap[t+1] - qmap[t]) / dt
        inverse = inv(qmap[t])
        push!(gen, derivative * inverse)
    end

    return (gen = gen, time = time[1:end-1])
end

function _changebasis!(
    Φ::Matrix{ComplexF64},
    from_basis::Vector{Matrix{ComplexF64}},
    to_basis::Vector{Matrix{ComplexF64}},
)
    # Hilbert-Schmidt inner product: <A,B> = Tr[A†B]
    Λ = [tr(A' * B) for A in from_basis, B in to_basis]
    return Λ' * Φ * Λ
end

function krausdecomposition(gen::Matrix{ComplexF64})
    # Generator computed with `qmap_to_gen` is in the pauli basis, I want to work
    # in the computational basis (here named `T_basis`).
    gen = _changebasis!(gen, pauli_basis, T_basis)
    # The matrix representation of the generator is in the A-form, that is not Hermitian.
    # Reshuffle its elements to obtain the B-form matrix representation: 
    # the matrix representing the generator is now Hermitian, thus it admits a spectral decomposition.
    # The B-form matrix representation of the generator coincides with the Choi matrix.
    choimatrix = reshape(A2Bmatrix * vec(gen), (4, 4))

    # Eigensolver
    E = eigen(choimatrix)
    # E.values contains the eigenvalues
    eigvals = E.values
    # E.vectors contains the eigenvectors (stored as columns)
    eigcols = E.vectors

    # Kraus decomposition of the generator: L(t)[̢ρ(t)] = ∑k λk Ek(t) ρ(t) Ek†(t).
    # - λk eigenvalues of the choi matrix
    # - Ek Kraus operator is obtained as Ek = ∑a v(k) Ta where v(k) is a eigenvectors of the Choi matrix.
    krausoperators = Vector{Matrix{ComplexF64}}()
    for k = 1:length(eigvals)
        val = sum(eigcols[:, k] .* T_basis)
        push!(krausoperators, val)
    end
    return eigvals, krausoperators
end

# ─────────────────────────────────────────────────────────────
# Thermodynamics quantities
# - Effective Hamiltonian Ks
# - Internal energy Us(t) = Tr[Ks(t)ρs(t)]
# ─────────────────────────────────────────────────────────────

function computeKs(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})
    total =
        sum([eigvals[k] * (tr(E[k]) * E[k]' - tr(E[k]') * E[k]) for k = 1:length(eigvals)])
    return -im / 4.0 * total
end

function computeKs(dirdata::String)
    # Execute map tomography given the density matrix of the states composing a tomographic basis.
    qmap = qmaptomography(dirdata)
    # Change dynamics description from dynamical map Φ to generator L = dΦ/dt Φ^(-1).
    gen = qmap_to_gen(qmap.qmap, qmap.time)

    # Ks(t)
    time = gen.time
    L = gen.gen
    Ks = Vector{Matrix{ComplexF64}}()
    for t = 1:length(time)
        # Kraus decomposition of the generator at time t
        eigvals, eigvecs = krausdecomposition(L[t])
        push!(Ks, computeKs(eigvals, eigvecs))
    end
    return (Ks = Ks, time = time)
end

# function effective_hamiltonian_check(map_tomo_path)

#     config = load_system_config(dirdata * "/config.JSON")
#     # compute the map and its qmap_to_gen
#     qmap = qmaptomography(map_tomo_path)
#     myGen = qmap_to_gen(qmap, config.dt)
#     # choi matrix of the qmap_to_gen and then kraus decomposition
#     choi_eigvals, EE = kraus_decomposition(myGen)
#     # effective hamiltonian
#     Ks_check = Vector{Matrix{ComplexF64}}()
#     # push iniziale con Ks al tempo zero?
#     for t in axes(EE, 1)
#         push!(
#             Ks_check,
#             -im / 4.0 * sum([
#                 T_basis[a]' * kraus_term(EE[t], choi_eigvals[t], T_basis[a]) -
#                 kraus_term(EE[t], choi_eigvals[t], T_basis[a]) * T_basis[a]' for
#                 a in axes(T_basis, 1)
#             ]),
#         )
#     end
#     return Ks_check
# end


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


# ─────────────────────────────────────────────────────────────
# Environment
# - Occupation number of the chain modes
# - Occupation number of the normal modes
# ─────────────────────────────────────────────────────────────

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

function _normalmodes(filename::String; extdendN = false)
    # Chain Hamiltonian for the env is H_E = c† A c, where A is a tridiagonal matrix:
    # major diagonal formed by frequencies ωn (n = 0...N-1),
    # sub-diagonal and super-diagonal formed by couplings κn (n = 1...N-1).
    params = JSON.parsefile(filename)
    # I fictitiously enlarge the chain, see Riva et al.: PRB 108, 195138 (2023), Appendix E.
    if extdendN == true
        params["PolyChaos_nquad"] = 10000
        params["chain_length"] = 5000
    end
    coeff = chain_coefficients(params)
    f = coeff.frequencies
    # Couplings κn with n = 0...N-1. But κ0, the one for the TLS-env interaction,
    # must be excluded because it doesn't belong with H_E.
    # Thus, I consider κn with n = 1...N-1.
    k = coeff.couplings[2:end]
    A = Tridiagonal(k, f, k)

    # Any Hermitian matrix can be diagonalized by a unitary matrix.
    # A = U† D U, with D diagonal matrix and U unitary matrix. Let P = U†.
    # Eigensolver:
    # E.values contains the eigenvalues (diagonal entries of D);
    # E.vectors contains the eigenvectors (columns of P).
    E = eigen(A)
    D = Diagonal(E.values)
    P = E.vectors
    return D, P
end

function envmodes_occupation(dirdata::String)
    # Measurements performed for each sites for each time in meas.time.
    # MPS of sys⊗env: environment sites from {2} to {N+1}
    # with N number of chain sites.
    # Let b{i} be the environment modes, occupation number is known, 
    # i.e <b†{i}b{i}> = meas["N{i}_re"] with i from 2 to N+1.
    meas = get_measurements(dirdata * "/measurements_N.dat", "N")
    # Compute chain coefficients and find normal modes.
    # Set extendN to `true` for computing normal modes for a longer chain
    # with NN modes.
    D, P = _normalmodes(dirdata * "/config.json"; extdendN = true)
    modes = diag(D)
    NN = length(modes)
    # The chain Hamiltonian of the environment is quadratic in b{i} 
    # with matrix A being tridiagonal, i.e H_E = (b†{i})_i A (b{i})_i.
    # Diagonalize A, that is  A = PDP^(-1) with D diagonal;
    # equivalent notation A = U†DU with U†:=P.
    # Normal modes decomposition of H_E is obtained defining
    # (t{i})_i:= U(b{i})_i. Thus, H_E = (t†{i})_i D (t{i})_i. 
    # Occupation number: 
    # <t†{i}t{i}> = ∑j P[i-1,j-1]^2 * <b†{j}b{j}>.
    # Modes numbering is i = 2...NN+1, P is a square matrix with indexes 1...NN.
    # Only <b†{i}b{i}>, with i from 1 to N+1, are measured, fill the others,
    # i.e. i from N+2 to NN+1, with 0s entries.
    data = [
        [
            sum(
                P[j-1, i-1]^2 * get(meas.result, "N{$j}_re", zeros(length(meas.time)))[t] for j = 2:NN+1
            ) for i = 2:NN+1
        ] for t = 1:length(meas.time)
    ]

    open("$dirdata/envmodes_data.dat", "w") do io
        # Write header
        println(io, join(["time"; ["mode_$(i)" for i = 1:length(modes)]], ","))
        # Write each row
        for (i, row) in enumerate(data)
            println(io, join([meas.time[i]; row], ","))
        end
    end

    open("$dirdata/envmodes_modes.dat", "w") do io
        # Write header
        println(io, "frequency")
        # Write each frequency
        for f in modes
            println(io, f)
        end
    end

end
