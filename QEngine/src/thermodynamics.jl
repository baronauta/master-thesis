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


# function internal_energy(Ks, rho)
#     return real(tr.(Ks .* rho[2:end])) # Ks(t=0) not defined
# end
