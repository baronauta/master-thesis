function effective_hamiltonian(choi_eigvals::Vector{Vector{ComplexF64}}, EE::Vector{Vector{Matrix{ComplexF64}}})
    Ks = Vector{Matrix{ComplexF64}}()
    # push iniziale con Ks al tempo zero?
    for t in axes(EE, 1)
        push!(
            Ks,
            -im / 4.0 * sum([
                choi_eigvals[t][k] * (tr(EE[t][k]) * EE[t][k]' - tr(EE[t][k]') * EE[t][k])
                for k = 1:size(choi_eigvals[t], 1)
            ]),
        )
    end
    return Ks
end


function effective_hamiltonian(map_tomo_path)
    # read config file
    config_path = map_tomo_path * "/config.JSON"
    config = load_config(config_path)
    # compute the map and its generator
    myMap = quantum_map(map_tomo_path)
    myGen = generator(myMap, config.dt)
    # choi matrix of the generator and then kraus decomposition
    choi_eigvals, EE = kraus_decomposition(myGen)
    # effective hamiltonian
    return effective_hamiltonian(choi_eigvals, EE)
end


function effective_hamiltonian_check(map_tomo_path)
    # read config file
    config_path = map_tomo_path * "/config.JSON"
    config = load_config(config_path)
    # compute the map and its generator
    myMap = quantum_map(map_tomo_path)
    myGen = generator(myMap, config.dt)
    # choi matrix of the generator and then kraus decomposition
    choi_eigvals, EE = kraus_decomposition(myGen)
    # effective hamiltonian
    Ks_check = Vector{Matrix{ComplexF64}}()
    # push iniziale con Ks al tempo zero?
    for t in axes(EE, 1)
        push!(
            Ks_check,
            -im / 4.0 * sum([
                T_basis[a]' * kraus_term(EE[t], choi_eigvals[t], T_basis[a]) -
                kraus_term(EE[t], choi_eigvals[t], T_basis[a]) * T_basis[a]' for
                a in axes(T_basis, 1)
            ]),
        )
    end
    return Ks_check
end


function internal_energy(Ks, rho)
    return real(tr.(Ks .* rho[2:end])) # Ks(t=0) not defined
end
