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

"Read measurements from file and return the evolved state for each time."
function get_evolved_states(filename)
    (rawdata, header) = readdlm(filename, ',', Float64, header = true)
    meas = Dict(header[i] => rawdata[:, i] for i in eachindex(header))
    nmeas = length(meas["time"])
    rho_t = [
        0.5 * (
            σ0 +
            2 * meas["Sx{1}_re"][i] * σ1 +
            2 * meas["iSy{1}_im"][i] * σ2 +
            2 * meas["Sz{1}_re"][i] * σ3
        ) for i = 1:nmeas
    ]
    return rho_t
end

"Read the normalization of the time evolved state: if everything works good it should be 1."
function get_norm(filename)
    (rawdata, header) = readdlm(filename, ',', Float64, header = true)
    meas = Dict(header[i] => rawdata[:, i] for i in eachindex(header))
    return meas["Norm_re"]
end

"Compute the quantum map in the Pauli basis."
function quantum_map(dir_path)
    # Read meaurements data from file
    evolvedUp = get_evolved_states(dir_path * "/measurements_Up.dat")
    evolvedDown = get_evolved_states(dir_path * "/measurements_Dn.dat")
    evolvedPlus = get_evolved_states(dir_path * "/measurements_+.dat")
    evolvedTrans = get_evolved_states(dir_path * "/measurements_i.dat")
    # Evolved Pauli basis is retrieved from the known time evolution of the tomographic basis
    evolved_σ0 = evolvedUp + evolvedDown
    evolved_σ1 = 2 * evolvedPlus - evolvedUp - evolvedDown
    evolved_σ2 = 2 * evolvedTrans - evolvedUp - evolvedDown
    evolved_σ3 = evolvedUp - evolvedDown
    # Define and compute the quantum map
    myMap = Vector{Matrix{ComplexF64}}(undef, length(evolvedUp))
    for t in 1:length(evolvedUp)
        evolved_basis = [evolved_σ0[t], evolved_σ1[t], evolved_σ2[t], evolved_σ3[t]]
        myMap[t] = zeros(ComplexF64, 4, 4)
        for (i, σi) in enumerate(pauli_basis)
            for (j, evolved_σj) in enumerate(evolved_basis)
                myMap[t][i, j] = 0.5 * tr(σi' * evolved_σj)  # A' is the hermitian transpose of A
            end
        end
    end
    return myMap
end

"Compute the generator of the dynamics."
function generator(myMap::Vector{Matrix{ComplexF64}}, dt)
    # derivative of the map
    derivative = Vector{Matrix{ComplexF64}}()
    for t = 1:size(myMap, 1)-1
        push!(derivative, (myMap[t+1] - myMap[t]) / dt) # forward difference
    end
    # inverse of the map
    inverse = Vector{Matrix{ComplexF64}}()
    for t = 1:size(derivative, 1)
        push!(inverse, inv(myMap[t]))
    end
    # generator
    return derivative .* inverse
end

"Change basis and compute the choi matrix of the dynamics."
function choi_matrix(myGen)
    # change of basis (Λ)ij = << σi | T >>
    # ATT! Inner product tr(A'*B); pauli/sqrt(2) is orthonormal, T is orthonormal
    Λ = zeros(Complex{Float64}, 4, 4)
    for (i, σi) in enumerate(pauli_basis)
        for (j, Tj) in enumerate(T_basis)
            Λ[i, j] = 1 / sqrt(2) * tr(σi' * Tj)
        end
    end
    # I have to reshuffle the representation of the generator in order to obtain the Choi matrix
    P = Matrix{ComplexF64}(
        [
            #   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 # 1
            0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 # 2
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 # 3
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 # 4
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 # 5
            0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 # 6
            0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 # 7
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 # 8
            0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 # 9
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 # 10
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 # 11
            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 # 12
            0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 # 13
            0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 # 14
            0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 # 15
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1  # 16
        ],
    )
    # choi matrix
    choi = Vector{Matrix{ComplexF64}}()
    for t in axes(myGen, 1)
        basis_change = Λ' * myGen[t] * Λ
        change_representation = reshape(P * vec(basis_change), (4, 4))
        push!(choi, change_representation)
    end
    return choi
end


function kraus_decomposition(myGen)
    # choi matrix from the generator
    choi = choi_matrix(myGen)
    # eigenvectors decomposition of the choi matrix
    choi_eigvals = Vector{Vector{ComplexF64}}()
    choi_eigvecs = Vector{Matrix{ComplexF64}}()
    for t in axes(choi, 1)
        eig_solver = eigen(choi[t])
        push!(choi_eigvals, eig_solver.values)
        push!(choi_eigvecs, eig_solver.vectors)
    end
    # kraus decomposition
    EE = Vector{Vector{Matrix{ComplexF64}}}()
    for t = 1:size(choi, 1)
        E = Vector{Matrix{ComplexF64}}()
        for k = 1:size(choi_eigvecs[t], 2)
            push!(E, sum(choi_eigvecs[t][:, k] .* T_basis))
        end
        push!(EE, E)
    end
    return choi_eigvals, EE
end

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



function matrix2vector(mat::Matrix{Complex{Float64}}, basis)
    vec = Complex{Float64}[]
    for elem in basis
        push!(vec, 0.5 * tr(mat' * elem)) # ATT! 0.5 * tr(A' * B) inner product for pauli basis to be orthonormal
    end
    vec
end


function vector2matrix(vec::Vector{Complex{Float64}}, basis)
    sum(vec .* basis)
end


function map_dynamics(map_tomo_path, init_state)
    myMap = quantum_map(map_tomo_path)
    vec_evolved_map = Vector{Vector{Complex{Float64}}}()
    for i = 1:size(myMap)[1]
        push!(vec_evolved_map, myMap[i] * matrix2vector(init_state, pauli_basis))
    end
    return map(x -> vector2matrix(x, pauli_basis), vec_evolved_map)
end


function me_dynamics(map_tomo_path, init_state, dt)
    myGen = generator(quantum_map(map_tomo_path), dt)
    vec_evolved_me = Vector{Vector{Complex{Float64}}}()
    push!(vec_evolved_me, matrix2vector(init_state, pauli_basis))
    for i = 1:size(myGen)[1]
        push!(vec_evolved_me, vec_evolved_me[i] + dt * myGen[i] * vec_evolved_me[i]) # eulero integration
    end
    return map(x -> vector2matrix(x, pauli_basis), vec_evolved_me)
end


function kraus_term(kraus_operator, eigvals, state)
    sum([
        eigvals[k] * kraus_operator[k] * state * kraus_operator[k]' for
        k in axes(eigvals, 1)
    ])
end


function kd_dynamics(map_tomo_path, init_state, dt)
    choi_eigvals, EE = kraus_decomposition(generator(quantum_map(map_tomo_path), dt))
    evolved_kd = Vector{Matrix{ComplexF64}}()
    push!(evolved_kd, init_state)
    for t in axes(choi_eigvals, 1)
        push!(
            evolved_kd,
            evolved_kd[t] + dt * kraus_term(EE[t], choi_eigvals[t], evolved_kd[t]),
        )
    end
    return evolved_kd
end
