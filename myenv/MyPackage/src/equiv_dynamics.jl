function threshold_complex(z::ComplexF64, threshold::Float64)
    re = abs(real(z)) < threshold ? 0.0 : real(z)
    im = abs(imag(z)) < threshold ? 0.0 : imag(z)
    return Complex(re, im)
end


function quantum_map_threshold_zero!(myMap; threshold=1e-5)
    for M in myMap
        for i in axes(M, 1), j in axes(M, 2)
            if i != j 
                M[i,j] = threshold_complex(M[i,j], threshold)
            end
        end
    end
    return myMap
end


function quantum_map(dir_path)
    # read meaurements data from file
    evolvedUp = get_evolved_states(dir_path * "/measurements_Up.dat")
    evolvedDown = get_evolved_states(dir_path * "/measurements_Dn.dat")
    evolvedPlus = get_evolved_states(dir_path * "/measurements_+.dat")
    evolvedTrans = get_evolved_states(dir_path * "/measurements_i.dat")
    # evolved pauli basis is retrieved from the known time evolution of the tomographic basis
    evolved_σ0 = evolvedUp + evolvedDown
    evolved_σ1 = 2 * evolvedPlus - evolvedUp - evolvedDown
    evolved_σ2 = 2 * evolvedTrans - evolvedUp - evolvedDown
    evolved_σ3 = evolvedUp - evolvedDown
    # define and compute the quantum map
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
    myMap = quantum_map_threshold_zero!(myMap)
    return myMap
end


function generator(myMap, dt)
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
