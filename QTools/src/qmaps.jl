include("constants.jl")

using LinearAlgebra

# ─────────────────────────────────────────────────────────────
# DENSITY MATRIX
#
# ρ = 1/2 ( σ0 + ∑ₖ rₖ ⋅ σₖ ), rₖ = Tr[ρ⋅σₖ]
# ─────────────────────────────────────────────────────────────

"Given the bloch vector, return the density matrix"
function densitymatrix(rx::Float64, ry::Float64, rz::Float64)
    return 0.5 * (paulimatrices[1] + rx * paulimatrices[2] + ry * paulimatrices[3] + rz * paulimatrices[4])
end


# ─────────────────────────────────────────────────────────────
# CHANGE BASIS
#
# From the pauli basis to the computational basis.
# ─────────────────────────────────────────────────────────────

"""
Change the representation of a matrix between orthonormal operator bases:
- FROM the Pauli basis {σ₀/√2, σ₁/√2, σ₂/√2, σ₃/√2};
- TO the canonical basis {|0⟩⟨0|, |0⟩⟨1|, |1⟩⟨0|, |1⟩⟨1|}.

Both bases are orthonormal with respect to the Hilbert-Schmidt inner product ⟨A, B⟩ = Tr[A†B].
"""
function changebasis(Φ::Matrix{ComplexF64})
    # Given Φ in Pauli basis, its representation in the canonical basis is
    # Λ† Φ Λ = 1/2 M† Φ M.
    # I define Λ = 1/√2 ⋅ M to avoid numerical roundings of √2.
    M = Matrix{ComplexF64}([1 0 0 1; 0 1 1 0; 0 im -im 0; 1 0 0 -1])
    Φnew = 1 / 2 * adjoint(M) * Φ * M

    # Sanity check
    # Φnew = Λ† Φ Λ → Λ Φnew Λ† = Λ Λ† Φ Λ Λ† = Φ → 1/2 M Φnew M† = Φ
    @assert 0.5 * M * Φnew * adjoint(M) ≈ Φ

    return Φnew
end


# ─────────────────────────────────────────────────────────────
# QUANTUM MAP TOMOGRAPHY
#
# The Pauli matrices (basis for 2×2 Hermitian operators) can be
# be written in term of the states of the the tomographic states,
# i.e. a state of state for whom the time evolution is known. 
# ─────────────────────────────────────────────────────────────

function qmaptomography(
    Up::Vector{Matrix{ComplexF64}},
    Down::Vector{Matrix{ComplexF64}},
    Plus::Vector{Matrix{ComplexF64}},
    Trans::Vector{Matrix{ComplexF64}},
)
    # Time-dependent Pauli matrices reconstructed from tomographic states
    σ₀_t = Up .+ Down
    σ₁_t = 2 .* Plus .- Up .- Down
    σ₂_t = 2 .* Trans .- Up .- Down
    σ₃_t = Up .- Down

    evolved_pauli = [σ₀_t, σ₁_t, σ₂_t, σ₃_t]

    qmap = [zeros(ComplexF64, 4, 4) for _ in eachindex(Up)]
    for t in axes(qmap, 1)
        # Φᵢⱼ(t) = 1/2 Tr[σᵢ† σⱼ(t)]
        for i = 1:4, j = 1:4
            qmap[t][i,j] = 0.5 * tr( adjoint(paulimatrices[i]) * evolved_pauli[j][t] )
        end
    end
    return qmap
end


# ─────────────────────────────────────────────────────────────
# GENERATOR (from quantum dynamical map)
#
# The generator L(t) is related to the dynamical map Φ(t) by 
# means of L(t) = (∂ₜΦₜ) Φₜ⁻¹.
# ─────────────────────────────────────────────────────────────

function qmap2gen(qmap::Vector{Matrix{ComplexF64}}, dt::Float64)

    # Compute Φ⁻¹
    inv_qmap = inv.(qmap)
    # Sanity check: Φ⁻¹ ⋅ Φ = 1
    I4 = Matrix{ComplexF64}(I, 4, 4)
    for t in axes(qmap, 1)
        @assert inv_qmap[t] * qmap[t] ≈ I4 "Failed at t = $(t * dt)"
    end

    # Compute ∂ₜ Φₜ
    d_qmap = [zeros(ComplexF64, 4, 4) for _ in axes(qmap, 1)]
    for t in axes(qmap, 1)
        if t == 1
            # First point: forward difference
            d_qmap[t] = (qmap[t+1] - qmap[t]) / dt
        elseif t == size(qmap, 1)
            # Last point: backward difference
            d_qmap[t] = (qmap[t] - qmap[t-1]) / dt
        else
            # Interior points: centered difference
            d_qmap[t] = (qmap[t+1] - qmap[t-1]) / (2 * dt)
        end
    end

    # Compute L(t) = (∂ₜΦₜ) Φₜ⁻¹
    gen = d_qmap .* inv_qmap
    # Sanity check: L(t) ⋅ Φₜ = ∂ₜΦₜ
    for t in axes(qmap, 1)
        @assert gen[t] * qmap[t] ≈ d_qmap[t] "Failed at t = $(t * dt)"
    end

    return gen
end


# ─────────────────────────────────────────────────────────────
# KRAUS DECOMPOSITION
#
# M is Hermitian ⇒ admits spectral decomposition:
# i.e. there exists a unitary Q and real eigenvalues λₖ such that
#       M = Q Λ Q† = ∑ₖ λₖ q(k) q(k)†,
# where the k-th column of Q is the eigenvector q(k) of M,
# and λₖ are the corresponding eigevalues.
#
# Kraus decomposition of a superoperator S: S(H) → S(H):
#       S(ρ) = ∑ₖ λₖ Eₖ ρ Eₖ†, 
# where
# - λₖ eigenvalues of matrix;
# - Eₖ = ∑ₐ qₐ(k) τₐ, {τₐ} canonical operator basis. 
# ─────────────────────────────────────────────────────────────

function krausdecomposition(M::Matrix{ComplexF64})

    @assert M ≈ adjoint(M) "Input matrix is not Hermitian"
    
    E = eigen(M)
    eigvals = E.values    # eigenvalues
    Q = E.vectors         # eigenvectors (stored as columns): Q[:, k]

    @assert all(isreal, eigvals) "Eigenvalues are not real"
    @assert M ≈ sum([ eigvals[k] * Q[:,k] * adjoint(Q[:, k]) for k in axes(eigvals, 1) ]) 

    krausoperators = [zeros(ComplexF64, 4, 4) for _ in axes(eigvals, 1)]
    for k = axes(eigvals, 1)
        val = sum([ Q[a, k] * canonical_op_basis[a] for a in axes(Q, 2) ])
        krausoperators[k] = val
    end

    return eigvals, krausoperators
end


function krausform(
    eigvals::Vector{Float64},
    E::Vector{Matrix{ComplexF64}},
    ρ::AbstractMatrix,
)
    # S(ρ) = ∑ₖ λₖ Eₖ ρ Eₖ†
    return sum([eigvals[k] * E[k] * ρ * adjoint(E[k]) for k in axes(eigvals, 1)])
end


# ─────────────────────────────────────────────────────────────
# EFFECTIVE HAMILTONIAN Ks
#
# Effective Hamiltonian Ks in computed from the pseudo-Kraus
# decomposition of the time-local generator of the dynamics:
# Ks(t) = - (i/4) * ∑ₖ λₖ(t) [Tr(Eₖ) * Eₖ† − Tr(Eₖ†) * Eₖ].
#
# Other way of calculating Ks:
# Ks(t) = - (i/4) * ∑ₐ [ τₐ† * S(τₐ) − S(τₐ) * τₐ† ],
# where S(ρ) = ∑ₖ λₖ Eₖ ρ Eₖ†, {τₐ} canonical operator basis. 
# ─────────────────────────────────────────────────────────────

function effective_hamiltonian(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})
    total =
        sum([
            eigvals[k] * (
                tr(E[k]) * adjoint(E[k]) - 
                tr(adjoint(E[k])) * E[k]
            )
            for k = axes(eigvals, 1)
        ])
    return -im / 4.0 * total
end


function effective_hamiltonian_check(eigvals::Vector{Float64}, E::Vector{Matrix{ComplexF64}})
    total = 
    sum([
        adjoint(canonical_op_basis[a]) * krausform(eigvals, E, canonical_op_basis[a]) -
        krausform(eigvals, E, canonical_op_basis[a]) * adjoint(canonical_op_basis[a])
        for a = axes(canonical_op_basis, 1)
    ])
    return -im / 4.0 * total
end


function effective_hamiltonian_pipeline(
    Up::Vector{Matrix{ComplexF64}},
    Dn::Vector{Matrix{ComplexF64}},
    Plus::Vector{Matrix{ComplexF64}},
    Trans::Vector{Matrix{ComplexF64}},
    dt::Float64;
)
   # Map tomography: from data obtain the dynamical map Φₜ expressed in the Pauli basis.
    qmap = qmaptomography(Up, Dn, Plus, Trans)

    # Kraus decomposition of the map
    # ...

    # Change dynamics description from dynamical map Φ to generator L(t) = (∂ₜΦₜ) Φₜ⁻¹
    gen = qmap2gen(qmap, dt)

    Ks = [ zeros(ComplexF64, 2, 2) for _ in axes(gen, 1)]
    Ks_check = [ zeros(ComplexF64, 2, 2) for _ in axes(gen, 1)]
    # Psuedo-Kraus decomposition of the generator
    for t in axes(gen, 1)
        # Change to computational basis
        gen_canonical = changebasis(gen[t])
        # Reshuffle the matrix from A-form to B-form (Choi matrix)
        choimatrix = reshape(A2Bmatrix * vec(gen_canonical), (4, 4))
        # Pseudo-Kraus decomposition of the generator at time t
        eigvals, pseudokraus = krausdecomposition(choimatrix)

        # Sanity check: ∑ₖ λₖ(t) Eₖ†(t) Eₖ(t) = 0 ∀t
        @assert isapprox(
            sum([
                eigvals[k] * adjoint(pseudokraus[k]) * pseudokraus[k]
                for k in axes(eigvals, 1)
            ]),
            zeros(ComplexF64, 2, 2);
            atol = 1e-6
        ) "Failed at t = $(t * dt)"

        Ks[t] = effective_hamiltonian(eigvals, pseudokraus)
        Ks_check[t] = effective_hamiltonian_check(eigvals, pseudokraus)

        @assert Ks[t] ≈ Ks_check[t]
    end

    return Ks
end


function write_Ks(dir::AbstractString, Ks::Vector{Matrix{ComplexF64}})
    open(joinpath(dir, "Ks_matrices.dat"), "w") do io
        for K in Ks
            flat_K = vec(K)                 # Flatten column-wise
            writedlm(io, [flat_K'], ',')    # Transpose to write as one row
        end
    end
end


"""
Given a 2×2 Hamiltonian compute the transition frequency,
i.e. ω = (E₊ - E₋)/ħ where E± are its eigenvalues
"""
function transition_freqs(H::Matrix{ComplexF64})
    E = eigen(H)
    Eminus = E.values[1]
    Eplus = E.values[2]
    # ħ is 1
    return Eplus - Eminus
end

function computeUs(Ks::Vector{Matrix{Float64}}, state::Vector{Matrix{Float64}})
    # I may have a discontinuous behaviour in t=0 since Ks evaluate effective system while
    # Up[1] is the TLS prepared in a system that has not already interact with the env.
    Us = real(LinearAlgebra.tr.(Ks .* state[1:end-1]))
end



