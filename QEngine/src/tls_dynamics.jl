"""
    densitymatrix(rx::Float64, ry::Float64, rz::Float64)

Given the component of the Bloch vector (real values), returns the denisty matrix.
"""
function densitymatrix(rx::Float64, ry::Float64, rz::Float64)
    return 0.5 * (σ0 + rx * σ1 + ry * σ2 + rz * σ3)
end

function qmaptomography(Up::Vector{Matrix{ComplexF64}}, Down::Vector{Matrix{ComplexF64}}, Plus::Vector{Matrix{ComplexF64}}, Trans::Vector{Matrix{ComplexF64}})
    # Compute evolved Pauli basis
    evolved_basis =
        [Up .+ Down, 2 .* Plus .- Up .- Down, 2 .* Trans .- Up .- Down, Up .- Down]

    # Compute the quantum map as a vector of matrices
    qmap = map(eachindex(Up)) do t
        [0.5 * tr(σi' * evolved_basis[j][t]) for σi in paulimatrices, j = 1:4]
    end
    qmap = [reshape(q, 4, 4) for q in qmap]
    return qmap
end

"""
    qmap_to_gen(qmap::Vector{Matrix{ComplexF64}}, dt::Float64)

Computes the time-local generator from a time-dependent quantum map, 
using finite difference methods. Integration step is required.
"""
function qmap_to_gen(qmap::Vector{Matrix{ComplexF64}}, dt::Float64)
    # The generator L(t) is related to the dynamical map Φ(t)
    # by means of L(t) = dΦ/dt Φ^(-1).
    N = size(qmap, 1)
    gen = Vector{Matrix{ComplexF64}}(undef, N)

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

    return gen
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
Returns eigenvalues and the corresponding Kraus operators.
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
    computeKs()

Compute the effective hamiltonian Ks(t) of a time-local quantum generator from measurement data in `dirdata`.  
Returns a named tuple with time-dependent matrices `Ks` and associated `time` vector.
"""
function computeKs(Up::Vector{Matrix{ComplexF64}}, Down::Vector{Matrix{ComplexF64}}, Plus::Vector{Matrix{ComplexF64}}, Trans::Vector{Matrix{ComplexF64}}, dt::Float64; check = false)
    # Map tomography: from data obtain the dynamical map Φ for every t in the Pauli basis.
    qmap = qmaptomography(Up, Down, Plus, Trans)
    # Change dynamics description from dynamical map Φ to generator L = dΦ/dt Φ^(-1).
    gen = qmap_to_gen(qmap, dt)

    # Ks(t)
    Ks = Vector{Matrix{ComplexF64}}()
    for t = 1:size(qmap, 1)
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
    effective_freqs()

Compute frequency transitions ωₛ′ = (E₊ - E₋)/ħ from the effective Hamiltonian Ks at specified time points.
"""
function effective_freqs(Ks::Vector{Matrix{Float64}})
    # Frequency transition of the effective Hamiltonian Ks is
    #   ωs' = (E₊ - E₋)/ħ 
    # with E± eigvals of Ks.
    E = eigen.(Ks)
    Eminus = [eig.values[1] for eig in E]
    Eplus = [eig.values[2] for eig in E]
    effcetive_freqs = Eplus - Eminus
    return effcetive_freqs
end

function computeUs(Ks::Vector{Matrix{Float64}}, state::Vector{Matrix{Float64}})
    # I may have a discontinuous behaviour in t=0 since Ks evaluate effective system while
    # Up[1] is the TLS prepared in a system that has not already interact with the env.
    Us = real(tr.(Ks .* state[1:end-1]))
    return
end


