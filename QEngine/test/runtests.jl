using QEngine
using LinearAlgebra
using Test

@testset "Depolarizing channel" begin

    # Depolarizing channel: Φ(ρ) = (1-p) * ρ + p * I/2
    function depolarizing_channel(ρ, p)
        (1 - p) * ρ + p * I / 2
    end

    # Depolarizing channel expressed in terms of Pauli matrices:
    # Φ(ρ) = (1 - 3p/4) * ρ + p/4 * (XρX + YρY + ZρZ).
    # Bloch vector is trasformed as
    # rₖ = Tr[ρσₖ] → rₖ' = Tr[Φ(ρ)σₖ] = (1-p)rₖ.
    function depolarizing_channel_pauli(p)
        D = Diagonal([1.0, 1 - p, 1 - p, 1 - p])
        return Matrix{ComplexF64}(D)
    end

    p = 0.1
    Φ_pauli = depolarizing_channel_pauli(p)
    Φ_canonical = QEngine.changebasis(Φ_pauli)

    # Compare action on a test state ρ
    ρ₊ = Matrix{ComplexF64}([0.5 0.5; 0.5 0.5])
    vecρ₊ = vec(ρ₊)
    vecρ_out = Φ_canonical * vecρ₊
    ρ_out = reshape(vecρ_out, 2, 2)
    ρ_expected = depolarizing_channel(ρ₊, p)
    @test norm(ρ_out - ρ_expected) ≈ 0

    # Kraus decomposition
    # Φ_canonical is an A-represention matrix of the quantum operation Φ,
    # kraus decomposition should be applied to B-represention matrix (Choi matrix).
    choimatrix = reshape(QEngine.A2Bmatrix * vec(Φ_canonical), (4, 4))
    eigvals, krausoperators = QEngine.krausdecomposition(choimatrix)
    # Check: λₖ ∈ Real 
    @test all(isreal, eigvals)
    # Check: ∑ₖ λₖ Eₖ† Eₖ = I
    @test sum([
        eigvals[k] * krausoperators[k]' * krausoperators[k] for k in eachindex(eigvals)
    ]) ≈ I
    # Kraus operator for depolarizing channel are known:
    # Φ(ρ) = ∑ₖ Mₖ ρ Mₖ†, with M₀=√(1-3p/4)σ₀, Mⱼ=√(p/4)σⱼ.
    # Check: √(λₖ) Eₖ = Mₖ...problems in finding appropriate matching.
    # @test eigvals[1] * krausoperators[1] ≈ sqrt(p/4) * paulimatrices[2]
    # @test eigvals[2] * krausoperators[2] ≈ sqrt(p/4) * paulimatrices[3]
    # @test eigvals[3] * krausoperators[3] ≈ sqrt(p/4) * paulimatrices[4]
    # @test sqrt(eigvals[4]) * krausoperators[4] ≈ sqrt(1-3*p/4) * paulimatrices[1]
end

@testset "Normal modes occupation" begin
    # Consider a chain of 4 sites with the given frequencies and couplings
    freqs = [
        0.7499999986477128,
        0.5833333202609593,
        0.5416666134711797,
        0.524999851201351,
    ]
    coups = [
        0.2303294331885381,
        0.1936491697544971,
        0.22537448194777276,
        0.23622784163992785,
    ]
    # Measurements of N on each sites for t = 1, 2, 3,
    # at t=1 each site is in the vaccum state.
    # measN is a matrix of dimension time × sites, i.e. 3 × 4.
    measN = [
        0.0 0.0 0.0 0.0;
        10.0 0.0 0.0 0.0;
        20.0 3.0 0.1 0.0
    ]

    # A = Tridiagonal(coups, freqs, coups) where coups = coups[2:end],
    # i.e. exclude the coupling sys-env.
    k = coups[2:end]
    A = Tridiagonal(k, freqs, k)
    # Diagonalize A, i.e. A = P D P†.
    # D is a diagonal matrix whose entries are the eigenvals of A.
    # P is a unitary matrix whose columns are the eigenvectors of A.
    D, P = QEngine.normalmodes(freqs, coups)
    # Check A = P D P†
    @test A ≈ P * D * P'
    # A eigvec = eigval eigvec
    @test all([A * P[:,i] ≈ D[i,i] * P[:,i] for i in 1:length(freqs)])
    # P is unitary, so Pᵀ=P⁻¹=P†
    @test transpose(P) ≈ P'
    @test inv(P) ≈ P'
    # I prefer to see A = U† D U where U = P†.
    U = transpose(P)
    @test A ≈ U' * D * U

    modes, occupation = QEngine.envmodes_occupation(freqs, coups, measN)
    @test modes ≈ diag(D)
    # t=1 (first row) the sites are in the vacuum state,
    # then also normal modes are in the vacuum state
    @test occupation[1,:] ≈ [0.0; 0.0; 0.0; 0.0]
    # t=2: <tₙ†tₙ> = ∑ₖ |Uₙₖ|² <cₖ†cₖ> = |Uₙ₁|² <c₁†c₁> because only <c₁†c₁>≠0
    Un1 = U[:,1]
    @test occupation[2,:] ≈ (Un1 .^2) * measN[2,1]
    # t=3: <tₙ†tₙ> = ∑ₖ |Uₙₖ|² <cₖ†cₖ> for k=1,2,3
    @test occupation[3,:] ≈ (U[:,1] .^2) * measN[3,1] + (U[:,2] .^2) * measN[3,2] + (U[:,3] .^2) * measN[3,3]
end
