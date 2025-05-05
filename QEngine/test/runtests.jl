using QEngine
using LinearAlgebra
using Test

import QEngine: changebasis, A2Bmatrix, krausdecomposition

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
    Φ_canonical = changebasis(Φ_pauli)

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
    choimatrix = reshape(A2Bmatrix * vec(Φ_canonical), (4, 4))
    eigvals, krausoperators = krausdecomposition(choimatrix)
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
