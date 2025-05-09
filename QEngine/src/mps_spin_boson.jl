"""
    defineSystem(sys_state::String, chain_length::Int64, local_dim::Int64)

Define the MPS of a system composed by a two-level system prepared in the `sys_state` (e.g. "Up")
and a chain of harmonic oscillators, with length `chain_length=250` and local dimension `local_dim=6`,
prepared in the vacuum state.

# Returns
- `sysenv::Vector{<:Index}`: a vector that combines the system and environment indices;
- `psi0`: MPS representation of the initial state.
"""
function defineSystem(
    sys_state::String,
    chain_length::Int64,
    local_dim::Int64
)
    # Indexes needed for defining the MPS
    sys = siteinds("S=1/2", 1)
    env = siteinds("Boson", dim = local_dim, chain_length)
    sysenv = vcat(sys, env)
    # State of the system prepared in the `sys_state`
    stateSys = [sys_state]
    # Environment prepared in the vacuum state
    stateEnv = ["0" for n = 1:chain_length]
    # MPS representation
    stateSE = vcat(stateSys, stateEnv)
    psi0 = productMPS(ComplexF64, sysenv, stateSE)
    return (sysenv, psi0)
end

"""
    createMPO()

Create the MPO of the Hamiltonian H = Hs + H + H_I with Hamiltonian of interaction given by σ_k ⊗ ( b†{2} + b{2} )

Returns:
    - MPO(mpo,sysenv) i.e. MPO representation of the overall hamiltonian
"""
function createMPO(
    sysenv,
    epsilon::Float64,
    Delta::Float64,
    intHsysSide::String,
    freqfile::String,
    coupfile::String,
)::MPO
    # Read TEDOPA coefficients
    coups = readdlm(coupfile)
    freqs = readdlm(freqfile)
    # Retrieve chain size
    NN::Int64 = size(sysenv)[1]
    NChain::Int64 = NN - 1

    # From ITensors documentation:
    # - "S+" (alises: "S⁺", "Splus") Raising operator;
    # - "S-" (aliases: "S⁻", "Sminus") Lowering operator.
    # ATT! Sₓ = 0.5 σx, so whenever we need to use σᵢ we need to multiply by 2.

    # Start building the MPO:
    # H += constant, operator, site, operator, site, ...
    mpo = OpSum()

    # System Hamiltonian: H_S = ϵ σz + Δ σx
    mpo += 2 * epsilon, "Sz", 1
    mpo += 2 * Delta, "Sx", 1

    # Interaction Hamiltonian: H_I = As ⊗ κ₀(c₀† + c₀)    
    mpo += 2 * coups[1], intHsysSide, 1, "Adag", 2
    mpo += 2 * coups[1], intHsysSide, 1, "A", 2

    # Environment Hamiltonian is composed by two terms:
    #   - local oscillators: ∑ₙ ωₙ cₙ†cₙ, with n=0...NChain
    #   - couplings: ∑ₙ κₙ(cₙ†cₙ₋₁ + cₙ₋₁†cₙ), with n=1...NChain
    for j = 2:NChain
        mpo += freqs[j-1], "N", j
    end
    for j = 2:NChain-1
        mpo += coups[j], "A", j, "Adag", j + 1
        mpo += coups[j], "Adag", j, "A", j + 1
    end

    return MPO(mpo, sysenv)
end
