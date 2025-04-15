"""
    defineSystem(; sys_type::String="S=1/2", sys_istate::String="Up", chain_length::Int64=250, local_dim::Int64=6)

Standard approach: the initial system-environment state is factorized, and the chain is always in the vacuum state.

By default, this function considers:
- A two-level system (TLS) (`sys_type="S=1/2"`) in the excited state (`sys_istate="Up"`).
- A chain of "Boson" with length `chain_length=250` and local dimension `local_dim=6`.

# Returns
- `sysenv::Vector{<:Index}`: A vector that combines the system and environment indices.
- `psi0`: An MPS (Matrix Product State) representation of the initial state.
"""
function defineSystem(;
    sys_type::String = "S=1/2",
    sys_istate::String = "Up",
    chain_length::Int64 = 250,
    local_dim::Int64 = 6,
)
    sys = siteinds(sys_type, 1)
    env = siteinds("Boson", dim = local_dim, chain_length)
    sysenv = vcat(sys, env)
    stateSys = [sys_istate]
    stateEnv = ["0" for n = 1:chain_length] # chain in the vacuum state
    stateSE = vcat(stateSys, stateEnv)
    psi0 = productMPS(ComplexF64, sysenv, stateSE)
    return (sysenv, psi0)
end

"""
Create the MPO of the Hamiltonian H = H_S + H_E + H_I with Hamiltonian of interaction given by S_k ⊗ ( b†{2} + b{2} )

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

    # From ITensors documentation
    # - "S+" (alises: "S⁺", "Splus") Raising operator 
    # - "S-" (aliases: "S⁻", "Sminus") Lowering operator
    # ATT! S_x = 0.5 σx, so whenever we need to use σ_i we need to multiply by 2

    # start building up the MPO
    # H += constant, operator, site, operator, site, ...
    mpo = OpSum()

    # system hamiltonian
    mpo += 2 * epsilon, "Sz", 1
    mpo += 2 * Delta, "Sx", 1

    # system-env interaction    
    mpo += 2 * coups[1], intHsysSide, 1, "Adag", 2
    mpo += 2 * coups[1], intHsysSide, 1, "A", 2

    # chain local hamiltonians
    for j = 2:NChain
        mpo += freqs[j-1], "N", j
    end

    # chain coupling hamiltonians
    for j = 2:NChain-1
        mpo += coups[j], "A", j, "Adag", j + 1
        mpo += coups[j], "Adag", j, "A", j + 1
    end

    return MPO(mpo, sysenv)
end
