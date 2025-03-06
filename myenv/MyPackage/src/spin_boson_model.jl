function defineSystem(;
    sys_type::String="S=1/2",sys_istate::String="Up",
    chain_size::Int64=250,local_dim::Int64=6
)

"""
Standard approach: initial sysenv state factorized, chain always in the vacuum state.
By dafault consider:
    - TLS system (sys_type::String="S=1/2") in the excited state (sys_istate::String="Up")
    - chain of "Boson" with length `chain_size` and local dimension `local_dim`

Returns:
    - `sysenv` vector that combines the system and the environment indices
    - `psi0` MPS representation of the initial state
"""

sys = siteinds(sys_type,1)
env = siteinds("Boson", dim=local_dim, chain_size)
sysenv = vcat(sys,env)

stateSys = [sys_istate]

stateEnv = ["0" for n=1:chain_size] # chain in the vacuum state

stateSE = vcat(stateSys,stateEnv)

psi0 = productMPS(ComplexF64,sysenv,stateSE)

return (sysenv,psi0)
end


function createMPO(
    sysenv, 
    eps::Float64, delta::Float64,
    intHsysSide::String,
    freqfile::String, coupfile::String
    )::MPO

"""
Create the MPO of the overall hamiltonian with hamiltonian of interaction given by S_k ⊗ ( b†{2} + b{2} )

Returns:
    - MPO(mpo,sysenv) i.e. MPO representation of the overall hamiltonian
"""

# to check what interaction I am considering calling this function
println("system-bath interaction: ", intHsysSide, " ⊗ ( b†{2} + b{2} )")

# read TEDOPA coefficients
coups = readdlm(coupfile)
freqs = readdlm(freqfile)

NN::Int64 = size(sysenv)[1]
NChain::Int64 = NN-1

# From ITensors documentation
# - "S+" (alises: "S⁺", "Splus") Raising operator 
# - "S-" (aliases: "S⁻", "Sminus") Lowering operator
# ATT! S_x = 0.5 σx, so whenever we need to use σ_i we need to multiply by 2

# start building up the MPO
# H += constant, operator, site, operator, site, ...
mpo = OpSum()

# system hamiltonian
mpo += 2*eps,"Sz",1
mpo += 2*delta,"Sx",1

# system-env interaction    
mpo += 2*coups[1],intHsysSide,1,"Adag",2
mpo += 2*coups[1],intHsysSide,1,"A",2

# chain local hamiltonians
for j=2:NChain
  mpo += freqs[j-1],"N",j
end

# chain coupling hamiltonians
for j=2:NChain-1
  mpo += coups[j],"A",j,"Adag",j+1
  mpo += coups[j],"Adag",j,"A",j+1
end

return MPO(mpo,sysenv)
end
