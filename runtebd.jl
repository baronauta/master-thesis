# using Revise
using Pkg
Pkg.activate("./TEBD")

using MKL
using JSON
using DelimitedFiles
using TEBD

# Save to directory
dir = "hc_2_2000"

state = "Up"         # Choose the initial state for the qubit
trackN = true        # Track env occupation
trackBondDim = true  # Track bond dimension

# Write here all the requirements to run the simulation
p = Dict(
    "qubit" => Dict(
        "epsilon" => 0.2,
        "delta" => 0.0
    ),
    "chain" => Dict(
        "size" => 100,          # N = 2 * κ∞ * tmax, with κ∞ = ωmax / 2
        "local_dim" => 6
    ),
    "integration" => Dict(
        "time_step" => 1 / 25,  # dt = 1 / (κ∞ * 50) = 1 / ωmax * 25
        "tmax" => 100,
        "measurement_step" => 1
    ),
    "mps" => Dict(
        "min_bond_dim" => 10,
        "max_bond_dim" => 120,
        "cutoff" => 1e-15
    ),
    "interaction" => "X"
)
# Save configuration to json file
open(joinpath(dir, "config_simul.json"), "w") do io
    JSON.print(io, p, 4)
end

# Extract parameters
ϵ = p["qubit"]["epsilon"]
Δ = p["qubit"]["delta"]
chain_size = p["chain"]["size"]
local_dim = p["chain"]["local_dim"]
τ = p["integration"]["time_step"]
tmax = p["integration"]["tmax"]
mStep = p["integration"]["measurement_step"]
minBondDim = p["mps"]["min_bond_dim"]
maxBondDim = p["mps"]["max_bond_dim"]
cutoff = p["mps"]["cutoff"]
sysenvInt = p["interaction"]

# Load chain coefficients
freqs = readdlm(joinpath(dir, "freqs.csv"), '\n', Float64)
coups = readdlm(joinpath(dir, "coups.csv"), '\n', Float64)
println("Chain coefficients loaded: $(length(freqs)) freqs and $(length(coups)) coups")

println("\nAvailable threads: $(Threads.nthreads())")

# Prepare state
sysenv, psi0 = prepareState(chain_size, local_dim; state = state)
psi = psi0
Lambdas, Gammas = convert_to_Vidal(psi)
println("Conversion to Vidal done: $(length(Lambdas)) Lambdas and $(length(Gammas)) Gammas")

# Run simulation
SpinBoson_evolution_TEBD(
    Gammas,
    Lambdas,
    sysenv;
    ϵ = ϵ,
    Δ = Δ,
    sysenvInt = sysenvInt,
    ChainLength = chain_size,
    tau = τ,
    ttotal = tmax,
    measStep = mStep,
    freqs = freqs,
    coups = coups,
    minBondDim = minBondDim,
    maxBondDim = maxBondDim,
    cutoff = cutoff,
    trackBondDim = trackBondDim,
    trackN = trackN,
    state = state,
    dir = dir
);