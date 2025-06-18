using Revise
using Pkg
Pkg.activate("QTools")

include("QTools/src/io.jl")
include("QTools/src/jOhmic.jl")
include("QTools/src/chaincoeffs.jl")

import SpecialFunctions: gamma

using JSON
using MPSDynamics


# ─────────────────────────────────────────────────────────────
# SIMULATION SETUP
# ─────────────────────────────────────────────────────────────

dir = "expc_s1_beta2000"
mkpath(dir)


# qubit: Hs = ϵ⋅σz + Δ⋅σx
ϵ = 0.0
Δ = 0.2
qubit_params = Dict("epsilon" => ϵ, "delta" => Δ)


# J(ω) = π/2 * α * κ * ωc * (ω/ωc)^s * exp[-s*(ω/ωc)]
# κ = [s^(s+1)] / Γ(s+1)

const myOhmic = "pi/2 * a[1] * a[4] * a[3] * (x/a[3])^a[2] * exp( -a[2] * x/a[3] )"

# a[1]: α, a[2]: s, a[3] = ωc, a[4]: κ
α = 0.1
s = 1
ωc = 1
κ = s^(s+1) / gamma(s+1)

# domain: [0, ωmax]
ωmax = 10

# temperature: β=2.0 (T=0.5), β=2000.0 (T=0.0005)
β = 2000

environment_params = Dict(
    "spectral_density_parameters" => [α, s, ωc, κ],
    "spectral_density_function" => myOhmic,
    "domain" => [0, ωmax],
    "β" => β,
    "nquad" => 10000,
)


# chain and integration:
#   - N = 2 * κ∞ * tmax = ωmax * tmax
#   - dt = 1 / κ∞ * 50 = 1 / ωmax * 25
# Note that κ∞ = ωmax / 2

tmax = 50

κ∞ = ωmax / 2
N = Int(ωmax * tmax)
dt = 1. / (ωmax * 25)

chain_params = Dict("size" => N, "local_dim" => 6)

integration_params = Dict(
    "time_step" => dt,
    "tmax" => tmax,
    "measurement_step" => 1,
    "env_measurement_step" => 200,
)

mps_params = Dict(
    "minBondDim" => 5,
    "maxBondDim" => 100,
    "cutoff" => 1e-14
)


# Store all into a JSON file
p = Dict(
    "qubit" => qubit_params,
    "chain" => chain_params,
    "environment" => environment_params,
    "integration" => integration_params,
    "interaction" => "Z",
    "mps" => mps_params,
    "notes" => Dict(
    "H_system" => "Hs = ε * Z + Δ * X ",
    "H_interaction" => "Z ⊗ ∑ (b† + b)",
)
)

open(joinpath(dir, "config.json"), "w") do io
    JSON.print(io, p, 4)
end


# ─────────────────────────────────────────────────────────────
# CHAIN COEFFICIENTS (via MPSDynamics)
# ─────────────────────────────────────────────────────────────

α, s, ωc, κ = p["environment"]["spectral_density_parameters"]
sdf_string = p["environment"]["spectral_density_function"]
ωmax = Float64(p["environment"]["domain"][2])
β = p["environment"]["β"]
chain_length = p["chain"]["size"]
nquad = p["environment"]["nquad"]

# using `MPSDynamics.chaincoeffs_finiteT`

# Spectral density base function (positive frequencies)
sdf(x) = (π / 2) * α * κ * ωc .* (x ./ ωc) .^ s .* exp.(- s.* x ./ ωc)
# Full spectral density function over two intervals (negative and positive frequencies)
function J(x, i)
    if i == 1
        return 0
    elseif i == 2
        return -0.5 .* (1 .+ coth.(0.5 .* x .* β)) .* sdf.(abs.(x))
    elseif i == 3
        return 0.5 .* (1 .+ coth.(0.5 .* x .* β)) .* sdf.(abs.(x))
    elseif i == 4
        return 0
    else
        println("Invalid interval index i = $i. Expected 1, 2, 3 or 4.")
    end
end

# Compute chain coefficients: [frequencies, couplings, [κ₀]]
println("Begin calculation of the chain coefficients")
coeff = chaincoeffs_finiteT(
    chain_length,
    β,
    false;
    J = J,                        # custom J(ω)
    mc = 4,                       # intervals
    AB = [[-Inf -ωmax]; [-ωmax 0]; [0 ωmax]; [ωmax Inf]],    # defines the domain of J
    Mmax = nquad,
    save = false,
)
println("Calculation finished")

freqs = coeff[1]
coups = vcat(coeff[3][1], coeff[2])

write_vector(dir, "freqs", freqs)
write_vector(dir, "coups", coups)

# Check correctness of the first chain coefficients by integral calculus
my_freqs, my_coups = check_chaincoeff(p, ωmax)
println("\nCheck first coefficients with integral calculus")
println("Frequencies")
println("\tMe:          $my_freqs")
println("\tMPSDynamics: $(freqs[1:4])")
println()
println("Couplings")
println("\tMe:          $my_coups")
println("\tMPSDynamics: $(coups[1:4])")
