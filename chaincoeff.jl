using Revise
using Pkg
Pkg.activate("./ChainCoefficients")

using JSON
using ChainCoefficients

# Save to directory
dir = "hc_0.5_2"
mkpath(dir)

s = 0.5                     # s=1 (Ohmic), s<1 (Sub-Ohmic), s>1 (Super-Ohmic)
ϵ, Δ = 0.2, 0.0             # Hs = ϵ/2 σz + Δ/2 σx, transition freq ωs=√ϵ²+Δ²
α = 0.01 / (π * sqrt(ϵ^2+Δ^2)^s)  # Rescale: α → α/(π ⋅ ωs^s) 

# I am interested into the following sdf (a = [α, s, ωc], see below)
hardcutOhmic = "2 * a[1] * a[3] * (x/a[3])^a[2]"
expcutOhmic = "pi/2 * a[1] * a[3] * (x/a[3])^a[2] * exp(-x/a[2])"

# Write here all the requirements to run the simulation
# Note: N = 2 * κ∞ * tmax, with κ∞ = ωmax / 2
p = Dict(
    "environment" => Dict(
        "spectral_density_parameters" => [α, s, 1.0],
        "spectral_density_function" => hardcutOhmic,
        "domain" => [0, 1],
        "β" => 2.0, # β=2.0 (T=0.5), β=2000.0 (T=0.0005)
    ),
    "chain_length" => 100,
    "nquad" => 5000
)
# Save configuration to a json file
open(joinpath(dir, "config_sdf.json"), "w") do io
    JSON.print(io, p, 4)
end

# Extract parameters
α, s, ωc = p["environment"]["spectral_density_parameters"]
sdf_func = p["environment"]["spectral_density_function"]
domain = p["environment"]["domain"]
β = p["environment"]["β"]
chain_length = p["chain_length"]
nquad = p["nquad"]

# Chain coefficients using `MPSDynamics.chaincoeffs_finiteT`
if sdf_func == hardcutOhmic
    freqs, coups = chaincoeff_jOhmic_hc(chain_length, β, α, s; nquad = nquad)
elseif sdf_func == expcutOhmic
    freqs, coups = chaincoeff_jOhmic_expc(chain_length, β, α, s; ωc=1.0, ωmax=Float64.(domain[2]), nquad = nquad)
else
    println("J(ω) not recognized")
end
write_chain_coefficients(dir, freqs, coups)

# Compute first chain coefficients by integral calculus
xmax = Float64.(domain[2])
my_freqs, my_coups = my_chain_coefficients(p, xmax)
println("\nCheck first coefficients with integral calculus")
println("Frequencies")
println("\tMe:          $my_freqs")
println("\tMPSDynamics: $(freqs[1:4])")
println()
println("Couplings")
println("\tMe:          $my_coups")
println("\tMPSDynamics: $(coups[1:4])")
