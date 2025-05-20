using Revise
using Pkg
Pkg.activate("./ChainCoefficients")

using DelimitedFiles
using JSON
using ChainCoefficients

# Specify the directory
dir = "hc_0.5_2000"
# Extend the chain length
extended_chain_length = 300
# Read the sdf
p = JSON.parsefile(joinpath(dir, "config_sdf.json"))

# I am interested into the following sdf (a = [α, s, ωc], see below)
hardcutOhmic = "2 * a[1] * a[3] * (x/a[3])^a[2]"
expcutOhmic = "pi/2 * a[1] * a[3] * (x/a[3])^a[2] * exp(-x/a[2])"

# Extract parameters
α, s, ωc = p["environment"]["spectral_density_parameters"]
sdf_func = p["environment"]["spectral_density_function"]
domain = p["environment"]["domain"]
β = p["environment"]["β"]
nquad = p["nquad"]

# Chain coefficients using `MPSDynamics.chaincoeffs_finiteT`
if sdf_func == hardcutOhmic
    freqs, coups = chaincoeff_jOhmic_hc(extended_chain_length, β, α, s; nquad = nquad)
elseif sdf_func == expcutOhmic
    freqs, coups = chaincoeff_jOhmic_expc(extended_chain_length, β, α, s; ωc=1.0, ωmax=10., nquad = nquad)
else
    println("J(ω) not recognized")
end

freqs_file = joinpath(dir, "freqs_$(extended_chain_length).csv")
coups_file = joinpath(dir, "coups_$(extended_chain_length).csv")
# Write frequencies
open(freqs_file, "w") do io
    for x in freqs
        writedlm(io, x, ',')
    end
end
# Write couplings
open(coups_file, "w") do io
    for x in coups
        writedlm(io, x, ',')
    end
end
println("Chain coefficients saved in directory: '$dir'")


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
