using Revise
using Pkg
Pkg.activate("./QEngine")

using DelimitedFiles
using JSON
using QEngine

# Choose the directory with data
dir = "hc_0.5_2000"

# Parse configuration information
config_sdf = JSON.parsefile(joinpath(dir, "config_sdf.json"))
config_simul = JSON.parsefile(joinpath(dir, "config_simul.json"))

# Read chain coefficients
freqs = readdlm(joinpath(dir, "freqs.csv"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(dir, "coups.csv"), '\n', Float64)[:, 1]
# Read time vector
time = readdlm(joinpath(dir, "time.dat"), '\n', Float64)[:, 1]
# Read bloch vectors
measUp = read_state(dir, "Up")
measDown = read_state(dir, "Dn")
measPlus = read_state(dir, "+")
measTrans = read_state(dir, "i")
# ... and corresponding density matrices
Up = densitymatrix.(measUp[:,1], measUp[:,2], measUp[:,3])
Down = densitymatrix.(measDown[:,1], measDown[:,2], measDown[:,3])
Plus = densitymatrix.(measPlus[:,1], measPlus[:,2], measPlus[:,3])
Trans = densitymatrix.(measTrans[:,1], measTrans[:,2], measTrans[:,3])
# Read measurements of chain sites occupation
measN = read_measN(dir)

# Now plots!
mkpath(joinpath("figs", dir))
outdir = joinpath("figs", dir) # outdir = nothing
# Parse some information
β = config_sdf["environment"]["β"]
α, s, ωc = config_sdf["environment"]["spectral_density_parameters"]
ϵ, Δ = config_simul["qubit"]["epsilon"], config_simul["qubit"]["delta"]
ωs = sqrt(ϵ^2+Δ^2) # Frequency transition of Hs
J = read_thermalized_sdf(config_sdf)

# State
plot_state(time, measUp[:,1], measUp[:,2], measUp[:,3], "Up"; outdir=outdir)
plot_state(time, measDown[:,1], measDown[:,2], measDown[:,3], "Down"; outdir=outdir)
plot_state(time, measPlus[:,1], measPlus[:,2], measPlus[:,3], "Plus"; outdir=outdir)
plot_state(time, measTrans[:,1], measTrans[:,2], measTrans[:,3], "Trans"; outdir=outdir)

# Effective Hamiltonian Ks
dt = config_simul["integration"]["time_step"]
Ks = computeKs(Up, Down, Plus, Trans, dt)
# write_Ks(dir, Ks)
plot_Ks(time, Ks, 1, 1, β, s; outdir=outdir)
plot_Ks(time, Ks, 1, 2, β, s; outdir=outdir)

# Chain occupation
animate_chain(time, measN, β, s, joinpath("figs", dir))

# Normal modes occupation
extended_chain = 300 # Consider a fictitiously enelarged chain
freqs = readdlm(joinpath(dir, "freqs_$(extended_chain).csv"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(dir, "coups_$(extended_chain).csv"), '\n', Float64)[:, 1]
modes, occupations = envmodes_occupation(freqs, coups, measN)
animate_envmodes(time, modes, occupations, J, ωs, β, s, joinpath("figs", dir))


