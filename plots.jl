using Revise
using Pkg
Pkg.activate("./QEngine")
using DelimitedFiles
using JSON
using LinearAlgebra
using QEngine
Pkg.activate("./QPlot")
using QPlot

# Choose the directory with data
dir = "hc_2_2_ld_6"

# Parse configuration information
config_sdf = JSON.parsefile(joinpath(dir, "config_sdf.json"))
config_simul = JSON.parsefile(joinpath(dir, "config_simul.json"))

# Read chain coefficients
freqs = readdlm(joinpath(dir, "freqs.csv"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(dir, "coups.csv"), '\n', Float64)[:, 1]
# Read time vector
ts = readdlm(joinpath(dir, "time_Up.dat"), '\n', Float64)[:, 1]
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
plot_state(ts, measUp[:,1], measUp[:,2], measUp[:,3], "Up"; outdir=outdir)
plot_state(ts, measDown[:,1], measDown[:,2], measDown[:,3], "Down"; outdir=outdir)
plot_state(ts, measPlus[:,1], measPlus[:,2], measPlus[:,3], "Plus"; outdir=outdir)
plot_state(ts, measTrans[:,1], measTrans[:,2], measTrans[:,3], "Trans"; outdir=outdir)

# Effective Hamiltonian Ks
dt = config_simul["integration"]["time_step"]
Ks = computeKs(Up, Down, Plus, Trans, dt)
eigenvalues = [eigvals(K) for K in Ks]
Eminus = [v[1] for v in eigenvalues]
Eplus = [v[2] for v in eigenvalues]
θs = Eplus - Eminus
# write_Ks(dir, Ks)
plot_Ks(ts, Ks, 1, 1, β, s; outdir=outdir)
plot_Ks(ts, Ks, 1, 2, β, s; outdir=outdir)


# Chain occupation
animate_chain(ts, measN, β, s, joinpath("figs", dir))

# Normal modes occupation
# extended_chain = 300 # Consider a fictitiously enlarged chain
# freqs = readdlm(joinpath(dir, "freqs_$(extended_chain).csv"), '\n', Float64)[:, 1]
# coups = readdlm(joinpath(dir, "coups_$(extended_chain).csv"), '\n', Float64)[:, 1]
modes, occupations = envmodes_occupation(freqs, coups, measN)
animate_envmodes(ts, modes, occupations, J, ωs, β, s, joinpath("figs", dir))
# animate_envmodes(ts, modes, occupations, J, ωs, θs,  β, s, joinpath("figs", dir))

