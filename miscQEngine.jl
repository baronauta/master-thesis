using Revise
using Pkg
Pkg.activate("./QEngine")

using QEngine

# Choose the directory with data
dir = "hc_2_2000"

# Parse configuration information
config_sdf = read_config_sdf(dir)
config_simul = read_config_simul(dir)

# Read chain coefficients
freqs = read_freqs(dir)
coups = read_coups(dir)
# Read time vector
time = read_time(dir)
# Read bloch vectors
measUp = read_state(dir, "Up")
measDown = read_state(dir, "Dn")
measPlus = read_state(dir, "+")
measTrans = read_state(dir, "i")
# Read measurements of chain sites occupation
measN = read_measN(dir)

# Compute effective Hamiltonian Ks
Up = densitymatrix.(measUp[:,1], measUp[:,2], measUp[:,3])
Down = densitymatrix.(measDown[:,1], measDown[:,2], measDown[:,3])
Plus = densitymatrix.(measPlus[:,1], measPlus[:,2], measPlus[:,3])
Trans = densitymatrix.(measTrans[:,1], measTrans[:,2], measTrans[:,3])
dt = config_simul["integration"]["time_step"]
Ks = computeKs(Up, Down, Plus, Trans, dt)
write_Ks(dir, Ks)

# Compute normal modes and theirs occupation
modes, occupations = envmodes_occupation(freqs, coups, measN)
write_envmodes(dir, modes, occupations)

# Now plots!
mkpath(joinpath("figs", dir))
outdir = joinpath("figs", dir)
outdir = nothing
β = config_sdf["environment"]["β"]
α, s, ωc = config_sdf["environment"]["spectral_density_parameters"]
# State
plot_state(time, measUp[:,1], measUp[:,2], measUp[:,3], "Up"; outdir=outdir)
plot_state(time, measDown[:,1], measDown[:,2], measDown[:,3], "Down"; outdir=outdir)
plot_state(time, measPlus[:,1], measPlus[:,2], measPlus[:,3], "Plus"; outdir=outdir)
plot_state(time, measTrans[:,1], measTrans[:,2], measTrans[:,3], "Trans"; outdir=outdir)
# Ks
myKs = read_Ks(dir)
plot_Ks(time, myKs, 1, 1, β, s; outdir=outdir)
plot_Ks(time, myKs, 1, 2, β, s; outdir=outdir)
# Chain occupation and normal mode occupations
animate_chain(time, measN, β, s, joinpath("figs", dir))
animate_envmodes(time, modes, occupations, β, s, joinpath("figs", dir))


