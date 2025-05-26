using Pkg
Pkg.activate(".")

using DelimitedFiles
using MPSDynamics
using SpinBosonDynamics

dir = "SB_hc_2_2000"
mkpath(dir)

# Spin-Boson Hamiltonian (see https://shareloqs.github.io/MPSDynamics.jl/methods/#MPSDynamics.spinbosonmpo-NTuple{5,%20Any})

d = 6               # number of Fock states of the chain modes
N = 100             # length of the chain

Δ = 0.0             # tunneling 
ω0 = 0.2            # TLS gap

s = 2               # ohmicity
α = 0.01 / ω0^2     # coupling strength
β = 2000.0          # inverse temperature

# Ohmic spectral density function `J(ω) = 2 α ωc (ω/ωc)^s θ(ω-ωc)`.
cpars = chaincoeffs_finiteT(N, β, true; α = α, s = s, ωc = 1, save = false)
# Chain coefficients: [frequencies, couplings, [κ₀]].
freqs = cpars[1]
coups = vcat(cpars[3][1], cpars[2])
write_chain_coefficients(dir, freqs, coups)

# Spin-Boson Hamiltonian as MPO
H = spinbosonmpo(ω0, Δ, d, N, cpars)

# System in the Up state and environment in the vacuum state
ψ = unitcol(1, 2)
A = productstatemps(physdims(H), state = [ψ, fill(unitcol(1, d), N)...])

# Observables to measure
ob1 = OneSiteObservable("sz", sz, 1)
ob2 = OneSiteObservable("chain mode occupation", numb(d), (2, N + 1))

# Run the simulation

dt = 0.05
T = 50.0
prec1 = 1e-10
A, dat = runsim(
    dt,
    T,
    A,
    H;
    name = "ohmic spin boson model",
    method = :DTDVP,
    obs = [ob2],
    convobs = [ob1],
    params = @LogParams(N, d, α, Δ, ω0, s),
    convparams = [prec1],
    verbose = false,
    save = false
)

measUp = dat["data/sz"]
measN = dat["data/chain mode occupation"]   # site x time
time = dat["data/times"]

open(joinpath(dir, "measUp.dat"), "w") do io
    for x in measUp
        writedlm(io, x, ',')
    end
end

measN = transpose(measN) # time x site
open(joinpath(dir, "measN.dat"), "w") do io
    writedlm(io, measN, ',')
end

open(joinpath(dir, "time.dat"), "w") do io
    for x in time
        writedlm(io, x, ',')
    end
end