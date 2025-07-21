using Revise
using Pkg
Pkg.activate("QTools")

include("QTools/src/io.jl")
include("QTools//src/envmodes.jl")

include("TEBD/src/twoSiteMeas.jl")

using JSON
using DelimitedFiles


# ─────────────────────────────────────────────────────────────
# DATA
# ─────────────────────────────────────────────────────────────

dir = "expc_s1_beta50_d0.1_alpha0.05_t300_n1000_5"
state = "Up"

config = JSON.parsefile(joinpath(dir, "config.json"))
chainsize = config["chain"]["size"]


# ─────────────────────────────────────────────────────────────
# NORMAL MODES OCCUPATION
# ─────────────────────────────────────────────────────────────

# Read chain coefficients
freqs = readdlm(joinpath(dir, "freqs.dat"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(dir, "coups.dat"), '\n', Float64)[:, 1]

# Diagonalize tridiagonal matrix
modes, U = diagonalize_tridiagonal(freqs, coups)

# Save normal modes frequencies to .dat file
write_vector(dir, "envmodes_$state", modes)

# While reading chain site occupation measurements,
# compute normal modes occupation and save to .dat file

# Note: first column in 'rawData' is time
rawData = readdlm(joinpath(dir, "measOcc_$state.dat"), ',', Complex{Float64}, '\n')

# using 'indDoubleMeas' to list indexes for two site measurements
measIndices = indDoubleMeas(chainsize)


open(joinpath(dir, "envmodesOcc_$state.dat"), "w") do io

    for (t, row) in enumerate(eachrow(rawData))

        println("$t of $(size(rawData,1))")
        
        # Fill a n_sites × n_sites matrix with the chain sites occupation,
        # needed 'indDoubleMeas' to proper mapping from measurements file 
        measMatr = zeros(ComplexF64, chainsize, chainsize)
        for i = 1:length(measIndices)
            measMatr[measIndices[i][1]-1, measIndices[i][2]-1] = rawData[t, i+1]
            measMatr[measIndices[i][2]-1, measIndices[i][1]-1] = rawData[t, i+1]'
        end

        @assert measMatr[end, end] == 0. + 0. im "Excitation reaches the end"
        # Compute normal modes occupation and save to .dat file as a single row,
        # each row of file will be the normal modes occupation for each measurements time
        occ = envmodes_occupation(U, measMatr)
        writedlm(io, permutedims(real.(occ)), ',')
    end
end


# ─────────────────────────────────────────────────────────────
# NORMAL MODES OCCUPATION - LONGER CHAN
# ─────────────────────────────────────────────────────────────

dir = "expc_s1_beta50_d0.1_alpha0.03_t100_n500"
config = JSON.parsefile(joinpath(dir, "config.json"))
chainsize = config["chain"]["size"]
larger_chainsize = 800

# Read chain coefficients
freqs = readdlm(joinpath(dir, "freqs$larger_chainsize.dat"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(dir, "coups$larger_chainsize.dat"), '\n', Float64)[:, 1]

# Diagonalize tridiagonal matrix
modes, U = diagonalize_tridiagonal(freqs, coups)

# Save normal modes frequencies to .dat file
write_vector(dir, "envmodes$larger_chainsize", modes)

# While reading chain site occupation measurements,
# compute normal modes occupation and save to .dat file

# Note: first column in 'rawData' is time
rawData = readdlm(joinpath(dir, "measOcc_Up.dat"), ',', Complex{Float64}, '\n')

# using 'indDoubleMeas' to list indexes for two site measurements
measIndices = indDoubleMeas(chainsize)

open(joinpath(dir, "envmodesOcc$larger_chainsize.dat"), "w") do io

    for (t, row) in enumerate(eachrow(rawData))

        println("$t of $(size(rawData,1))")
        
        # Fill a n_sites × n_sites matrix with the chain sites occupation,
        # needed 'indDoubleMeas' to proper mapping from measurements file 
        measMatr = zeros(ComplexF64, larger_chainsize, larger_chainsize)
        for i = 1:length(measIndices)
            measMatr[measIndices[i][1]-1, measIndices[i][2]-1] = rawData[t, i+1]
            measMatr[measIndices[i][2]-1, measIndices[i][1]-1] = rawData[t, i+1]'
        end

        # Compute normal modes occupation and save to .dat file as a single row,
        # each row of file will be the normal modes occupation for each measurements time
        occ = envmodes_occupation(U, measMatr)
        writedlm(io, permutedims(real.(occ)), ',')
    end
end