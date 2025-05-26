using Revise
using Pkg
Pkg.activate("C:/Users/andre/Desktop/master-thesis/QEngine")

using LinearAlgebra
using DelimitedFiles
using QEngine

basedir = @__DIR__  # folder where the script is located

# Read chain coefficients
freqs = readdlm(joinpath(basedir, "freqs.dat"), '\n', Float64)[:, 1]
coups = readdlm(joinpath(basedir, "coups.dat"), '\n', Float64)[:, 1]

# Me
measN = real.(readdlm(joinpath(basedir, "popMeas.dat"), ',', Complex{Float64}, '\n'))[:,2:end]
modes, me_occu = envmodes_occupation(freqs, coups, measN)

# Tamascelli
datiPop = readdlm(joinpath(basedir, "popMeas.dat"), ',', Complex{Float64}, '\n')
popMatrix = real(datiPop)
ham = Array(Tridiagonal(coups[2:end],freqs,coups[2:end]))
eigensystem = eigen(ham)
uu = eigensystem.vectors
normalModesPop=zeros(Float64,size(popMatrix)...)
n_rows = size(popMatrix, 1)
n_cols = size(popMatrix, 2)
n_modes = size(uu, 2)
for i in 1:n_rows
    normalModesPop[i, 1] = popMatrix[i, 1]
    for j in 2:n_modes
        appo = map(x -> x * x, transpose(uu[:, j]))
        normalModesPop[i, j + 1] = appo * popMatrix[i, 2:end]
    end
end
nonme_occu = normalModesPop[:,2:end]

me_occu ≈ nonme_occu

diff = abs.(me_occu - nonme_occu)
maximum(diff)