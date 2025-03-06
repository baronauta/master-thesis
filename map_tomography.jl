#!/usr/bin/env julia
using Pkg
Pkg.activate("./myenv")
using MyPackage

# Define sdf parameters
sdf_type = "debye"
a = 0.1
T = 0.0
sdf = "./sdf/$(sdf_type)_a_" * string(a) * "_T_" * string(T) * ".json"

# Define TLS parameters
ϵ = 0.0
Δ = 0.2

map_tomography(sdf, ϵ, Δ; dt=0.005, tmax=30, growMPSval=10, local_dim=6)
