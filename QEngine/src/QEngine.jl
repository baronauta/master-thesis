module QEngine

using DelimitedFiles    # write and read files
using JSON              # read JSON file
using LinearAlgebra
using CairoMakie        # plot
using GLMakie           # animated plot

export σ0, σ1, σ2, σ3
export paulimatrices, canonical_op_basis
export densitymatrix, computeKs, write_Ks
export envmodes_occupation, write_envmodes

export read_state, read_measN
export read_modes, read_occupations
export read_Ks
export read_sdf, read_thermalized_sdf

export plot_state
export plot_Ks
export animate_chain, animate_envmodes

export read_thermalized_sdf

include("constants.jl")
include("read.jl")
include("tls_dynamics.jl")
include("env_dynamics.jl")
include("tls_plot.jl")
include("env_plot.jl")

end
