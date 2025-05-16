module ChainCoefficients

using DelimitedFiles
using Integrals
using MPSDynamics

export chaincoeff_jOhmic_hc, chaincoeff_jOhmic_expc
export write_chain_coefficients
export jOhmic_hc, therm_jOhmic
export my_chain_coefficients

include("chaincoeffs.jl")
include("write.jl")
include("jOhmic.jl")
include("check.jl")

end
