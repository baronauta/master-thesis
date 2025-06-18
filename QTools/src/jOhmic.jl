import SpecialFunctions: gamma
import JSON

"""
    myJ(ω; α=1, s=1, ωc=1.0)

Spectral density function from the Ohmic family:

- Peaks at ω = ωc for any s
- Normalized to keep ∫J(ω)dω constant across s

# Arguments
- `ω`: Frequency
- `α`: Coupling strength
- `s`: Ohmicity (1 = Ohmic, <1 sub-Ohmic, >1 super-Ohmic)
- `ωc`: Characteristic frequency
"""
function myJ(ω; α=1, s=1, ωc=1.0)
    κ = s^(s+1) / gamma(s+1)
    return π /2. * α * κ * ωc * (ω / ωc)^s * exp(- s * ω / ωc)
end

"Given J(ω), compute the thermalized spectral density function."
function thermalize_sdf(x, sdf, β)
    return 0.5 * (1 + coth(0.5 * x * β)) * sign(x) * sdf(abs(x))
end


# --- Read J(ω) from dictionary ---

"Read J(ω) from dictionary and return a callable function."
function read_sdf(dict::Dict{String,Any})
    a = Float64.(dict["environment"]["spectral_density_parameters"])
    # Retrieves a string representing the spectral density function formula
    fn = dict["environment"]["spectral_density_function"]
    # Creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    # Defining a function sdf of one variable x. Calls the dynamically created function tmp
    sdf = x -> Base.invokelatest(tmp, a, x)
    return sdf
end

"Read J(ω) from dictionary and return a callable function for the thermalized J(ω,β)."
function read_thermalized_sdf(dict::Dict{String,Any})
    sdf = read_sdf(dict)
    β = dict["environment"]["β"]
    therm_sdf = x -> thermalize_sdf(x, sdf, β)
    return therm_sdf
end