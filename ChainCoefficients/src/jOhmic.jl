"Given J(ω), compute the thermalized spectral density function."
function thermalize_sdf(x, sdf, β)
    return 0.5 * (1 + coth(0.5 * x * β)) * sign(x) * sdf(abs(x))
end

function read_sdf(dict::Dict{String, Any})
    a = Float64.(dict["environment"]["spectral_density_parameters"])
    # retrieves a string representing the spectral density function formula
    fn = dict["environment"]["spectral_density_function"]
    # creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    # defining a function sdf of one variable x. Calls the dynamically created function tmp
    sdf = x -> Base.invokelatest(tmp, a, x)
    return sdf
end

function read_thermalized_sdf(dict::Dict{String, Any})
    sdf = read_sdf(dict)
    β = dict["environment"]["β"]
    therm_sdf = x -> thermalize_sdf(x, sdf, β)
    return therm_sdf
end