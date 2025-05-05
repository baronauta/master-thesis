# ─────────────────────────────────────────────────────────────
# Spectral density function
# - sdf
# - thermalized sdf
# - chain coefficients
# ─────────────────────────────────────────────────────────────

"""
    sdf(filename::String)

Compute the spectral density J(ω) over its domain from JSON parameters in `filename`.
"""
function sdf(filename::String)
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    a = Float64.(p["environment"]["spectral_density_parameters"])
    support = Float64.(p["environment"]["domain"])
    fn = p["environment"]["spectral_density_function"]
    # Creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    f = x -> Base.invokelatest(tmp, a, x)
    xs = collect(range(support..., 1000))
    ys = f.(xs)
    return xs, ys
end

function boson_thermal_factor(x, T)
    if T == 0
        # In this case we could write f=1 and be done with it, but instead we choose
        # this more articulated way so that even if the caller doesn't already
        # exclude (-∞,0) from the support, we do it ourselves now.
        f = x > 0 ? one(x) : zero(x)
    else
        f = 1 / 2 * (1 + coth(0.5 * x / T))
    end
    return f
end

"""
    thermal_sdf(filename::String)

Construct the thermal spectral density J_β(ω) at temperature T from JSON in `filename`.
"""
function thermal_sdf(filename::String)
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    a = Float64.(p["environment"]["spectral_density_parameters"])
    T = Float64(p["environment"]["temperature"])
    support = Float64.(p["environment"]["domain"])
    extended_support = (-support[2], support[2])
    fn = p["environment"]["spectral_density_function"]
    # Creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, a, x)
    # Defining the extended function
    thermal_sdf(x) = boson_thermal_factor(x, T) * sign(x) * sdf(abs(x))
    xs = collect(range(extended_support..., 1000000))
    ys = thermal_sdf.(xs)
    return xs, ys, a, T
end

"""
    plot_sdf(filename::String)

Plot the spectral density J(ω) from JSON file `filename`.
Returns a `Figure`.
"""
function plot_sdf(filename::String)
    xs, ys = sdf(filename)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys)
    f
end

"""
    plot_thermal_sdf(filename::String) -> Figure

Plot the thermal spectral density J_β(ω) with temperature and parameters annotated.
Returns a `Figure`.
"""
function plot_thermal_sdf(filename::String)
    xs, ys, a, temp = thermal_sdf(filename)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J_\beta(\omega)")
    lines!(ax, xs, ys, label = L"T=%$(temp),\,\alpha=%$(a[1]),\,s=%$(a[3])")
    axislegend(position = :rt)
    f
end

"""
    plot_chain_coefficients(filename::String)

Plot chain transformation frequencies and couplings from `filename` data.
Returns a `Figure`.
"""
function plot_chain_coefficients(filename::String)
    coefficients = chain_coefficients(filename)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"n", ylabel = "Chain coefficients")
    lines!(
        ax,
        collect(1:length(coefficients.frequencies)),
        coefficients.frequencies,
        label = "Freqs",
    )
    lines!(
        ax,
        collect(1:length(coefficients.couplings)),
        coefficients.couplings,
        label = "Coups",
    )
    axislegend(position = :rb)
    f
end
