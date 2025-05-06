# ─────────────────────────────────────────────────────────────
# Spectral density function
# - sdf
# - thermalized sdf
# - chain coefficients
# ─────────────────────────────────────────────────────────────

"""
    plot_sdf(filename::String)

Plot the spectral density J(ω) from JSON file `filename`.
Returns a `Figure`.
"""
function plot_sdf(filename::String)
    # Parse filename with sdf paramters
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    # Read a and support
    a = Float64.(p["environment"]["spectral_density_parameters"])
    support = Float64.(p["environment"]["domain"])
    # Create and a callable object to compute sdf
    fn = p["environment"]["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, a, x)

    f = Figure()
    xs = collect(range(support..., 1000))
    ys = sdf.(xs)
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys)
    f
end

"Copied from pharrex repository."
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
    plot_thermal_sdf(filename::String) -> Figure

Plot the thermal spectral density J_β(ω) with temperature and parameters annotated.
Returns a `Figure`.
"""
function plot_thermal_sdf(filename::String)
    # Parse filename with sdf paramters
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    # Read a and support
    a = Float64.(p["environment"]["spectral_density_parameters"])
    support = Float64.(p["environment"]["domain"])
    therm_support = [-support[2], support[2]]
    # Create a callable object to compute sdf
    fn = p["environment"]["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, a, x)
    # Defining the thermalized sdf
    therm_sdf(x) = boson_thermal_factor(x, T) * sign(x) * sdf(abs(x))

    f = Figure()
    xs = collect(range(therm_support..., 1000000))
    ys = therm_sdf.(xs)
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
