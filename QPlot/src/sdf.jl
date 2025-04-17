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

function thermal_sdf(filename::String)
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    a = Float64.(p["environment"]["spectral_density_parameters"])
    T = Float64(p["environment"]["temperature"])
    if T == 0
        println("ERROR: environment at zero temperature.")
        return
    end
    support = Float64.(p["environment"]["domain"])
    extended_support = (-support[2], support[2])
    fn = p["environment"]["spectral_density_function"]
    # Creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    f = x -> Base.invokelatest(tmp, a, x)
    # Defining the extended function
    thermal_f = x -> (1 / 2) * (1 + coth(0.5 * x / T)) * sign(x) * f(abs(x))
    xs = collect(range(extended_support..., 1000))
    ys = thermal_f.(xs)
    return xs, ys, T
end

function plot_sdf(filename::String)
    xs, ys = sdf(filename)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys)
    f
end

function plot_thermal_sdf(filename::String)
    xs, ys, temp = thermal_sdf(filename)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J_\beta(\omega)")
    lines!(ax, xs, ys, label = L"T=%$(temp)")
    axislegend(position = :rt)
    f
end

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
