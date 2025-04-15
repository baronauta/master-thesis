function plot_sdf(sdf_filename::String)
    input = open(sdf_filename)
    s = read(input, String)
    p = JSON.parse(s)
    a = Float64.(p["environment"]["spectral_density_parameters"])
    support = Float64.(p["environment"]["domain"])
    fn = p["environment"]["spectral_density_function"]
    # Creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, a, x)
    xs = collect(range(support..., 1000))
    ys = sdf.(xs)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label = L"\alpha=%$(a[1])")
    axislegend(position = :rt)
    f
end

function plot_extended_sdf(sdf_filename::String)
    p = open(sdf_filename) do input
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
    sdf = x -> Base.invokelatest(tmp, a, x)
    # Defining the extended function
    extended_sdf = x -> (x > 0 ? sdf(x) : -sdf(-x)) * (1 / 2) * (1 + coth(x / (T * 2.0)))
    xs = collect(range(extended_support..., 1000))
    ys = extended_sdf.(xs)
    # Create and display the plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"J_\beta(\omega)")
    lines!(ax, xs, ys, label = L"\alpha=%$(a[1]),\,T=%$(T)")
    axislegend(position = :rt)
    f
end

function plot_chain_coefficients(sdf_filename::String)
    coefficients = chain_coefficients(sdf_filename)
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
