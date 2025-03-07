function load_sdf(sdf_filename)
    # Load spectral density function parameters
    input = open(sdf_filename)
    s = read(input, String)
    sdf_dict = JSON.parse(s)
    return sdf_dict
end


function sdf_params(sdf_filename)
    sdf_dict = load_sdf(sdf_filename)
    α, ωc = Float64.(sdf_dict["environment"]["spectral_density_parameters"])
    T = Float64(sdf_dict["environment"]["temperature"])
    sdf_eq = String(sdf_dict["environment"]["spectral_density_function"])
    chain_size = Int64(sdf_dict["chain_length"])
    return α, ωc, T, sdf_eq, chain_size
end


function sdf_naming(sdf_eq::String)
    if sdf_eq == "pi/2* a[1] * x * exp(-x/a[2])"
        return "ohmic"
    elseif sdf_eq == "pi/2 * a[1] * x * a[2]^2 / ( x^2 + a[2]^2 )"
        return "debye"
    else
        return "unknown"
    end
end


function save_tedopa_coefficients(coefficients, sdf_filename)
    # Store coefficients inside directory sdf/chain_coefficients
    mkpath("./sdf/chain_coefficients")
    # Extract base filename without directory and extension
    base_filename = splitext(basename(sdf_filename))[1]
    freqs_file = "./sdf/chain_coefficients/$(base_filename)_freqs.csv"
    coups_file = "./sdf/chain_coefficients/$(base_filename)_coups.csv"
    # frequencies
    io = open(freqs_file, "w")
    for x in coefficients.frequencies
        writedlm(io, x, ',')
    end
    close(io)
    # couplings
    io = open(coups_file, "w")
    for x in coefficients.couplings
        writedlm(io, x, ',')
    end
    close(io)
    return freqs_file, coups_file
end


function figs_sdf(sdf_filename)
    # open and read the JSON file
    input = open(sdf_filename)
    s = read(input, String)
    # parsing JSON data: parse the string s as a JSON object and assigns the resulting dictionary-like structure to p
    p = JSON.parse(s)
    # extracting parameters
    a = Float64.(p["environment"]["spectral_density_parameters"])
    support = Float64.(p["environment"]["domain"])
    # retrieves a string representing the spectral density function formula
    fn = p["environment"]["spectral_density_function"]
    # creating a callable function object tmp
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    # defining a function sdf of one variable x. Calls the dynamically created function tmp
    sdf = x -> Base.invokelatest(tmp, a, x)
    # creating an ohmic spectral density function
    jOhmic(ω) = sdf(ω)
    # plot the spectral density function
    xs = collect(range(support..., 1000))
    ys = jOhmic.(xs)
    # Create and display the plot
    p = plot(
        xs,
        ys,
        label = L"\alpha=%$(a[1]),\,\omega_c=%$(a[2])",
        xlabel = L"\omega",
        ylabel = L"J(\omega)",
        xlabelfontsize = 16,
        ylabelfontsize = 16,
        grid = false,
    )
    mkpath("./sdf/figs")
    savefig("./sdf/figs/$(sdf_type)_a_" * string(a) * "_T_" * string(T) * ".json")
end


function figs_tedopa_coefficients(sdf_filename)
    sdf_dict = load_sdf(sdf_filename)
    T = Float64(sdf_dict["environment"]["temperature"])
    # TEDOPA or T-TEDOPA according to temperature T
    coefficients =
        T == 0 ? chainmapping_tedopa(sdf_filename) : chainmapping_ttedopa(sdf_filename)
    # Plot
    freqs_x = collect(1:length(coefficients.frequencies))
    coups_x = collect(1:length(coefficients.couplings))
    p = plot(freqs_x, coefficients.frequencies, label = L"\omega_n")
    plot!(coups_x, coefficients.couplings, label = L"\kappa_n")
    plot!(ylabel = L"\textbf{Chain coefficients}", xlabel = L"n")
    # Save plot
    mkpath("./sdf/figs")
    savefig("./sdf/figs/$(sdf_type)_a_" * string(a) * "_T_" * string(T) * "_coeff.json")
end

