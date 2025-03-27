function figs_sdf(sdf_filename)
    # open and read the JSON file
    input = open(sdf_filename)
    s = read(input, String)
    # parsing JSON data: parse the string s as a JSON object and assigns the resulting dictionary-like structure to p
    p = JSON.parse(s)
    # extracting parameters
    a = Float64.(p["environment"]["spectral_density_parameters"])
    T = Float64(p["environment"]["temperature"])
    support = Float64.(p["environment"]["domain"])
    # retrieves a string representing the spectral density function formula
    fn = p["environment"]["spectral_density_function"]
    sdf_type = sdf_naming(fn)
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
    savefig("./sdf/figs/$(sdf_type)_a_" * string(a[1]) * "_T_" * string(T) * ".png")
end


function figs_tedopa_coefficients(sdf_filename)
    α, ωc, T, sdf_eq, chain_size = sdf_params(sdf_filename)
    sdf_type = sdf_naming(sdf_eq)
    # TEDOPA or T-TEDOPA according to temperature T
    coefficients =
        T == 0 ? chainmapping_tedopa(sdf_filename) : chainmapping_ttedopa(sdf_filename)
    # Plot
    freqs_x = collect(1:length(coefficients.frequencies))
    coups_x = collect(1:length(coefficients.couplings))
    p = plot(freqs_x, coefficients.frequencies, label = L"\omega_n")
    plot!(coups_x, coefficients.couplings, label = L"\kappa_n")
    plot!(ylabel = L"\textbf{Chain\,\,\,coefficients}", xlabel = L"n")
    # Save plot
    mkpath("./sdf/figs")
    savefig("./sdf/figs/$(sdf_type)_a_" * string(α) * "_T_" * string(T) * "_coeff.png")
end

function makedir_figs(map_tomo_path)
    # Create the directory figs
    figs_path = map_tomo_path * "/figs"
    mkpath(figs_path)
    # Read config file
    config_path = map_tomo_path * "/config.json"
    config = load_config(config_path)
    return figs_path, config
end

function time_filter(ts, ys; t_max = Inf)
    ts_filtered = ts[ts.<=t_max]
    ys_filtered = ys[ts.<=t_max]
    return ts_filtered, ys_filtered
end


function fig_tomo(map_tomo_path)
    figs_path, config = makedir_figs(map_tomo_path)
    p = plot(legend = :topright, layout = (4, 1), size = (600, 800))

    measurement_files = [
        ("/measurements_Up.dat", "Up"),
        ("/measurements_Dn.dat", "Down"),
        ("/measurements_+.dat", "Plus"),
        ("/measurements_i.dat", "Trans"),
    ]

    for (i, (file, title)) in enumerate(measurement_files)
        rho = get_evolved_states(map_tomo_path * file)
        norm = get_norm(map_tomo_path * file)
        ts = config.dt * (1:length(rho))

        plot!(p[i], title = title, xlabel = L"t")
        plot!(p[i], ts, real(map(x -> tr(σ1 * x), rho)), label = L"\langle\sigma_x\rangle")
        plot!(p[i], ts, real(map(x -> tr(σ2 * x), rho)), label = L"\langle\sigma_y\rangle")
        plot!(p[i], ts, real(map(x -> tr(σ3 * x), rho)), label = L"\langle\sigma_z\rangle")
        plot!(p[i], ts, norm, label = L"Tr(\rho)", lc = :gray)
    end
    display(p)
    savefig(figs_path * "/tomostates_dynamics.png")
end


function fig_Ks(map_tomo_path, row_idx, col_idx; t_max = Inf)
    figs_path, config = makedir_figs(map_tomo_path)
    p = plot(legend = :topright, layout = (2, 1))

    Ks = effective_hamiltonian(map_tomo_path)
    ReKs = [real(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ImKs = [imag(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ts = config.dt * collect(1:length(Ks))
    ts_filtered, ReKs_filtered = time_filter(ts, ReKs; t_max = t_max)
    _, ImKs_filtered = time_filter(ts, ImKs; t_max = t_max)

    plot!(
        p[1],
        ts_filtered,
        ReKs_filtered,
        xlabel = L"t",
        ylabel = L"\textbf{Re}(K_{%$(row_idx)%$(col_idx)})",
        xlabelfontsize = 16,
        ylabelfontsize = 16,
        grid = false,
        label = L"\alpha=%$(config.α),\,T=%$(config.T)",
    )
    plot!(
        p[2],
        ts_filtered,
        ImKs_filtered,
        xlabel = L"t",
        ylabel = L"\textbf{Im}(K_{%$(row_idx)%$(col_idx)})",
        xlabelfontsize = 16,
        ylabelfontsize = 16,
        grid = false,
        label = L"\alpha=%$(config.α),\,T=%$(config.T)",
    )

    plot!(p, plot_title = L"\epsilon=%$(config.ϵ),\,\Delta=%$(config.Δ)")
    display(p)
    savefig(figs_path * "/Ks" * string(row_idx) * string(col_idx) * ".png")
end


function fig_Us(map_tomo_path; t_max = Inf)
    figs_path, config = makedir_figs(map_tomo_path)
    p = plot(legend = :topright, layout = (4, 1), size = (600, 700))

    measurement_files = [
        ("/measurements_Up.dat", "Up"),
        ("/measurements_Dn.dat", "Down"),
        ("/measurements_+.dat", "Plus"),
        ("/measurements_i.dat", "Trans"),
    ]

    Ks = effective_hamiltonian(map_tomo_path)
    ts = config.dt * collect(1:length(Ks))

    for (i, (file, title)) in enumerate(measurement_files)
        rho = get_evolved_states(map_tomo_path * file)
        Us = internal_energy(Ks, rho)
        ts_filtered, Us_filtered = time_filter(ts, Us; t_max = t_max)

        plot!(p[i], title = title, xlabel = L"t", ylabel = L"Tr(\rho K_s)")
        plot!(p[i], ts_filtered, Us_filtered, legend = false)
    end
    display(p)
    savefig(figs_path * "/internal_energy_tomostates.png")
end
