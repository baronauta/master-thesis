"Prepare the directory to store figures. Return path of this directory and return the dictionary with all the parameters of the simulation."
function makedir_figs(dir_data)
    # Create a subfolder figs inside the directory with data
    figs_path = dir_data * "/figs"
    mkpath(figs_path)
    # Read config file
    config_path = dir_data * "/config.json"
    config = load_config(config_path)
    return figs_path, config
end

"Plot the spectral density function."
function fig_sdf(dir_data)
    figs_path, config = makedir_figs(dir_data)
    # Extracting parameters
    a = config.α
    T = config.T
    support = config.domain
    # Retrieves a string representing the spectral density function formula
    fn = config.sdf_eq
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
    display(p)
    savefig(figs_path * "/sdf.png")
end

"Plot tedopa coefficients. Read from sdf_filename. Display the graph and returns interesting values for setting simulation."
### Rifare leggendo da directory dove dovrebbero esserci già file da leggere con i coefficienti
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

"Plots time evolved states of the tomographic basis."
function fig_tomo(dir_data)
    figs_path, config = makedir_figs(dir_data)
    measurement_files = [
        ("/measurements_Up.dat", "Up"),
        ("/measurements_Dn.dat", "Down"),
        ("/measurements_+.dat", "Plus"),
        ("/measurements_i.dat", "Trans"),
    ]
    p = plot(legend = :topright, layout = (4, 1), size = (600, 800))
    for (i, (file, title)) in enumerate(measurement_files)
        rho = get_evolved_states(dir_data * file)
        norm = get_norm(dir_data * file)
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

"Compute trace-distance between density matrix. Take advantage of the fact that they are Hemritian."
function trace_distance(rho1, rho2)
    # Compute the difference matrix (Hermitian)
    delta = rho1 - rho2
    # Compute eigenvalues (real since delta is Hermitian)
    lambda = eigvals(delta)
    # Compute trace distance using absolute eigenvalues
    return 0.5 * sum(abs.(lambda))
end

"Plot trace-distance between density matrix from tomographic basis."
function fig_trdistance(dir_data)
    figs_path, config = makedir_figs(dir_data)
    # Read meaurements data from file
    evolvedUp = get_evolved_states(dir_data * "/measurements_Up.dat")
    evolvedDown = get_evolved_states(dir_data * "/measurements_Dn.dat")
    evolvedPlus = get_evolved_states(dir_data * "/measurements_+.dat")
    evolvedTrans = get_evolved_states(dir_data * "/measurements_i.dat")
    p = plot(legend = :topright)
    ts = config.dt * (1:length(evolvedUp))
    plot!(p, title = L"\epsilon=%$(config.ϵ),\,\Delta=%$(config.Δ)\,-\,\alpha=%$(config.α),\,T=%$(config.T)", xlabel = L"t", ylabel=L"\textbf{d}_\textbf{Tr}\left(\rho_1,\rho_2\right)")
    plot!(p, ts, trace_distance.(evolvedUp, evolvedDown), label = L"\textbf{Up,Down}")
    plot!(p, ts, trace_distance.(evolvedUp, evolvedPlus), label = L"\textbf{Up,Plus}")
    plot!(p, ts, trace_distance.(evolvedUp, evolvedTrans), label = L"\textbf{Up,Trans}")
    plot!(p, ts, trace_distance.(evolvedDown, evolvedPlus), label = L"\textbf{Down,Plus}")
    plot!(p, ts, trace_distance.(evolvedDown, evolvedTrans), label = L"\textbf{Down,Trans}")
    plot!(p, ts, trace_distance.(evolvedPlus, evolvedTrans), label = L"\textbf{Plus,Trans}")
    plot!(p, ylims=(0,1))
    display(p)
    savefig(figs_path * "/tr_distance.png")
end

"Shorten time domain for plots."
function time_filter(ts, ys; t_max = Inf)
    ts_filtered = ts[ts.<=t_max]
    ys_filtered = ys[ts.<=t_max]
    return ts_filtered, ys_filtered
end

"Plot Effective Hamiltonian Ks."
function fig_Ks(dir_data, row_idx, col_idx; t_max = Inf)
    figs_path, config = makedir_figs(dir_data)
    p = plot(legend = :topright, layout = (2, 1))
    Ks = effective_hamiltonian(dir_data)
    ReKs = [real(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ImKs = [imag(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ts = config.dt * collect(1:length(Ks))
    # Shorten time domain if t_max is not Inf
    ts_filtered, ReKs_filtered = time_filter(ts, ReKs; t_max = t_max)
    _, ImKs_filtered = time_filter(ts, ImKs; t_max = t_max)
    # Plot Real part
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
    # Plot Imag part
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

"Plot internal energy."
function fig_Us(dir_data; t_max = Inf)
    figs_path, config = makedir_figs(dir_data)
    measurement_files = [
        ("/measurements_Up.dat", "Up"),
        ("/measurements_Dn.dat", "Down"),
        ("/measurements_+.dat", "Plus"),
        ("/measurements_i.dat", "Trans"),
    ]
    p = plot(legend = :topright, layout = (4, 1), size = (600, 700))
    Ks = effective_hamiltonian(dir_data)
    ts = config.dt * collect(1:length(Ks))
    for (i, (file, title)) in enumerate(measurement_files)
        rho = get_evolved_states(dir_data * file)
        Us = internal_energy(Ks, rho)
        ts_filtered, Us_filtered = time_filter(ts, Us; t_max = t_max)
        plot!(p[i], title = title, xlabel = L"t", ylabel = L"Tr(\rho K_s)")
        plot!(p[i], ts_filtered, Us_filtered, legend = false)
    end
    display(p)
    savefig(figs_path * "/internal_energy.png")
end
