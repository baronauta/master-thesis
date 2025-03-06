function makedir_checks(map_tomo_path)
    # Inside the directory with tomograpich states dynamics, create the directory checks
    checks_path = map_tomo_path * "/checks"
    mkpath(checks_path)
    # Read config file
    config_path = map_tomo_path * "/config.JSON"
    config = load_config(config_path) 
    return checks_path, config
end


function check_tomostates_dynamics(map_tomo_path)
    checks_path, config = makedir_checks(map_tomo_path)
    # Evolve the same initial state using
    # - Quantum map (map)
    # - Master equation with generator (me)
    # - Kraus decomposition of the generator (kd)
    # Verify consistency taking as benchmark the time evolution of the state "Up" with TDVP1.
    init_state = Complex{Float64}[1.0 0.0; 0.0 0.0] # Up state
    # TDVP1
    evolvedUp = get_evolved_states(map_tomo_path * "/measurements_Up.dat")
    # QUANTUM MAP
    evolved_map = map_dynamics(map_tomo_path, init_state)
    # MASTER EQUATION
    evolved_me =  me_dynamics(map_tomo_path, init_state, config.dt)
    # KRAUS DECOMPOSITION
    evolved_kd = kd_dynamics(map_tomo_path, init_state, config.dt)

    p = plot(legend=:topright, layout=(3,1), size=(800, 1000))
    # Helper function to plot comparison with TDVP1!
    function compare!(subplot, tdvp_data, data, title)
        ts = config.dt * collect(1:length(data))
        plot!(subplot, title=title, xlabel=L"t")
        plot!(subplot, ts, real(map(x -> tr(σ1 * x), data)), label="σ_x")
        plot!(subplot, ts, real(map(x -> tr(σ2 * x), data)), label="σ_y")
        plot!(subplot, ts, real(map(x -> tr(σ3 * x), data)), label="σ_z")
        plot!(subplot, ts, real(map(x -> tr(σ1 * x), tdvp_data)), label="σ_x (TDVP1)", linestyle=:dash, color=:black)
        plot!(subplot, ts, real(map(x -> tr(σ2 * x), tdvp_data)), label="σ_y (TDVP1)", linestyle=:dash, color=:black)
        plot!(subplot, ts, real(map(x -> tr(σ3 * x), tdvp_data)), label="σ_z (TDVP1)", linestyle=:dash, color=:black)
    end
    # Plot
    compare!(p[1], evolvedUp, evolved_map, "Quantum map vs TDVP1")
    compare!(p[2], evolvedUp, evolved_me, "Master equation vs TDVP1")
    compare!(p[3], evolvedUp, evolved_kd, "Kraus decomposition vs TDVP1")
    # Save image
    savefig(checks_path * "/check_tomostates_dynamics.png")
end


function check_Ks(map_tomo_path, row_idx, col_idx; t_max=Inf)
    figs_path, config = makedir_checks(map_tomo_path)
    p = plot(legend=:topright, layout=(2,1))

    Ks = effective_hamiltonian(map_tomo_path)
    Ks_check = effective_hamiltonian_check(map_tomo_path)
    ReKs = [real(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ImKs = [imag(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ReKs_check = [real(Ks_check[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ImKs_check = [imag(Ks_check[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ts = config.dt * collect(1:length(Ks))
    ts_filtered, ReKs_filtered = time_filter(ts, ReKs; t_max=t_max)
    _, ImKs_filtered = time_filter(ts, ImKs; t_max=t_max)
    _, ReKs_check_filtered = time_filter(ts, ReKs_check; t_max=t_max)
    _, ImKs_check_filtered = time_filter(ts, ImKs_check; t_max=t_max)

    plot!(p[1], ts_filtered, ReKs_filtered, xlabel=L"t", ylabel=L"\textbf{Re}(K_{%$(row_idx)%$(col_idx)})", 
        xlabelfontsize=16, ylabelfontsize=16, grid=false, legend=false
    )
    plot!(p[1], ts_filtered, ReKs_check_filtered, linestyle=:dash, color=:black)
    plot!(p[2], ts_filtered, ImKs_filtered, xlabel=L"t", ylabel=L"\textbf{Im}(K_{%$(row_idx)%$(col_idx)})", 
        xlabelfontsize=16, ylabelfontsize=16, grid=false, legend=false
    )
    plot!(p[2], ts_filtered, ImKs_check_filtered, linestyle=:dash, color=:black)
    
    plot!(p, plot_title=L"\epsilon=%$(config.ϵ),\,\Delta=%$(config.Δ),\,\alpha=%$(config.α),\,T=%$(config.T)")

    savefig(figs_path * "/Ks" * string(row_idx) * string(col_idx) * ".png")
end