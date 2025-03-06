function tedopa_coefficients(sdf_filename)
    sdf_dict = load_sdf(sdf_filename)
    T = Float64(sdf_dict["environment"]["temperature"])
    # TEDOPA or T-TEDOPA according to temperature T
    coefficients =
        T == 0 ? chainmapping_tedopa(sdf_filename) : chainmapping_ttedopa(sdf_filename)
    # save coefficients
    freqs_file, coups_file = save_tedopa_coefficients(coefficients, sdf_filename)
    return freqs_file, coups_file
end


function evolve_tomographic_states(sdf_filename, config_file, dir_name)
    # load TEDOPA coefficients
    freqs_file, coups_file = tedopa_coefficients(sdf_filename)
    # load config
    config = load_config(config_file)
    dt_step = config.dt
    dt_meas = dt_step
    # map tomography
    tomoStates = ["Up", "Dn", "+", "i"]
    for case in tomoStates
        # define overall system initial state and overall hamiltonian
        (sysenv, psi0) = defineSystem(
            sys_type = "S=1/2",
            sys_istate = case,
            chain_size = config.chain_size,
            local_dim = config.local_dim,
        )
        H = createMPO(sysenv, config.ϵ, config.Δ, "Sz", freqs_file, coups_file)
        # measurements
        # - Sx, Sy, Sz on reduced system
        operators = [
            LocalOperator(Dict(1 => "Sx")),
            LocalOperator(Dict(1 => "iSy")),
            LocalOperator(Dict(1 => "Sz")),
        ]
        cb = ExpValueCallback(operators, sysenv, dt_meas)
        # grow initial state MPS for tdvp1!
        psi0ext = growMPS(psi0, config.growMPSval)[1]
        # run tdvp1!
        tmpfile = tempname()
        tdvp1!(
            psi0ext,
            H,
            dt_step,
            config.tmax;
            hermitian = true,
            normalize = false,
            callback = cb,
            progress = true,
            store_psi0 = false,
            io_file = tmpfile,
            io_ranks = "NUL",
            io_times = "NUL",
        )
        # store results
        cp(tmpfile, dir_name * "/measurements_" * case * ".dat", force = true)
    end
end


function get_evolved_states(filename)
    (rawdata, header) = readdlm(filename, ',', Float64, header = true)
    meas = Dict(header[i] => rawdata[:, i] for i in eachindex(header))
    nmeas = length(meas["time"])
    rho_t = [
        0.5 * (
            σ0 +
            2 * meas["Sx{1}_re"][i] * σ1 +
            2 * meas["iSy{1}_im"][i] * σ2 +
            2 * meas["Sz{1}_re"][i] * σ3
        ) for i = 1:nmeas
    ]
    return rho_t
end


function get_norm(filename)
    (rawdata, header) = readdlm(filename, ',', Float64, header = true)
    meas = Dict(header[i] => rawdata[:, i] for i in eachindex(header))
    return meas["Norm_re"]
end


function map_tomography(
    sdf_filename,
    ϵ,
    Δ;
    dt = 0.01,
    tmax = 30,
    growMPSval = 10,
    local_dim = 6,
)
    # Create directory to store data of map tomography, add config file
    dir_name, config_file =
        makedir_simul(sdf_filename, ϵ, Δ, dt, tmax, growMPSval, local_dim)
    # Evolve tomographic states
    evolve_tomographic_states(sdf_filename, config_file, dir_name)
end
