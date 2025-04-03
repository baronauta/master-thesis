"Create and fill a file with all the parameters that defines a simulation."
function save_config(
    config_file,
    α,
    ωc,
    sdf_eq,
    domain,
    T,
    chain_length,
    ϵ,
    Δ,
    dt,
    tmax,
    growMPSval,
    local_dim,
)
    config = Dict(
        "α" => α,
        "ωc" => ωc,
        "T" => T,
        "ϵ" => ϵ,
        "Δ" => Δ,
        "dt" => dt,
        "tmax" => tmax,
        "growMPSval" => growMPSval,
        "chain_length" => chain_length,
        "local_dim" => local_dim,
        "sdf_eq" => sdf_eq,
        "domain" => domain
    )
    # Write to JSON file
    open(config_file, "w") do file
        JSON.print(file, config, 4)
    end
end

"Read the configuration file and returns parameters in a dctionary."
function load_config(filename)
    # Read and parse the JSON file
    config = JSON.parsefile(filename)
    # Extract values and return them
    return (
        α = config["α"],
        ωc = config["ωc"],
        T = config["T"],
        ϵ = config["ϵ"],
        Δ = config["Δ"],
        dt = config["dt"],
        tmax = config["tmax"],
        growMPSval = config["growMPSval"],
        chain_length = config["chain_length"],
        local_dim = config["local_dim"],
        sdf_eq = config["sdf_eq"],
        domain = config["domain"]
    )
end

""
function save_tedopa_coefficients(coefficients, save_to_dir, sdf_name)
    freqs_file = save_to_dir * "/$(sdf_name)_freqs.csv"
    coups_file = save_to_dir * "/$(sdf_name)_coups.csv"
    # Frequencies
    io = open(freqs_file, "w")
    for x in coefficients.frequencies
        writedlm(io, x, ',')
    end
    close(io)
    # Couplings
    io = open(coups_file, "w")
    for x in coefficients.couplings
        writedlm(io, x, ',')
    end
    close(io)
    return freqs_file, coups_file
end

"Execute map tomography."
function map_tomography(
    sdf_filename,
    ϵ,
    Δ,
    dt,
    tau;
    growMPSval = 10,
    local_dim = 6,
    tomoStates = ["Up", "Dn", "+", "i"]
)
    # Create directory to store data
    dir = "./runs/" * Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkpath(dir)
    # Create a configuration file with all the parameters
    config_file = dir * "/config.json"
    # - from sdf.json
    input = open(sdf_filename)
    s = read(input, String)
    p = JSON.parse(s)
    α, ωc = Float64.(p["environment"]["spectral_density_parameters"])
    sdf_eq = String(p["environment"]["spectral_density_function"])
    domain = Float64.(p["environment"]["domain"])
    T = Float64(p["environment"]["temperature"])
    chain_length = Int64(p["chain_length"])
    tmax=Int(round(tau*π/Δ))
    save_config(
        config_file,
        α,
        ωc,
        sdf_eq,
        domain,
        T,
        chain_length,
        ϵ,
        Δ,
        dt,
        tmax,
        growMPSval,
        local_dim,
    )
    # Compute TEDOPA or T-TEDOPA coefficients according to temperature T
    coefficients =
        T == 0 ? TEDOPA.chainmapping_tedopa(sdf_filename) : TEDOPA.chainmapping_ttedopa(sdf_filename)
    # Save coefficients into files
    coeff_dir = dir * "/chain_coefficients"
    mkpath(coeff_dir)
    # Extract base filename without directory and extension
    sdf_name = splitext(basename(sdf_filename))[1]
    freqs_file, coups_file = save_tedopa_coefficients(coefficients, coeff_dir, sdf_name)
    # Evolve tomographic states
    for case in tomoStates
        # Define overall system initial state and overall hamiltonian
        (sysenv, psi0) = defineSystem(
            sys_type = "S=1/2",
            sys_istate = case,
            chain_length = chain_length,
            local_dim = local_dim,
        )
        H = createMPO(sysenv, ϵ, Δ, "Sz", freqs_file, coups_file)
        # Measurements
        # - Sx, Sy, Sz on reduced system
        operators = [
            LocalOperator(Dict(1 => "Sx")),
            LocalOperator(Dict(1 => "iSy")),
            LocalOperator(Dict(1 => "Sz")),
        ]
        # Measurument time step is equal to integration time step
        dt_meas = dt
        cb = ExpValueCallback(operators, sysenv, dt_meas)
        # Grow initial state MPS for tdvp1!
        psi0ext = growMPS(psi0, growMPSval)[1]
        # Run TDVP1
        tmpfile = tempname()
        tdvp1!(
            psi0ext,
            H,
            dt,
            tmax;
            hermitian = true,
            normalize = false,
            callback = cb,
            progress = true,
            store_psi0 = false,
            io_file = tmpfile,
            io_ranks = "NUL",
            io_times = "NUL",
        )
        # Store results
        cp(tmpfile, dir * "/measurements_" * case * ".dat", force = true)
    end
end