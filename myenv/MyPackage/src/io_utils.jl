function save_config(
    config_file,
    α,
    ωc,
    T,
    ϵ,
    Δ,
    dt,
    tmax,
    growMPSval,
    chain_size,
    local_dim,
    sdf_eq,
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
        "chain_size" => chain_size,
        "local_dim" => local_dim,
        "sdf_eq" => sdf_eq,
    )
    # Write to JSON file
    open(config_file, "w") do file
        JSON.print(file, config, 4)
    end
end


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
        chain_size = config["chain_size"],
        local_dim = config["local_dim"],
        sdf_eq = config["sdf_eq"],
    )
end


function load_sdf(sdf_filename)
    # load sdf parameters
    input = open(sdf_filename)
    s = read(input, String)
    sdf_dict = JSON.parse(s)
    return sdf_dict
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

function sdf_naming(sdf_eq::String)
    if sdf_eq == "pi/2* a[1] * x * exp(-x/a[2])"
        return "ohmic"
    elseif sdf_eq == "pi/2 * a[1] * x * a[2]^2 / ( x^2 + a[2]^2 )"
        return "debye"
    else
        return "unknown"
    end
end


function makedir_simul(sdf_filename, ϵ, Δ, dt, tmax, growMPSval, local_dim)
    # Load spectral density parameters and simulation settings from SDF file
    sdf_dict = load_sdf(sdf_filename)
    α, ωc = Float64.(sdf_dict["environment"]["spectral_density_parameters"])
    T = Float64(sdf_dict["environment"]["temperature"])
    sdf_eq = String(sdf_dict["environment"]["spectral_density_function"])
    chain_size = Int64(sdf_dict["chain_length"])
    # Create directory to store data of simulation
    # ./runs/sdf_type/eps_x.x_Delta_x.x/a_x.x_T_x.x/timestamp
    sdf_type = sdf_naming(sdf_eq)
    dir_name =
        "./runs/" *
        sdf_type *
        "/eps_$(ϵ)_Delta_$(Δ)" *
        "/a_$(α)_T_$(T)/" *
        Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkpath(dir_name)
    println("Simulation results will be stored in $(dir_name)")
    # Create and save the configuration file inside dir_name
    config_file = joinpath(dir_name, "config.json")
    save_config(
        config_file,
        α,
        ωc,
        T,
        ϵ,
        Δ,
        dt,
        tmax,
        growMPSval,
        chain_size,
        local_dim,
        sdf_eq,
    )
    return dir_name, config_file
end
