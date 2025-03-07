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


function makedir_simul(sdf_filename, ϵ, Δ, dt, tmax, growMPSval, local_dim)
    # Load spectral density function parameters from SDF file
    sdf_params(sdf_filename)
    α, ωc, T, sdf_eq, chain_size = Float64.(sdf_dict["environment"]["spectral_density_parameters"])
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
