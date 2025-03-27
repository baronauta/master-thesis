"Return all the parameters about environment and chain length from sdf JSON file."
function sdf_params(sdf_filename)
    # Open the JSON file and parse the dictionary
    input = open(sdf_filename)
    s = read(input, String)
    sdf_dict = JSON.parse(s)
    # Read all parameters
    # Environment
    α, ωc = Float64.(sdf_dict["environment"]["spectral_density_parameters"])
    sdf_eq = String(sdf_dict["environment"]["spectral_density_function"])
    domain = Float64.(sdf.dict["environment"]["domain"])
    T = Float64(sdf_dict["environment"]["temperature"])
    # Chain length
    chain_length = Int64(sdf_dict["chain_length"])
    return α, ωc, sdf_eq, domain, T, chain_length
end

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

"Prepare the directory to store data from the simulation. Desidered format is ./runs/sdf_type/eps_x.x_Delta_x.x/a_x.x_T_x.x/timestamp."
function makedir_simul(sdf_filename, ϵ, Δ, dt, tmax, growMPSval, local_dim)
    # Load spectral density function parameters from SDF file
    α, ωc, sdf_eq, domain, T, chain_length = sdf_params(sdf_filename)
    # Choose the name associated to the sdf
    if sdf_eq == "pi/2* a[1] * x * exp(-x/a[2])"
        sdf_type = "ohmic"
    elseif sdf_eq == "pi/2 * a[1] * x * a[2]^2 / ( x^2 + a[2]^2 )"
        sdf_type = "debye"
    else
        sdf_type = "unknown"
    end
    # Directory name
    dir_name =
        "./runs/" *
        sdf_type *
        "/eps_$(ϵ)_Delta_$(Δ)" *
        "/a_$(α)_T_$(T)/" *
        Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkpath(dir_name)
    println("Simulation results will be stored in $(dir_name)")
    # Create and save the configuration file inside dir_name
    config_file = dir_name * "\config.json"
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
    return dir_name
end
