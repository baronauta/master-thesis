"""
    setup(dirname::String)

Interactively prompts the user to configure simulation parameters, including spectral density function,
environment settings, and two-level system parameters. Saves the configuration to a JSON file.

Arguments:
- `dirname::String`: Directory name where the configuration file will be saved.

Prompts:
- Spectral density function (custom or predefined "ohmic").
- Spectral density parameters (`a` array).
- Spectral domain, temperature, and chain length.
- Two-level system parameters (`epsilon`, `Delta`).
- Simulation settings (time step `dt`, Hilbert space dimension, simulation time, and MPS growth).
- Suggests optimal time step (`dt`) and chain length (`N`), with user approval.

The configuration is saved to `runs/<dirname>/config.json`.
"""
function setup(dirname::String)
    println("=== SPECTRAL DENSITY FUNCTION SETUP ===")
    println("You can enter a custom expression, e.g. pi/2 * a[1] * x * exp(-x/a[2]).")
    println(
        "Note: 'a' is a parameter array that you will be asked to provide in the next step.",
    )
    println("Or use a predefined model by name: ohmic, superohmic, subohmic.")
    print("Enter spectral density function: ")
    sdf = readline()

    if sdf == "ohmic"
        sdf = "pi/2* a[1] * x * exp(-x/a[2])"
        println("→ Using predefined Ohmic spectral density:")
        println("  $sdf")
        println("Now enter the spectral density parameters as comma-separated values:")
        println("  a[1] (α): Coupling strength")
        println("  a[2] (ωc): Cut-off frequency")
        print("Enter values (e.g., 0.1,1.0): ")
        a_str = readline()
        a = parse.(Float64, split(a_str, ","))
    elseif sdf == "superohmic"
        sdf = "pi/2 * a[1] * x^(a[3]) / a[2]^(a[3]-1) * exp(-x/a[2])"
        println("→ Using predefined Super-Ohmic spectral density:")
        println("  $sdf")
        println("Now enter the spectral density parameters as comma-separated values:")
        println("  a[1] (α): Coupling strength")
        println("  a[2] (ωc): Cut-off frequency")
        println("  a[3] (s): Spectral type — s>1 (Super-Ohmic)")
        print("Enter values (e.g., 0.1,1.0,1.2): ")
        a_str = readline()
        a = parse.(Float64, split(a_str, ","))
    elseif sdf == "subohmic"
        sdf = "pi/2 * a[1] * x^(a[3]) / a[2]^(a[3]-1) * exp(-x/a[2])"
        println("→ Using predefined Super-Ohmic spectral density:")
        println("  $sdf")
        println("Now enter the spectral density parameters as comma-separated values:")
        println("  a[1] (α): Coupling strength")
        println("  a[2] (ωc): Cut-off frequency")
        println("  a[3] (s): Spectral type — s<1 (Sub-Ohmic)")
        print("Enter values (e.g., 0.1,1.0,0.4): ")
        a_str = readline()
        a = parse.(Float64, split(a_str, ","))
    else
        println("→ Using custom spectral density function:")
        println("  $sdf")
    end

    print("Enter spectral domain maximum as an integer [0,max] (e.g. 10): ")
    domain_str = readline()
    domain_max = parse(Int, domain_str)
    domain = [0, domain_max]

    print("Enter environment temperature: ")
    temp = parse(Float64, readline())

    print("Enter environment chain length: ")
    chain_length = parse(Int, readline())

    # Default value
    PolyChaos_nquad = 5000

    println("\n=== TWO-LEVEL SYSTEM SETUP ===")
    println("The system Hamiltonian is: H = ϵσz + Δσx")
    print("Enter ϵ (bias): ")
    epsilon = parse(Float64, readline())
    print("Enter Δ (tunneling amplitude): ")
    Delta = parse(Float64, readline())

    println("\n=== SIMULATION SETTINGS ===")
    print("Enter time step dt: ")
    dt = parse(Float64, readline())

    print("Enter total simulation time (in units of tΔ/π): ")
    t_Delta = parse(Float64, readline())

    print("Enter local Hilbert space dimension per environment site: ")
    local_dim = parse(Int, readline())

    print("Enter MPS bond dimension growth: ")
    grow_mps = parse(Int, readline())

    data = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => a,
            "spectral_density_function" => sdf,
            "domain" => domain,
            "temperature" => temp,
        ),
        "chain_length" => chain_length,
        "PolyChaos_nquad" => PolyChaos_nquad,
        "TLS" => Dict("epsilon" => epsilon, "Delta" => Delta),
        "simulation" => Dict(
            "dt" => dt,
            "local_dim" => local_dim,
            "t_Delta_over_pi" => t_Delta,
            "grow_mps" => grow_mps,
        ),
    )

    suggested_dt, suggested_N = _suggest_params(data)
    println("\n=== PARAMETER SUGGESTIONS ===")
    println("Based on your input, we suggest:")
    println("  Time step    dt  = $suggested_dt")
    println("  Chain length  N  = $suggested_N")
    print("Do you want to apply these suggestions? [Y/N]: ")
    change = readline()
    if change == "Y"
        data["simulation"]["dt"] = suggested_dt
        data["chain_length"] = suggested_N
        println("→ Suggestions applied.")
    else
        println("→ Keeping original parameters.")
    end

    mkdir("runs")
    newdir = "runs/" * dirname
    mkdir(newdir)
    filename = newdir * "/config.json"
    open(filename, "w") do f
        write(f, JSON.json(data, 4))
    end

    println("\nConfiguration saved to '$filename'")
end

"""
    _suggest_params(parameters::Dict{<:AbstractString, Any})

Calculates suggested time step (`dt`) and chain length (`N`) based on input parameters.
"""
function _suggest_params(parameters::Dict{<:AbstractString,Any})
    # Extract necessary coefficients for parameter calculations
    coefficients = chain_coefficients(parameters)
    k∞ = coefficients.couplings[end] # Asymptotic coupling coefficient

    # Rule-of-thumb formulas to determine the optimal chain_length and time step dt:
    # - dt = 1 / (k∞ * 50)
    # - N = 2 * k∞ * tmax

    # Suggested time step (dt)
    dt = round(1 / (k∞ * 50), sigdigits = 1, base = 10)

    # Maximum simulation time in seconds unit
    Delta = parameters["TLS"]["Delta"]
    t_Delta_over_pi = parameters["simulation"]["t_Delta_over_pi"]
    tmax = t_Delta_over_pi * π / Delta
    # Suggested chain length N
    N = round(Int, 2 * k∞ * tmax)

    return dt, N
end

"""
    loadconfig(filename::String) -> NamedTuple

Load system and simulation parameters from a JSON configuration file.
Return a flat named tuple containing all relevant parameters.
"""
function loadconfig(filename::String)
    config = JSON.parsefile(filename)

    return (
        epsilon = config["TLS"]["epsilon"],
        Delta = config["TLS"]["Delta"],
        dt = config["simulation"]["dt"],
        t_Delta_over_pi = config["simulation"]["t_Delta_over_pi"],
        grow_mps = config["simulation"]["grow_mps"],
        local_dim = config["simulation"]["local_dim"],
        domain = config["environment"]["domain"],
        a = config["environment"]["spectral_density_parameters"],
        sdf = config["environment"]["spectral_density_function"],
        temperature = config["environment"]["temperature"],
        chain_length = config["chain_length"],
    )
end
