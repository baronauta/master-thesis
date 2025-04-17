"For my checks."
function sdf_constant(; alpha = 0.1, epsilon = 0.5, Delta = 0.2, omegac = 1, s = 1)
    omegas = 2 * sqrt(epsilon^2 + Delta^2)
    π / 2 * alpha / (omegas^s * omegac^(s - 1))
end

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
    println("\n=== TWO-LEVEL SYSTEM SETUP ===")
    epsilon, Delta = _input_tls_parameters()

    println("\n=== SPECTRAL DENSITY FUNCTION SETUP ===")
    sdf, a = _input_spectral_density(epsilon, Delta)

    println("\n=== ENVIRONMENT SETUP ===")
    domain, temp, chain_length = _input_environment()

    println("\n=== SIMULATION SETTINGS ===")
    t_Delta, local_dim, grow_mps = _input_simulation_settings()

    data = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => a,
            "spectral_density_function" => sdf,
            "domain" => domain,
            "temperature" => temp,
        ),
        "chain_length" => chain_length,
        "PolyChaos_nquad" => 5000,
        "TLS" => Dict("epsilon" => epsilon, "Delta" => Delta),
        "simulation" => Dict(
            "dt" => nothing,
            "local_dim" => local_dim,
            "t_Delta_over_pi" => t_Delta,
            "grow_mps" => grow_mps,
        ),
    )

    # Parameter suggestions
    _apply_suggestions!(data)

    _save_config(dirname, data)
end

function _input_tls_parameters()
    println("The system Hamiltonian is: H = ϵσz + Δσx")
    print("Enter ϵ (bias): ")
    epsilon = parse(Float64, readline())
    print("Enter Δ (tunneling amplitude): ")
    Delta = parse(Float64, readline())
    return epsilon, Delta
end

function _input_spectral_density(epsilon, Delta)
    println("You can enter a custom expression, e.g. pi/2 * a[1] * x * exp(-x/a[2]).")
    println(
        "Note: 'a' is a parameter array that you will be asked to provide in the next step.",
    )
    println("Or use a predefined model by name: ohmic.")
    print("Enter spectral density function: ")
    sdf = readline()

    if sdf == "ohmic"
        sdf, a = _handle_predefined_ohmic(epsilon, Delta)
    else
        println("→ Using custom spectral density function:\n  $sdf")
        print(
            "Enter spectral density parameters as comma-separated values (e.g., 0.1,1.0,...): ",
        )
        a = parse.(Float64, split(readline(), ","))
    end
    return sdf, a
end

function _handle_predefined_ohmic(epsilon, Delta)
    sdf = "pi/2 * a[1] * x^(a[3]) / a[2]^(a[3]-1) * exp(-x/a[2])"
    println("→ Using predefined ohmic spectral density:")
    println("  $sdf\n  a[1] (α), a[2] (ωc), a[3] (s)")
    print(
        "Enter the spectral density parameters as comma-separated values (e.g. 0.1,1.0,1.0): ",
    )
    a = parse.(Float64, split(readline(), ","))
    println("Possible rescaling: α → α / (ωs)^s. Dou want to apply? [Y/N]")
    scale = readline()
    if scale == "Y"
        # Rescaling of α with the frequency value of the system.
        # For the considered TLS: ωs = (E+ - E-)/ħ, 
        # where E± eigvals of Hs, is 2√(ϵ^2+Δ^2).
        ωs = 2 * sqrt(epsilon^2 + Delta^2)
        a[1] /= ωs^a[3]
        println("→ Rescaling applied.")
    elseif scale == "N"
        println("→ No rescaling.")
    else
        println("ERROR: invalid input!")
    end

    return sdf, a
end

function _input_environment()
    print("Enter spectral domain maximum as an integer [0,max] (e.g. 10): ")
    domain_max = parse(Int, readline())
    domain = [0, domain_max]

    print("Enter environment temperature: ")
    temp = parse(Float64, readline())

    print("Enter environment chain length: ")
    chain_length = parse(Int, readline())

    return domain, temp, chain_length
end

function _input_simulation_settings()
    print("Enter total simulation time (in units of tΔ/π): ")
    t_Delta = parse(Float64, readline())

    print("Enter local Hilbert space dimension per environment site: ")
    local_dim = parse(Int, readline())

    print("Enter MPS bond dimension growth: ")
    grow_mps = parse(Int, readline())

    return t_Delta, local_dim, grow_mps
end

function _apply_suggestions!(data)
    suggested_dt, suggested_N = _suggest_params(data)

    println("\n=== PARAMETER SUGGESTIONS ===")
    println("Based on your input, we suggest:")
    println("  Chain length  N  = $suggested_N")
    println("  Time step    dt  = $suggested_dt")
    print("Do you want to apply these suggestions? [Y/N]: ")

    change = readline()
    if change == "Y"
        data["simulation"]["dt"] = suggested_dt
        data["chain_length"] = suggested_N
        println("→ Suggestions applied.")
    elseif change == "N"
        println("→ Keeping original chain length N.")
        print("Enter time step dt: ")
        data["simulation"]["dt"] = parse(Float64, readline())
    else
        println("ERROR: invalid input!")
    end
end

function _save_config(dirname::String, data::Dict)
    mkpath("runs/$dirname")
    filename = "runs/$dirname/config.json"
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
