# ─────────────────────────────────────────────────────────────
# Given a sdf J(ω) compute the chain coefficients.
# ─────────────────────────────────────────────────────────────
"""
    chain_coefficients(parameters::Dict{<:AbstractString,Any})

Compute the chain coefficients based on the system's temperature.
This is a wrapper function for `chainmapping_tedopa` and `chainmapping_ttedopa` from the package `TEDOPA`.
Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by the `parameters` dictionary.
"""
function chain_coefficients(parameters::Dict{<:AbstractString,Any})
    T = Float64(parameters["environment"]["temperature"])
    # TEDOPA or T-TEDOPA according to temperature T
    coefficients =
        T == 0 ? TEDOPA.chainmapping_tedopa(parameters) :
        TEDOPA.chainmapping_ttedopa(parameters)
    return coefficients
end

"""
    chain_coefficients(filename::String)

Compute the chain coefficients based on the system's temperature.
This is a wrapper function for `chainmapping_tedopa` and `chainmapping_ttedopa` from the package `TEDOPA`.
Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file called `filename`.
"""
function chain_coefficients(filename::String)
    p = open(filename) do input
        JSON.parse(read(input, String))
    end
    chain_coefficients(p)
end

"""
    save_coefficients_to_csv(coefficients, save_to_dir)

Save the `coefficients` of a TEDOPA chain mapping, which contains `frequencies` and `couplings` fields, 
to CSV files. A folder named `chain_coefficients` will be created in the specified directory `save_to_dir`
to store the CSV files.
"""
function save_tedopa_coefficients(coefficients, save_to_dir)
    mkpath(save_to_dir * "/chain_coefficients")
    freqs_file = save_to_dir * "/chain_coefficients/freqs.csv"
    coups_file = save_to_dir * "/chain_coefficients/coups.csv"
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

# ─────────────────────────────────────────────────────────────
# Simulation
# ─────────────────────────────────────────────────────────────

"""
    tomodynamics(filename::String; tomoStates = ["Up", "Dn", "+", "i"])

# Arguments
- `filename::String`: Path to the JSON configuration file.
- `tomoStates::Vector{String}`: (Optional) List of initial states to evolve. Defaults to `["Up", "Dn", "+", "i"]`.
"""
function tomodynamics(filename::String; tomoStates = ["Up", "Dn", "+", "i"])
    config = loadconfig(filename)
    # Compute TEDOPA or T-TEDOPA coefficients according to temperature T
    coefficients = chain_coefficients(filename)
    # Store chain coefficients in CSV files;
    # these will be used in MPO definition.
    dir = dirname(filename)
    freqs_file, coups_file = save_tedopa_coefficients(coefficients, dir)

    # Evolve tomographic states
    for case in tomoStates
        # MPS of the initial state of the composite system TLS⊗bath
        (sysenv, psi0) = defineSystem(
            sys_type = "S=1/2",
            sys_istate = case,
            chain_length = config.chain_length,
            local_dim = config.local_dim,
        )
        # MPO of the Hamiltonian H = H_S + H_E + H_I
        H = createMPO(sysenv, config.epsilon, config.Delta, "Sz", freqs_file, coups_file)

        # Measure Sx, Sy, Sz on the TLS, i.e. determination of the density matrix ρ of the TLS.
        # Measuring "iSy" to workaround an issue with the ITensors (or ITensorMPS?) package.
        operators = [
            LocalOperator(Dict(1 => "Sx")),
            LocalOperator(Dict(1 => "iSy")),
            LocalOperator(Dict(1 => "Sz")),
        ]
        # Measurements performed at every integration step during the simulation
        dtmeas = config.dt
        cb = ExpValueCallback(operators, sysenv, dtmeas)

        # Increase manually the bond dimension of the initial state (which is factorized, so that χ=1 on all bonds)
        psi0ext = growMPS(psi0, config.grow_mps)[1]

        # Calculate tmax from tΔ/π: tmax = tΔ/π * π/Δ
        tmax = config.t_Delta_over_pi * π / config.Delta

        # Run TDVP1
        tmpfile = tempname()
        tdvp1!(
            psi0ext,
            H,
            config.dt,
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

"""
    envdynamics(filename::String)

# Arguments
- `filename::String`: Path to the JSON configuration file.
"""
function envdynamics(filename)
    config = loadconfig(filename)
    # Compute TEDOPA or T-TEDOPA coefficients according to temperature T
    coefficients = chain_coefficients(filename)
    # Store chain coefficients in CSV files;
    # these will be used in MPO definition.
    dir = dirname(filename)
    freqs_file, coups_file = save_tedopa_coefficients(coefficients, dir)

    # MPS of the initial state of the composite system TLS⊗bath;
    # consider the TLS to be in the "Up" state.
    (sysenv, psi0) = defineSystem(
        sys_type = "S=1/2",
        sys_istate = "Up",
        chain_length = config.chain_length,
        local_dim = config.local_dim,
    )
    # MPO of the Hamiltonian H = H_S + H_E + H_I
    H = createMPO(sysenv, config.epsilon, config.Delta, "Sz", freqs_file, coups_file)

    # Measure the ProjUp observable on the first site of the MPS, i.e. the one for the TLS.
    # Measure N on each site of the chain: these are the ones from 2 to chain_length + 1.
    operators = [LocalOperator(Dict(1 => "ProjUp"))]
    for i = 2:config.chain_length+1
        push!(operators, LocalOperator(Dict(i => "N")))
    end
    # Measurements are performed every n integration steps to reduce computational cost. 
    # This lower resolution is sufficient for the plots while speeding up the simulation.
    dtmeas = 10 * config.dt
    cb = ExpValueCallback(operators, sysenv, dtmeas)

    # Increase manually the bond dimension of the initial state (which is factorized, so that χ=1 on all bonds)
    psi0ext = growMPS(psi0, config.grow_mps)[1]

    # Calculate tmax from tΔ/π: tmax = tΔ/π * π/Δ
    tmax = config.t_Delta_over_pi * π / config.Delta

    # Run TDVP1
    tmpfile = tempname()
    tdvp1!(
        psi0ext,
        H,
        config.dt,
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
    cp(tmpfile, dir * "/measurements_N.dat", force = true)
end

"""
    get_measurements(filename::String, measuretype::String)
"""
function get_measurements(filename::String, measuretype::String)
    # DelimitedFiles.readdlm: 
    # Read a matrix from the source where each line (separated by eol) gives one row, with elements separated by the given delimiter.
    (rawdata, header) = readdlm(filename, ',', Float64, header = true)
    meas = Dict(header[i] => rawdata[:, i] for i in eachindex(header))

    # This function will return a named tuple (time = meas["time"], result = meas["..."] )
    ts = meas["time"]

    if measuretype == "densitymatrix"
        sx = 2 .* meas["Sx{1}_re"]
        sy = 2 .* meas["iSy{1}_im"]
        sz = 2 .* meas["Sz{1}_re"]
        ρt = [0.5 * (σ0 .+ sx[i] * σ1 .+ sy[i] * σ2 .+ sz[i] * σ3) for i in eachindex(sx)]
        return (time = ts, result = ρt)

    elseif measuretype == "norm"
        return (time = ts, result = meas["Norm_re"])

    elseif measuretype == "N"
        cols = filter(h -> occursin(r"^N\{\d+\}_re$", h), header)
        return (time = ts, result = Dict(h => meas[h] for h in cols))

    else
        error(
            "Unsupported measure type: $measuretype. Choose from \"densitymatrix\", \"norm\", or \"N\".",
        )
    end
end
