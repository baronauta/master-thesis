# ─────────────────────────────────────────────────────────────
# Internal functions
# ─────────────────────────────────────────────────────────────

"Scale an existing time array `ts` by Δ / π."
function _scalets(ts, Delta)
    scale = Delta / π
    return scale .* ts
end

"Shorten time domain for plots."
function _timefilter(ts, ys, tmax)
    mask = ts .<= tmax
    return ts[mask], ys[mask]
end

# ─────────────────────────────────────────────────────────────
# Density matrix
# - dynamics of the tomographic states
# - trace distance between tomographic states
# ─────────────────────────────────────────────────────────────

function plot_tomostates(dirdata)
    measurement_files = [
        ("/measurements_Up.dat", "Up"),
        ("/measurements_Dn.dat", "Down"),
        ("/measurements_+.dat", "Plus"),
        ("/measurements_i.dat", "Trans"),
    ]

    config = loadconfig(dirdata * "/config.json")

    f = Figure()

    for (i, (file, title)) in enumerate(measurement_files)
        # Create a new axis for each subplot (4 rows, 1 column)
        ax = Axis(
            f[i, 4],
            xlabel = L"t \Delta / \pi",
            ylabel = "$(title)",
            width = 800,
            height = 200,
        )
        ylims!(ax, -1, 1)
        # get_measurements returns a named tuple ( time = meas["time"], result = meas["..."] )
        meas = get_measurements(dirdata * file, "densitymatrix")
        rho = meas.result
        ts_scaled = _scalets(meas.time, config.Delta)
        lines!(
            ax,
            ts_scaled,
            real.(map(x -> tr(σ1 * x), rho)),
            label = L"\langle\sigma_x\rangle",
        )
        lines!(
            ax,
            ts_scaled,
            real.(map(x -> tr(σ2 * x), rho)),
            label = L"\langle\sigma_y\rangle",
        )
        lines!(
            ax,
            ts_scaled,
            real.(map(x -> tr(σ3 * x), rho)),
            label = L"\langle\sigma_z\rangle",
        )
        # plot the norm of the state
        measnorm = get_measurements(dirdata * file, "norm")
        ts_scaled = _scalets(measnorm.time, config.Delta)
        lines!(ax, ts_scaled, measnorm.result, label = L"Tr(\rho)")
        axislegend(position = :rt)
    end

    resize_to_layout!(f)
    f
end

function _trace_distance(rho1, rho2)
    # Compute the difference matrix (Hermitian)
    delta = rho1 - rho2
    # Compute eigenvalues (real since delta is Hermitian)
    lambda = eigvals(delta)
    # Compute trace distance using absolute eigenvalues
    return 0.5 * sum(abs.(lambda))
end

function plot_tomo_trdistance(dirdata)

    measUp = get_measurements(dirdata * "/measurements_Up.dat", "densitymatrix")
    measDown = get_measurements(dirdata * "/measurements_Dn.dat", "densitymatrix")
    measPlus = get_measurements(dirdata * "/measurements_+.dat", "densitymatrix")
    measTrans = get_measurements(dirdata * "/measurements_i.dat", "densitymatrix")

    # xs: time scaled by Δ / π
    config = loadconfig(dirdata * "/config.json")
    ts_scaled = _scalets(measUp.time, config.Delta)

    # ys: density Matrix
    Up = measUp.result
    Down = measDown.result
    Trans = measTrans.result
    Plus = measPlus.result

    f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"t \Delta / \pi",
        ylabel = L"\text{d}_\text{Tr}\left(\rho_1,\rho_2\right)",
    )
    ylims!(ax, 0, 1)

    # Plot the trace distances
    lines!(ax, ts_scaled, _trace_distance.(Up, Down), label = L"\text{d(Up,Down)}")
    lines!(ax, ts_scaled, _trace_distance.(Up, Plus), label = L"\text{d(Up,Plus)}")
    lines!(ax, ts_scaled, _trace_distance.(Up, Trans), label = L"\text{d(Up,Trans)}")
    lines!(ax, ts_scaled, _trace_distance.(Down, Plus), label = L"\text{d(Down,Plus)}")
    lines!(ax, ts_scaled, _trace_distance.(Down, Trans), label = L"\text{d(Down,Trans)}")
    lines!(ax, ts_scaled, _trace_distance.(Plus, Trans), label = L"\text{d(Plus,Trans)}")

    axislegend(position = :rt)

    f
end

# ─────────────────────────────────────────────────────────────
# Thermodynamics
# - Effective Hamiltonian Ks
# ─────────────────────────────────────────────────────────────

function plot_Ks(dirdata::String, row_idx, col_idx; tmax = nothing)
    effective_hamiltonian = computeKs(dirdata::String)
    # xs: time scaled by Δ / π
    config = loadconfig(dirdata * "/config.json")
    ts = _scalets(effective_hamiltonian.time, config.Delta)
    # ys: Real and Imag part of Effective Hamiltonian Ks
    Ks = effective_hamiltonian.Ks
    ReKs = [real(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]
    ImKs = [imag(Ks[t][row_idx+1, col_idx+1]) for t in eachindex(Ks)]

    if tmax !== nothing
        # Shorten time domain
        mask = ts .<= tmax
        ts = ts[mask]
        ReKs = ReKs[mask]
        ImKs = ImKs[mask]
    end

    f = Figure()
    # Real part of Ks
    ax1 = Axis(
        f[1, 1],
        xlabel = L"t \Delta / \pi",
        ylabel = L"\text{Re}(K_{%$(row_idx)%$(col_idx)})",
    )
    lines!(ax1, ts, ReKs, label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])")
    axislegend(position = :rt)
    # Imaginary part of Ks
    ax2 = Axis(
        f[2, 1],
        xlabel = L"t \Delta / \pi",
        ylabel = L"\text{Im}(K_{%$(row_idx)%$(col_idx)})",
    )
    lines!(ax2, ts, ImKs, label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])")
    axislegend(position = :rt)

    f
end
