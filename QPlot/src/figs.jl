"Shorten time domain for plots."
function _timefilter(ts, ys, tmax)
    mask = ts .<= tmax
    return ts[mask], ys[mask]
end

"""
    plot_state(dirdata; state="Up")

Plot the expectation values of σ_x, σ_y, σ_z and the trace norm for a single state.

# Arguments
- `dirdata::AbstractString`: Directory containing measurement files and config.json.
- `state::AbstractString`: One of "Up", "Down", "Plus", or "Down".
"""
function plot_state(dirdata; state = "Up")
    file = Dict(
        "Up" => "/measurements_Up.dat",
        "Down" => "/measurements_Dn.dat",
        "Plus" => "/measurements_+.dat",
        "Trans" => "/measurements_i.dat",
    )
    meas = get_measurements(dirdata * file[state], "densitymatrix")
    measnorm = get_measurements(dirdata * file[state], "norm")
    # xs: time 
    ts = meas.time
    # ys: density Matrix
    rho = meas.result

    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega_c t", ylabel = "$state")
    ylims!(ax, -1, 1)
    # Plot
    lines!(ax, ts, real.(map(x -> tr(σ1 * x), rho)), label = L"\langle\sigma_x\rangle")
    lines!(ax, ts, real.(map(x -> tr(σ2 * x), rho)), label = L"\langle\sigma_y\rangle")
    lines!(ax, ts, real.(map(x -> tr(σ3 * x), rho)), label = L"\langle\sigma_z\rangle")
    lines!(ax, ts, measnorm.result, label = L"Tr(\rho)")
    axislegend(position = :rt)

    return f
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

    # xs: time 
    ts = measUp.time

    # ys: density Matrix
    Up = measUp.result
    Down = measDown.result
    Trans = measTrans.result
    Plus = measPlus.result

    f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\text{d}_\text{Tr}\left(\rho_1,\rho_2\right)",
    )
    ylims!(ax, 0, 1)

    # Plot the trace distances
    lines!(ax, ts, _trace_distance.(Up, Down), label = L"\text{d(Up,Down)}")
    lines!(ax, ts, _trace_distance.(Up, Plus), label = L"\text{d(Up,Plus)}")
    lines!(ax, ts, _trace_distance.(Up, Trans), label = L"\text{d(Up,Trans)}")
    lines!(ax, ts, _trace_distance.(Down, Plus), label = L"\text{d(Down,Plus)}")
    lines!(ax, ts, _trace_distance.(Down, Trans), label = L"\text{d(Down,Trans)}")
    lines!(ax, ts, _trace_distance.(Plus, Trans), label = L"\text{d(Plus,Trans)}")

    axislegend(position = :rt)

    f
end

# ─────────────────────────────────────────────────────────────
# Thermodynamics
# - Effective Hamiltonian Ks
# ─────────────────────────────────────────────────────────────

function plot_Ks(dirdata::String, row_idx, col_idx; tmax = nothing)
    effective_hamiltonian = computeKs(dirdata::String)
    config = loadconfig(dirdata*"/config.JSON")
    # xs: time 
    ts = effective_hamiltonian.time
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
        xlabel = L"\omega_c t",
        ylabel = L"\text{Re}(K_{%$(row_idx)%$(col_idx)})",
    )
    lines!(ax1, ts, ReKs, label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])")
    axislegend(position = :rt)
    # Imaginary part of Ks
    ax2 = Axis(
        f[2, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\text{Im}(K_{%$(row_idx)%$(col_idx)})",
    )
    lines!(ax2, ts, ImKs, label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])")
    axislegend(position = :rt)

    f
end
