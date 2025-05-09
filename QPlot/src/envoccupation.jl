# ─────────────────────────────────────────────────────────────
# Environment Occupation Visualization
#
# This file provides functions for visualizing occupation numbers:
#
# • Fixed-Time Plots:
#   - Plot occupation number of chain sites at a chosen time
#   - Plot occupation number of normal modes at a chosen time
#
# • Time-Resolved Animations:
#   - Animate chain site occupations over time (.mp4)
#   - Animate normal mode occupations over time (.mp4)
#
# This file also contains functions for visualizing:
#   - Comparison between transition frequency of Hs and Ks
#
# Note
# In QEngine are provided the function `chain_occupation` and
# `envmodes_occupation`. Since the computation of normal modes
# is time consuming, `envmodes_occupation` saves data into .dat
# files. Plot functions for occupation number of normal modes
# will read data from this files.
# ─────────────────────────────────────────────────────────────

"""
    plot_chain_occupation(dirdata::String, t::Real)

Plot the chain site occupations ⟨nᵢ⟩ at a given scaled time `t` using data from `dirdata`.
"""
function plot_chain_occupation(dirdata::String, t::Real)
    # Collect data
    data = chain_occupation(dirdata)
    sites = data.sites
    ns = data.ns
    # Time (for labels)
    ts = data.time

    # Find the index of the closest ts to `t`
    diffs = abs.(ts .- t)
    t_idx = argmin(diffs)  # index of closest time

    # Extract corresponding values
    xs = sites[t_idx]
    ys = ns[t_idx]
    label_t = round(ts[t_idx]; digits = 1)

    # Plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\text{chain sites}\,i", ylabel = L"\langle n_i \rangle")
    lines!(
        ax,
        xs,
        ys,
        label = L"t\Delta / \pi = %$(label_t)\,T=%$(config.temperature),\,\text{s}=%$(config.a[3])",
    )
    axislegend(position = :rt)

    return f
end

"""
    plot_envmodes_occupation(dirdata::String, t::Real)

Plot environment mode occupations ⟨n_ω⟩ at a given scaled time `t`, with spectral density and frequency transitions.
"""
function plot_envmodes_occupation(dirdata::String, t::Real)
    # Collect data
    (rawdata, header) = readdlm("$dirdata/envmodes_data.dat", ',', Float64, header = true)
    meastime = rawdata[:, 1]
    ns = rawdata[:, 2:end]
    (rawdata, header) = readdlm("$dirdata/envmodes_modes.dat", ',', Float64, header = true)
    modes = rawdata[:]
    # Thermalized spectral density function for reference
    sdfx, sdfy = thermal_sdf(dirdata * "/config.json")
    # Frequency transition for Hs (bare system Hamiltonian):
    # Hs = ϵ σz + Δ σx → E± = √(ϵ^2+Δ^2) eigvals of Hs,
    # frequency transition is ωs = (E+ - E-)/ħ = 2√(ϵ^2+Δ^2).
    Delta = config.Delta
    epsilon = config.epsilon
    ωs = 2 * sqrt(epsilon^2 + Delta^2)
    # Frequency transition of Ks (effective Hamiltonian):
    # ωs' = (E+ - E-)/ħ with E± eigvals of Ks.
    freqs = _effective_freqs(dirdata, meastime)

    # Time (for labels)
    ts = data.time

    # Find the index of the closest ts to `t`
    diffs = abs.(ts .- t)
    t_idx = argmin(diffs)  # index of closest time

    # Extract corresponding values
    xs = sites[t_idx]
    ys = ns[t_idx]
    ωeff = freqs[t_idx]
    label_t = round(ts[t_idx]; digits = 1)

    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"\omega", ylabel = L"\langle n_\omega \rangle")
    # Modes occupation
    lines!(
        ax,
        xs,
        ys,
        label = L"t\Delta / \pi = %$(label_t)\,T=%$(config.temperature),\,\text{s}=%$(config.a[3])",
    )
    # Fill space under the plot
    band!(ax, xs, zeros(length(xs)), ys, color = (:darkorange1, 0.2))
    # Sdf for reference
    lines!(ax, sdfx, sdfy, color = :gray70)
    # Vertical lines indicating eigenvalues of the bare system Hs
    vlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\omega_s\quad (H_S)",
    )
    vlines!(ax, -ωs, color = (:purple, 0.8), linewidth = 1, linestyle = :dash)
    vlines!(
        ax,
        ωeff,
        color = (:green, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\omega_s'\quad (K_S)",
    )
    vlines!(ax, -ωeff, color = (:green, 0.8), linewidth = 1, linestyle = :dash)
    axislegend(position = :rt)

    f
end

"""
    plot_transitionfreqs(dirdata::String; tmax = nothing)

Plot a comparison between transition frequency of the bare Hamiltonian Hs and the effective Hamiltonian Ks.
"""
function plot_transitionfreqs(dirdata::String; tmax = nothing)
    config = loadconfig(dirdata * "/config.json")

    # Bare Hamiltonian Hs (time-independent frequency)
    freq_Hs = 2 * sqrt(config.epsilon^2 + config.Delta^2)

    # Effective Hamiltonian Ks (time-dependent frequency)
    (rawdata, header) = readdlm("$dirdata/envmodes_data.dat", ',', Float64, header = true)
    meastime = rawdata[:, 1]
    freqs_Ks = effective_freqs(dirdata, meastime)

    # xs: time
    ts = meastime

    if tmax !== nothing
        # Shorten time domain
        mask = ts .<= tmax
        ts = ts[mask]
        freqs_Ks = freqs_Ks[mask]
    end

    f = Figure()

    ax = Axis(
        f[1, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\text{Transition frequency}\;\omega",
    )
    # Hs eigenvalues
    hlines!(
        ax,
        freq_Hs,
        color = (:purple, 0.8),
        linestyle = :dash,
        label = L"\text{Bare Hamiltonian}\;H_s",
    )
    # Ks eigenvalues
    lines!(ax, ts, freqs_Ks, color = :green, label = L"\text{Effective Hamiltonian}\;K_s")
    axislegend(position = :rt)

    f
end



# ─────────────────────────────────────────────────────────────
# Animations using GLMakie
#
# See for reference 
# https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd
#
# %% 1. Collect data to be plotted. Prefer a stepping manner.
# The goal is to create a "step" function which
# once called it progresses the data for one animation frame.
#
# %% 2. Initialize the `Observable`s of the animation 
# If any plotted element is initialized as an observable, 
# then Makie.jl understands this. Updating the observable 
# would update the plot values.
#
# %% 3. Plot the `Observable`s and any other static elements
#
# %% 4. Create the "animation stepping function"
# Using the functions of step 1, we now define a function
# that updates the observables. Makie.jl understands observable
# updates and directly reflects this on the plotted elements.
#
# %% 5. Save in .mp4 file using `record`.
# ─────────────────────────────────────────────────────────────

function animate_chain_occupation(dirdata::String, outdir::String)

    config = loadconfig(dirdata * "/config.JSON")
    # 1. Collect data
    data = chain_occupation(dirdata)
    sites = data.sites
    ns = data.ns
    # Time
    ts = round.(data.time, digits = 1)
    # Stepping function that returns the new data, here simply
    # progress the index `i`.
    function progress_for_one_step!(i, ns, ts)
        i += 1
        return i, ns[i], ts[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    # Static data -> no `Observable`
    xs = sites
    # Dynamic data -> `Observable`
    obs_ys = Observable(ns[1])
    obs_time = Observable(ts[1])
    # Text for the tile, using `@lift` it will be updated runtine.
    text = @lift("t = $($obs_time)")

    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"\text{chain sites}\,i",
        ylabel = L"\langle n_i \rangle",
        title = text,
    )
    # Set ylims considering the maximum value in `ns`
    ymax = maximum(maximum(ns)) * 11 / 10 #  Add 10% buffer for visibility.
    ylims!(ax, 0, ymax)
    # Plot
    lines!(
        ax,
        xs,
        obs_ys,
        color = :dodgerblue,
        label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])",
    )
    # Fill space under the plot
    band!(ax, xs, zeros(length(xs)), obs_ys, color = (:dodgerblue, 0.2))

    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ns, ts, obs_ys, obs_time)
        i, newys, newtime = progress_for_one_step!(i, ns, ts)
        obs_ys[] = newys
        obs_time[] = newtime
    end

    # 5. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "$outdir/chain_occupation.mp4", frames; framerate = 60) do i
        animstep!(i, ns, ts, obs_ys, obs_time)
    end

end

function animate_envmodes_occupation(dirdata::String, outdir::String; xmin=-1, xmax=1)
    # 1. Collect data
    config = loadconfig(dirdata * "/config.JSON")
    # Normal modes (xs), modes occupation (ns) and time for labels (ts)
    modes, ns, meastime = read_envmodes_occupation(dirdata)
    ts = round.(meastime, digits = 1)
    # Get max value of the mode occupations
    max_ns = maximum([maximum(n) for n in ns])
    # Thermalized spectral density function for reference
    support = [-config.domain[2], config.domain[2]]
    sdfx = collect(range(support..., 1000))
    sdfy = [
        thermal_ohmic_sdf(
            x,
            SpectralDensityParams(config.a[1], config.a[3], config.a[2]),
            config.temperature,
        ) for x in sdfx
    ]
    # Scale sdfy to 50% of max_ns
    scaled_sdfy = sdfy ./ maximum(sdfy) .* (0.5 * max_ns)
    # Frequency transition for Hs (bare system Hamiltonian):
    # Hs = ϵ σz + Δ σx → E± = √(ϵ^2+Δ^2) eigvals of Hs,
    # frequency transition is ωs = (E+ - E-)/ħ = 2√(ϵ^2+Δ^2).
    ωs = 2 * sqrt(config.epsilon^2 + config.Delta^2)
    # Frequency transition of Ks (effective Hamiltonian):
    # ωs' = (E+ - E-)/ħ with E± eigvals of Ks.
    freq = effective_freqs(dirdata, meastime)
    negfreq = -freq

    # Stepping function that returns the new data, here simply
    # progress the index `i`.
    function progress_for_one_step!(i, ns, ts, freq, negfreq)
        i += 1
        return i, ns[i, :], ts[i], freq[i], negfreq[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    # Static data -> no `Observable`
    xs = modes
    # Dynamic data -> `Observable`
    obs_ys = Observable(ns[1, :])
    obs_time = Observable(ts[1])
    obs_freq, obs_negfreq = Observable(freq[1]), Observable(negfreq[1])
    # Text for the tile, using `@lift` it will be updated runtime
    text = @lift("t = $($obs_time)")

    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"\omega",
        ylabel = L"\langle n_\omega \rangle",
        title = text,
    )

    # Set y-axis limits
    ylims!(ax, 0, max_ns * 11 / 10) # add 10% for better visibility
    # Set x-axis limits
    xlims!(ax, xmin, xmax)

    lines!(
        ax,
        xs,
        obs_ys,
        color = :darkorange1,
        label = L"T=%$(config.temperature),\,\text{s}=%$(config.a[3])",
    )
    # Fill space under the plot
    band!(ax, xs, zeros(length(xs)), obs_ys, color = (:darkorange1, 0.2))
    # Sdf for reference
    lines!(ax, sdfx, scaled_sdfy, color = :gray, label=L"J(\omega)")
    # Vertical lines indicating eigenvalues of the bare system Hs
    vlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\text{Transition frequency } H_S",
    )
    vlines!(ax, -ωs, color = (:purple, 0.8), linewidth = 1, linestyle = :dash)
    # Vertical lines indicating frequency transition of Ks
    vlines!(
        ax,
        obs_freq,
        color = (:green, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\text{Transition frequency } K_S",
    )
    vlines!(ax, obs_negfreq, color = (:green, 0.8), linewidth = 1, linestyle = :dash)

    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ns, ts, freq, negfreq, obs_ys, obs_time, obs_freq, obs_negfreq)
        i, newys, newtime, newfreq, newnegfreq =
            progress_for_one_step!(i, ns, ts, freq, negfreq)
        obs_ys[] = newys
        obs_time[] = newtime
        obs_freq[], obs_negfreq[] = newfreq, newnegfreq
    end

    # 5. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "$outdir/envmodes_occupation_xmin_$(xmin)_xmax_$(xmax).mp4", frames; framerate = 60) do i
        animstep!(i, ns, ts, freq, negfreq, obs_ys, obs_time, obs_freq, obs_negfreq)
    end

end
