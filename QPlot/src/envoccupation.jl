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
# Note
# In QEngine are provided the function `chain_occupation` and
# `envmodes_occupation`. Since the computation of normal modes
# is time consuming, `envmodes_occupation` saves data into .dat
# files. Plot functions for occupation number of normal modes
# will read data from this files.
# ─────────────────────────────────────────────────────────────

function plot_chain_occupation(dirdata::String, t_scaled::Real)
    # Collect data
    data = chain_occupation(dirdata)
    sites = data.sites
    ns = data.ns
    # Time scaled by Δ / π (for labels)
    config = loadconfig(dirdata * "/config.json")
    ts = _scalets(data.time, config.Delta)

    # Find the index of the closest ts to t_scaled
    diffs = abs.(ts .- t_scaled)
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


function plot_envmodes_occupation(dirdata::String, t_scaled::Real)
    # Collect data
    (rawdata, header) = readdlm("$dirdata/envmodes_data.dat", ',', Float64, header = true)
    meastime = rawdata[:, 1]
    ns = rawdata[:, 2:end]
    (rawdata, header) = readdlm("$dirdata/envmodes_modes.dat", ',', Float64, header = true)
    modes = rawdata[:]
    # Thermalized spectral density function for reference
    sdfx, sdfy = thermal_sdf(dirdata * "/config.json")
    # Eigenvalues of Hs (bare system Hamiltonian)
    Delta = config.Delta
    epsilon = config.epsilon
    ωs = 2 * sqrt(epsilon^2 + Delta^2)

    # Time scaled by Δ / π (for labels)
    config = loadconfig(dirdata * "/config.json")
    ts = _scalets(data.time, config.Delta)

    # Find the index of the closest ts to t_scaled
    diffs = abs.(ts .- t_scaled)
    t_idx = argmin(diffs)  # index of closest time

    # Extract corresponding values
    xs = sites[t_idx]
    ys = ns[t_idx]
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
        color = (:purple, 0.5),
        linewidth = 1,
        linestyle = :dash,
        label = L"\text{eigs}(H_s)",
    )
    vlines!(ax, -ωs, color = (:purple, 0.5), linewidth = 1, linestyle = :dash)
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

    # 1. Collect data
    data = chain_occupation(dirdata)
    sites = data.sites
    ns = data.ns
    # Time scaled by Δ / π (for labels)
    config = loadconfig(dirdata * "/config.json")
    ts = round.(_scalets(data.time, config.Delta), digits = 1)
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
    text = @lift("tΔ/π = $($obs_time)")

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

function animate_envmodes_occupation(
    dirdata::String,
    outdir::String;
    xmin = nothing,
    xmax = nothing,
)

    # 1. Collect data
    (rawdata, header) = readdlm("$dirdata/envmodes_data.dat", ',', Float64, header = true)
    meastime = rawdata[:, 1]
    ns = rawdata[:, 2:end]
    (rawdata, header) = readdlm("$dirdata/envmodes_modes.dat", ',', Float64, header = true)
    modes = rawdata[:]
    # Time scaled by Δ / π (for labels)
    config = loadconfig(dirdata * "/config.json")
    ts = round.(_scalets(meastime, config.Delta), digits = 1)
    # Thermalized spectral density function for reference
    sdfx, sdfy = thermal_sdf(dirdata * "/config.json")
    # Eigenvalues of Hs (bare system Hamiltonian)
    Delta = config.Delta
    epsilon = config.epsilon
    ωs = 2 * sqrt(epsilon^2 + Delta^2)
    # Frequency transition of Ks
    Ks = computeKs(dirdata).Ks
    E = eigen.(Ks)
    eminus = [eig.values[1] for eig in E]
    eplus = [eig.values[2] for eig in E]
    # Stepping function that returns the new data, here simply
    # progress the index `i`.
    function progress_for_one_step!(i, ns, ts, eminus, eplus)
        i += 1
        return i, ns[i,:], ts[i], eminus[i], eplus[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    # Static data -> no `Observable`
    xs = modes
    # Dynamic data -> `Observable`
    obs_ys = Observable(ns[1,:])
    obs_time = Observable(ts[1])
    obs_eminus = Observable(eminus[1])
    obs_eplus = Observable(eplus[1])
    # Text for the tile, using `@lift` it will be updated runtine.
    text = @lift("tΔ/π = $($obs_time)")

    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"\omega",
        ylabel = L"\langle n_\omega \rangle",
        title = text,
    )

    # Set ylims considering the maximum value in `ns`
    ymax = maximum(ns) * 1.1  # Add 10% buffer for visibility
    ylims!(ax, 0, ymax)
    # Rescale thermalized sdf (it is just for reference)
    sdfy = sdfy ./ maximum(sdfy) .* 1 / 10 * ymax

    if xmin !== nothing && xmax !==nothing
        xlims!(ax, xmin, xmax)
    else
        xlims!(ax, minimum(xs), maximum(xs))
    end

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
    lines!(ax, sdfx, sdfy, color = :gray70)
    # Vertical lines indicating eigenvalues of the bare system Hs
    vlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\text{eigs}(H_s)",
    )
    vlines!(ax, -ωs, color = (:purple, 0.8), linewidth = 1, linestyle = :dash)
    # Vertical lines indicating frequency transition of Ks
    vlines!(
        ax,
        obs_eminus,
        color = (:green, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\text{eigs}(K_s)",
    )
    vlines!(ax, obs_eplus, color = (:green, 0.8), linewidth = 1, linestyle = :dash)

    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ns, ts, eminus, eplus, obs_ys, obs_time, obs_eminus, obs_eplus)
        i, newys, newtime, neweminus, neweplus = progress_for_one_step!(i, ns, ts, eminus, eplus)
        obs_ys[] = newys
        obs_time[] = newtime
        obs_eminus[] = neweminus
        obs_eplus[] = neweplus
    end

    # 5. Save in a .mp4 file
    if xmin !== nothing && xmax !==nothing
        outpath = "$outdir/envmodes_occupation_zoom_$(xmin)_$(xmax).mp4"
    else
        outpath = "$outdir/envmodes_occupation.mp4"
    end
    frames = 1:length(ts)-1
    record(fig, filepath, frames; framerate = 60) do i
        animstep!(i, ns, ts, eminus, eplus, obs_ys, obs_time, obs_eminus, obs_eplus)
    end

end
