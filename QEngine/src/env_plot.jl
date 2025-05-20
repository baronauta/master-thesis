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
# ─────────────────────────────────────────────────────────────


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

function animate_chain(time::Vector{Float64}, measN::Matrix{Float64}, β::Float64, s::Float64, outdir::AbstractString)

    # 1. Collect data
    sites = 1:size(measN,2)
    ts = round.(time, digits = 1)
    # Stepping function that returns the new data
    function progress_for_one_step!(i, ys, labels)
        i += 1
        return i, ys[i,:], labels[i]
    end

    # 2. Initialize the `Observable`s of the animation.

    # Static data -> no `Observable`
    xs = sites

    # Dynamic data -> `Observable`
    obs_ys = Observable(measN[1,:])
    obs_time = Observable(ts[1])
    # Text for the title, using `@lift` it will be updated runtine.
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
    ymax = maximum(maximum(measN)) * 11 / 10 #  Add 10% buffer for visibility.
    ylims!(ax, 0, ymax)
    # Plot
    lines!(
        ax,
        xs,
        obs_ys,
        color = :dodgerblue,
        label = L"\beta=%$(β),\,s=%$(s)",
    )
    # Fill space under the plot
    band!(ax, xs, zeros(length(xs)), obs_ys, color = (:dodgerblue, 0.2))

    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ys, labels, obs_ys, obs_labels)
        i, newys, newlabels = progress_for_one_step!(i, ys, labels)
        obs_ys[] = newys
        obs_labels[] = newlabels
    end

    # 5. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "$outdir/chain_occupation.mp4", frames; framerate = 60) do i
        animstep!(i, measN, ts, obs_ys, obs_time)
    end

end

function animate_envmodes(time::Vector{Float64}, modes::Vector{Float64}, occupations::Matrix{Float64}, J, ωs::Float64, β::Float64, s::Float64, outdir::AbstractString; xmin=-1, xmax=1)

    # 1. Collect data
    ts = round.(time, digits = 1)
    # Thermalized sdf for reference
    support = [minimum(modes), maximum(modes)]
    sdfx = collect(range(support..., 1000))
    sdfy = [ J(x) for x in sdfx ]
    scaled_sdfy = sdfy ./ maximum(sdfy) .* (0.5 * maximum(occupations[400,:]))
    
    # Stepping function that returns the new data
    function progress_for_one_step!(i, ys, labels)
        i += 1
        return i, ys[i,:], labels[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    # Dynamic data -> `Observable`
    obs_ys = Observable(occupations[1,:])
    obs_time = Observable(ts[1])
    # Text for the title, using `@lift` it will be updated runtine.
    text = @lift("t = $($obs_time)")

    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"\omega",
        ylabel = L"\langle n_\omega \rangle",
        title = text,
    )
    # Set ylims considering the maximum value in `ns`
    ymax = maximum(maximum(occupations)) * 11 / 10 #  Add 10% buffer for visibility.
    ylims!(ax, 0, ymax)
    # Plot
    lines!(
        ax,
        modes,
        obs_ys,
        color = :darkorange1,
        label = L"\beta=%$(β),\,s=%$(s)",
    )
    # Fill space under the plot
    band!(ax, modes, zeros(length(modes)), obs_ys, color = (:darkorange1, 0.2))
    # Sdf for reference
    lines!(ax, sdfx, scaled_sdfy, color = :gray, label = L"J(\omega)")
    # Vertical lines indicating eigenvalues of the bare system Hs
    vlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\pm\omega_s",
    )
    vlines!(ax, -ωs, color = (:purple, 0.8), linewidth = 1, linestyle = :dash)

    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ys, labels, obs_ys, obs_labels)
        i, newys, newlabels = progress_for_one_step!(i, ys, labels)
        obs_ys[] = newys
        obs_labels[] = newlabels
    end

    # 5. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "$outdir/envmodes_occupation_xmin_$(xmin)_xmax_$(xmax).mp4", frames; framerate = 60) do i
        animstep!(i, occupations, ts, obs_ys, obs_time)
    end

end


function animate_envmodes(time::Vector{Float64}, modes::Vector{Float64}, occupations::Matrix{Float64}, J, ωs::Float64, θs::Vector{Float64}, β::Float64, s::Float64, outdir::AbstractString; xmin=-1, xmax=1)

    # 1. Collect data
    # xs (modes), ys (occupations) to be updated
    ts = round.(time, digits = 1) # Labels, to be updated
    # ± ωs transition frequency for Hs
    freq = θs           # positive transition frquency Ks, to be updated
    negfreq = -freq     # negative transition frquency Ks, to be updated
    # Thermalized sdf for reference
    support = [minimum(modes), maximum(modes)]
    sdfx = collect(range(support..., 1000))
    sdfy = [ J(x) for x in sdfx ]
    scaled_sdfy = sdfy ./ maximum(sdfy) .* (0.5 * maximum(ns[10,:]))

    function progress_for_one_step!(i, ys, labels, vs, negvs)
        i += 1
        return i, ys[i, :], labels[i], vs[i], negvs[i]
    end

    # 2. Initialize the `Observable`s of the animation.

    # Dynamic data -> `Observable`
    obs_ys = Observable(occupations[1,:])
    obs_time = Observable(ts[1])
    obs_freq, obs_negfreq = Observable(freq[1]), Observable(negfreq[1])
    # Text for the title, using `@lift` it will be updated runtine.
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
    ylims!(ax, 0, maximum[end,:] * 11 / 10) # add 10% for better visibility
    # Set x-axis limits
    xlims!(ax, xmin, xmax)
    
    lines!(
        ax,
        xs,
        obs_ys,
        color = :darkorange1,
        label = L"\beta=%$(β),\,s=%$(s)",
    )

    # Fill space under the plot
    band!(ax, modes, zeros(length(modes)), obs_ys, color = (:darkorange1, 0.2))
    # Sdf for reference
    lines!(ax, sdfx, scaled_sdfy, color = :gray, label = L"J(\omega)")
    # Vertical lines indicating eigenvalues of the bare system Hs
    vlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\pm\omega_s",
    )
    vlines!(ax, -ωs, color = (:purple, 0.8), linewidth = 1, linestyle = :dash)
    # Vertical lines indicating frequency transition of Ks
    vlines!(
        ax,
        obs_freq,
        color = (:green, 0.8),
        linewidth = 1,
        linestyle = :dash,
        label = L"\pm\theta_s",
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
    record(
        fig,
        "$outdir/envmodes_occupation_xmin_$(xmin)_xmax_$(xmax).mp4",
        frames;
        framerate = 60,
    ) do i
        animstep!(i, occupations, ts, freq, negfreq, obs_ys, obs_time, obs_freq, obs_negfreq)
    end

end

