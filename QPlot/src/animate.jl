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


function animate_chain_occupation(dirdata::String)
    # 1. Collect the data to be plotted. Prefer a stepping manner
    data = chain_occupation(dirdata)
    sites = data.sites
    ns = data.ns
    # Time scaled by Δ / π (for labels)
    ts = data.time
    config = loadconfig(dirdata * "/config.json")
    labels = round.(_scalets(ts, config.Delta), digits=3)
    # Stepping function that returns the new data, here simply
    # progress the index `i`.
    function progress_for_one_step!(i, ns, ts)
        i += 1
        return i, ns[i], ts[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    xs = sites # This is static
    ys = Observable(ns[1])
    labels = Observable(ts[1])
    
    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure(); display(fig)
    ax = Axis(
        fig[1, 1],
        xlabel = L"\text{chain sites}\,i",
        ylabel = L"\langle n_i \rangle",
        title = L"\text{Ohmic}:\,\alpha=%$(config.a[1]),\,T=%$(config.temperature)"
    )
    # Set ylims considering the maximum value in `ns`.
    # Add 10% for better visualization.
    ymax = maximum(maximum(ns)) * 11/10
    ylims!(ax, 0, ymax)
    # Unfortunately string are immutable and label can't be update in such manner
    lines!(ax, xs, ys, color=:dodgerblue, label=L"t\Delta / \pi = %$(labels[])")
    band!(ax, xs, zeros(length(xs)), ys, color = (:dodgerblue, 0.2))
    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ns, ts, ys, labels)
        i, newys, newlabels = progress_for_one_step!(i, ns, ts)
        ys[] = newys
        labels[] = newlabels
    end

    # # 5. That's all to make the animated plot!
    # i = 1
    # while i < length(ts)
    #     animstep!(i, ns, ts, ys, labels)
    #     sleep(0.001)
    #     i += 1
    # end

    # 6. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "chain_occupation.mp4", frames; framerate = 60) do i
        animstep!(i, ns, ts, ys, labels)
    end

end

function animate_envmodes_occupation(dirdata::String)
    # 1. Collect the data to be plotted. Prefer a stepping manner.
    data = envmodes_occupation(dirdata)
    modes = data.modes
    ns = data.ns
    # Time scaled by Δ / π (for labels)
    ts = data.time
    config = loadconfig(dirdata * "/config.json")
    labels = round.(_scalets(ts, config.Delta), digits=3)
    # Stepping function that returns the new data, here simply
    # progress the index `i`.
    function progress_for_one_step!(i, ns, ts)
        i += 1
        return i, ns[i], ts[i]
    end

    # 2. Initialize the `Observable`s of the animation.
    xs = modes # This is static
    ys = Observable(ns[1])
    labels = Observable(ts[1])
    
    # 3. Plot the `Observable`s and any other static elements   
    fig = Figure(); display(fig)
    ax = Axis(
        fig[1, 1], 
        xlabel = L"\omega",
        ylabel = L"\langle n_\omega \rangle",
        title = L"\text{Ohmic}:\,\alpha=%$(config.a[1]),\,T=%$(config.temperature)"
    )
    # Set ylims considering the maximum value in `ns`.
    # Add 10% for better visualization.
    ymax = maximum(maximum(ns)) * 11/10
    ylims!(ax, 0, ymax)
    # Unfortunately string are immutable and label can't be update in such manner
    lines!(ax, xs, ys, color=:darkorange1, label=L"t\Delta / \pi = %$(labels[])")
    band!(ax, xs, zeros(length(xs)), ys, color = (:darkorange1, 0.2))
    axislegend(position = :rt)

    # 4. Create the "animation stepping function".
    function animstep!(i, ns, ts, ys, labels)
        i, newys, newlabels = progress_for_one_step!(i, ns, ts)
        ys[] = newys
        labels[] = newlabels
    end

    # # 5. That's all to make the animated plot!
    # i = 1
    # while i < length(ts)
    #     animstep!(i, ns, ts, ys, labels)
    #     sleep(0.001)
    #     i += 1
    # end

    # 6. Save in a .mp4 file
    frames = 1:length(ts)-1
    record(fig, "envmodes_occupation.mp4", frames; framerate = 60) do i
        animstep!(i, ns, ts, ys, labels)
    end

end

