function plot_state(ts::Vector{Float64}, rx::Vector{Float64}, ry::Vector{Float64}, rz::Vector{Float64}, title::AbstractString; outdir::Union{Nothing, String} = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\omega_c t", title = L"\mathrm{%$(title)}")
    ylims!(ax, -1, 1)
    # Plot
    lines!(ax, ts, rx, label = L"\langle\sigma_x\rangle")
    lines!(ax, ts, ry, label = L"\langle\sigma_y\rangle")
    lines!(ax, ts, rz, label = L"\langle\sigma_z\rangle")
    axislegend(position = :rt)

    if isnothing(outdir)
        display(fig)
    else
        save(joinpath(outdir, "tomo$title.png"), fig)
    end

    return fig
end

# function _trace_distance(rho1, rho2)
#     # Compute the difference matrix (Hermitian)
#     delta = rho1 - rho2
#     # Compute eigenvalues (real since delta is Hermitian)
#     lambda = eigvals(delta)
#     # Compute trace distance using absolute eigenvalues
#     return 0.5 * sum(abs.(lambda))
# end

# function plot_tomo_trdistance(dirdata)

#     measUp = get_measurements(dirdata * "/measurements_Up.dat", "densitymatrix")
#     measDown = get_measurements(dirdata * "/measurements_Dn.dat", "densitymatrix")
#     measPlus = get_measurements(dirdata * "/measurements_+.dat", "densitymatrix")
#     measTrans = get_measurements(dirdata * "/measurements_i.dat", "densitymatrix")

#     # xs: time 
#     ts = measUp.time

#     # ys: density Matrix
#     Up = measUp.result
#     Down = measDown.result
#     Trans = measTrans.result
#     Plus = measPlus.result

#     f = Figure()
#     ax = Axis(
#         f[1, 1],
#         xlabel = L"\omega_c t",
#         ylabel = L"\text{d}_\text{Tr}\left(\rho_1,\rho_2\right)",
#     )
#     ylims!(ax, 0, 1)

#     # Plot the trace distances
#     lines!(ax, ts, _trace_distance.(Up, Down), label = L"\text{d(Up,Down)}")
#     lines!(ax, ts, _trace_distance.(Up, Plus), label = L"\text{d(Up,Plus)}")
#     lines!(ax, ts, _trace_distance.(Up, Trans), label = L"\text{d(Up,Trans)}")
#     lines!(ax, ts, _trace_distance.(Down, Plus), label = L"\text{d(Down,Plus)}")
#     lines!(ax, ts, _trace_distance.(Down, Trans), label = L"\text{d(Down,Trans)}")
#     lines!(ax, ts, _trace_distance.(Plus, Trans), label = L"\text{d(Plus,Trans)}")

#     axislegend(position = :rt)

#     f
# end

function plot_Ks(ts::Vector{Float64}, Ks::Vector{Matrix{ComplexF64}}, i::Integer, j::Integer, β::Float64, s::Float64; outdir::Union{Nothing, String} = nothing)
    # Real and Imag part of Effective Hamiltonian Ks
    ReKs = [real(K[i, j]) for K in Ks]
    ImKs = [imag(K[i, j]) for K in Ks]

    fig = Figure()
    # Real part of Ks
    ax1 = Axis(
        fig[1, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\text{Re}(K_{%$(i-1)%$(j-1)})",
    )
    lines!(ax1, ts, ReKs, label = L"\beta=%$(β),\,s=%$(s)")
    axislegend(position = :rt)
    # Imaginary part of Ks
    ax2 = Axis(
        fig[2, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\text{Im}(K_{%$(i-1)%$(j-1)})",
    )
    lines!(ax2, ts, ImKs, label = L"\beta=%$(β),\,s=%$(s)")
    axislegend(position = :rt)

    if isnothing(outdir)
        display(fig)
    else
        save(joinpath(outdir, "Ks$(i)$(j).png"), fig)
    end

    return fig
end

"""
    plot_transitionfreqs(time::Vector{Float64}, ωs::Float64, θs::Vector{Float64})

Compare the transition frequencies of bare Hamiltonian Hs and the effective Hamiltonian Ks.
"""
function plot_transitionfreqs(time::Vector{Float64}, ωs::Float64, θs::Vector{Float64}; outdir::Union{Nothing, String} = nothing)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"\omega_c t",
        ylabel = L"\mathrm{Transition frequency}",
    )
    # Hs transition freequency
    hlines!(
        ax,
        ωs,
        color = (:purple, 0.8),
        linestyle = :dash,
        label = L"\omega s",
    )
    # Ks transition freequency
    lines!(ax, time, θs, color = :green, label = L"θs")
    axislegend(position = :rt)

    if isnothing(outdir)
        display(fig)
    else
        save(joinpath(outdir, "tomo$state.png"), fig)
    end

    return fig
end