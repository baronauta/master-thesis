function maybe_save(fig; path::String, save::Bool=false)
    if save 
        @info "Saving figure to $path"
        CairoMakie.save(path, fig)
    else
        fig
    end
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

