struct SpectralDensityParams
    α::Float64
    s::Float64
    ωc::Float64
end

struct TLSParams
    ϵ::Float64
    Δ::Float64
end

function ohmic_sdf(x, params::SpectralDensityParams)
    return π / 2 * params.α * x^params.s / params.ωc^(params.s - 1) * exp(-x / params.ωc)
end

function ohmic_sdf(x, params::SpectralDensityParams, tls::TLSParams)
    scale = 2 * sqrt(tls.ϵ^2 + tls.Δ^2)^params.s
    return π / 2 * params.α / scale * x^params.s / params.ωc^(params.s - 1) *
           exp(-x / params.ωc)
end

function thermal_ohmic_sdf(x, params::SpectralDensityParams, T)
    if T == 0
        return x > 0 ? ohmic_sdf(x, params) : 0
    else
        return 0.5 * (1 + coth(0.5 * x / T)) * sign(x) * ohmic_sdf(abs(x), params)
    end
end

function thermal_ohmic_sdf(x, params::SpectralDensityParams, T, tls::TLSParams)
    if T == 0
        return x > 0 ? ohmic_sdf(x, params, tls) : 0
    else
        return 0.5 * (1 + coth(0.5 * x / T)) * sign(x) * ohmic_sdf(abs(x), params, tls)
    end
end

function compare_sdf(s; α_gatto = 0.1, α_me = 0.01, xmax = nothing)

    ϵ, Δ = 0.1, 0.0
    trans_freq = 2 * sqrt(ϵ^2 + Δ^2)
    scale = trans_freq^s

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    # Angela
    xs = collect(range([0, 1]..., 1000))
    ys = [2 * 0.01 / scale * x^s for x in xs]
    lines!(
        ax,
        xs,
        ys,
        label = L"J(\omega)=2\frac{0.01}{\omega_s^s} \omega^s \theta(\omega-1)\quad s=%$(s)\quad\mathrm{[Riva]}",
    )
    # Gatto
    xs = collect(range([0, 10]..., 1000))
    ys = [π / 2 * α_gatto * x / (x^2 + 1) for x in xs]
    lines!(
        ax,
        xs,
        ys,
        label = L"J(\omega)= \frac{\pi}{2} \alpha \omega \frac{1}{\omega^2+1}\quad\alpha=%$(α_gatto)\quad\mathrm{[Gatto]}",
    )
    # Me
    xs = collect(range([0, 10]..., 1000))
    ys = [π / 2 * α_me / scale * x^s * exp(-x) for x in xs]
    lines!(
        ax,
        xs,
        ys,
        label = L"J(\omega)= \frac{\pi}{2} \frac{\alpha}{\omega_s^s} \omega^s e^{-\omega}\quad\alpha=%$(α_me)\;s=%$(s)\quad\mathrm{[Me]}",
    )

    vlines!(
        ax,
        trans_freq,
        label = L"\omega_s=%$(trans_freq)",
        color = :gray,
        linestyle = :dash,
    )
    if !isnothing(xmax)
        xlims!(ax, 0, xmax)
    end
    axislegend(position = :rt)
    return fig
end

function plot_ohmic_sdf(params::SpectralDensityParams, support)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [ohmic_sdf(x, params) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label = L"s = %$(params.s)")
    axislegend(position = :rt)
    return fig
end

function plot_ohmic_sdf(params::SpectralDensityParams, support, tls::TLSParams)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [ohmic_sdf(x, params, tls) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label = L"s = %$(params.s)")
    axislegend(position = :rt)
    return fig
end

function plot_thermal_ohmic_sdf(params::SpectralDensityParams, T, support)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [thermal_ohmic_sdf(x, params, T) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label = L"s = %$(params.s),\, T = %$T")
    axislegend(position = :rt)
    return fig
end

function plot_thermal_ohmic_sdf(params::SpectralDensityParams, T, support, tls::TLSParams)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [thermal_ohmic_sdf(x, params, T, tls) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label = L"s = %$(params.s),\, T = %$T")
    axislegend(position = :rt)
    return fig
end


"""
    plot_chain_coefficients(filename::String)

Plot chain transformation frequencies and couplings from `filename` data.
Returns a `Figure`.
"""
function plot_chain_coefficients(dirdata::String)
    coefficients = chain_coefficients(dirdata * "/config.json")
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"n", ylabel = "Chain coefficients")
    lines!(
        ax,
        collect(1:length(coefficients.frequencies)),
        coefficients.frequencies,
        label = "Freqs",
    )
    lines!(
        ax,
        collect(1:length(coefficients.couplings)),
        coefficients.couplings,
        label = "Coups",
    )
    axislegend(position = :rb)
    f
end
