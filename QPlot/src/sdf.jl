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
    return π/2 * params.α * x^params.s / params.ωc^(params.s - 1) * exp(-x / params.ωc)
end

function ohmic_sdf(x, params::SpectralDensityParams, tls::TLSParams)
    scale = 2*sqrt(tls.ϵ^2 + tls.Δ^2) ^ params.s
    return π/2 * params.α/scale * x^params.s / params.ωc^(params.s - 1) * exp(-x / params.ωc)
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

function plot_ohmic_sdf(params::SpectralDensityParams, support)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [ohmic_sdf(x, params) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label=L"s = %$(params.s)")
    axislegend(position = :rt)
    return fig
end

function plot_ohmic_sdf(params::SpectralDensityParams, support, tls::TLSParams)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [ohmic_sdf(x, params, tls) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label=L"s = %$(params.s)")
    axislegend(position = :rt)
    return fig
end

function plot_thermal_ohmic_sdf(params::SpectralDensityParams, T, support)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [thermal_ohmic_sdf(x, params, T) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label=L"s = %$(params.s),\, T = %$T")
    axislegend(position = :rt)
    return fig
end

function plot_thermal_ohmic_sdf(params::SpectralDensityParams, T, support, tls::TLSParams)
    fig = Figure()
    xs = collect(range(support..., 1000))
    ys = [thermal_ohmic_sdf(x, params, T, tls) for x in xs]

    ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
    lines!(ax, xs, ys, label=L"s = %$(params.s),\, T = %$T")
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
