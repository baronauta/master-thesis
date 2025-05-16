using Revise
using Pkg
Pkg.activate("./QPlot")

using JSON
using ChainCoefficients
using QPlot

β = 2.0             # β=2.0 (T=0.5), β=2.0 (T=0.0005)

# α
# rescaled by transition frequency of Hs = ϵ/2 σz + Δ/2 σx
ϵ = 0.2
Δ = 0.0
α = 0.01 / sqrt(ϵ^2+Δ^2)^s

# s
s = 2.0

J(x) = jOhmic_hc(x, α, s; ωc = 1)
tJ(x) = therm_jOhmic(x, J, β)

fig = Figure()
xs = collect(range([-1,1]..., 1000))
ys = [tJ(x) for x in xs]
ax = Axis(fig[1, 1], xlabel = L"\omega", ylabel = L"J(\omega)")
lines!(ax, xs, ys)
display(fig)