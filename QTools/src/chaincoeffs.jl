using Integrals

"""
    chaincoeff_jOhmic_hc(nummodes, β, α, s; nquad=20000)

Compute the frequency and coupling coefficients for a thermalized Ohmic quantum bath
based on the chain mapping of a spectral density `J(ω) = 2 α ωc (ω/ωc)^s θ(ω-ωc)`.

# Notes
This function wraps `MPSDynamics.chaincoeffs_finiteT` that assumes by default the Ohmic
spectral density function.
"""
function chaincoeff_jOhmic_hc(nummodes, β, α, s; nquad = 20000)
    # Compute chain coefficients: [frequencies, couplings, [κ₀]]
    println("Begin calculation of the chain coefficients")
    coeff = chaincoeffs_finiteT(
        nummodes,
        β,
        true; # consider Ohmic spectral density function `J(ω) = 2 α ωc (ω/ωc)^s θ(ω-ωc)`
        α = α,
        s = s,
        ωc = 1,
        Mmax = nquad,
        save = false,
    )
    println("Calculation finished")
    # Extract oscillator frequencies
    freqs = coeff[1]
    # Combine κ₀ (initial coupling) with subsequent chain couplings
    coups = vcat(coeff[3][1], coeff[2])
    return freqs, coups
end


function _integrateSDF(func, xmax::Float64)
    # Define a small value to approximate the domain boundaries near zero
    epsilon = 1e-15

    integrand(x, p) = func(x)

    # Right domain from zero to xmax
    right_domain = (epsilon, xmax)
    right_integral = solve(IntegralProblem(integrand, right_domain), QuadGKJL()).u

    # Left domain from -xmax to -epsilon
    left_domain = (-xmax, -epsilon)
    left_integral = solve(IntegralProblem(integrand, left_domain), QuadGKJL()).u

    # Return the sum of the left and right integrals
    return left_integral + right_integral
end

function check_chaincoeff(J, xmax::Float64)
    # κ0² 
    κ0_squared = _integrateSDF(J, xmax)
    # α0
    numα0(x) = J(x) * x
    denα0(x) = J(x)
    α0 = _integrateSDF(numα0, xmax) / _integrateSDF(denα0, xmax)

    # π1
    MOP1(x) = x - α0
    MOP1_squared(x) = MOP1(x)^2

    # α1
    numα1(x) = J(x) * x * MOP1_squared(x)
    denα1(x) = J(x) * MOP1_squared(x)
    α1 = _integrateSDF(numα1, xmax) / _integrateSDF(denα1, xmax)

    # β1
    numβ1(x) = J(x) * MOP1_squared(x)
    denβ1(x) = J(x)
    β1 = _integrateSDF(numβ1, xmax) / _integrateSDF(denβ1, xmax)

    # π2
    MOP2(x) = (x - α1) * MOP1(x) - β1
    MOP2_squared(x) = MOP2(x)^2

    # α2
    numα2(x) = J(x) * x * MOP2_squared(x)
    denα2(x) = J(x) * MOP2_squared(x)
    α2 = _integrateSDF(numα2, xmax) / _integrateSDF(denα2, xmax)

    # β2
    numβ2(x) = J(x) * MOP2_squared(x)
    denβ2(x) = J(x) * MOP1_squared(x)
    β2 = _integrateSDF(numβ2, xmax) / _integrateSDF(denβ2, xmax)

    # π3
    MOP3(x) = (x - α2) * MOP2(x) - β2 * MOP1(x)
    MOP3_squared(x) = MOP3(x)^2

    # α3
    numα3(x) = J(x) * x * MOP3_squared(x)
    denα3(x) = J(x) * MOP3_squared(x)
    α3 = _integrateSDF(numα3, xmax) / _integrateSDF(denα3, xmax)

    # β3
    numβ3(x) = J(x) * MOP3_squared(x)
    denβ3(x) = J(x) * MOP2_squared(x)
    β3 = _integrateSDF(numβ3, xmax) / _integrateSDF(denβ3, xmax)

    return [α0, α1, α2, α3], sqrt.([κ0_squared, β1, β2, β3])
end

function check_chaincoeff(dict::Dict{String,Any}, ωmax::Float64)
    # Parse sdf from dict and built the callable function
    J = read_thermalized_sdf(dict)
    # Compute first chain coefficients by integral calculus
    my_freqs, my_coups = check_chaincoeff(J, ωmax)
    return my_freqs, my_coups
end
