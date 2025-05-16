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

function chaincoeff_jOhmic_expc(nummodes, β, α, s; ωc=1.0, ωmax=10., nquad = 20000)
    # Spectral density base function (positive frequencies)
    sdf(x) = (π / 2) * α * ωc .* (x ./ ωc) .^ s .* exp.(-x ./ ωc)
    # Full spectral density function over two intervals (negative and positive frequencies)
    function J(x, i)
        if i == 1
            return -0.5 .* (1 .+ coth.(0.5 .* x .* β)) .* sdf.(-x)
        elseif i == 2
            return 0.5 .* (1 .+ coth.(0.5 .* x .* β)) .* sdf.(x)
        else
            println("Invalid interval index i = $i. Expected 1 or 2.")
        end
    end
    # Compute chain coefficients: [frequencies, couplings, [κ₀]]
    println("Begin calculation of the chain coefficients")
    coeff = chaincoeffs_finiteT(
        nummodes,
        β,
        false; 
        J = J,             # custom J(ω)
        mc = 2,            # interval
        AB = [[-ωmax 0];[0 ωmax]],     # defines the domain of J
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