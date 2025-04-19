function savefigs(
    dirdataname::String;
    tomostates = false,
    sdf = false,
    Ks = false,
    chain = false,
    envmodes = false,
    tmax = Inf,
    xmin = nothing,
    xmax = nothing
)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    filename = dirdata * "/config.json"
    if tomostates == true
        save("$outdir/tomostates.png", plot_tomostates(dirdata))
        save("$outdir/tomostates_trdistance.png", plot_tomo_trdistance(dirdata))
    end
    if sdf == true
        save("$outdir/sdf.png", plot_sdf(filename))
        save("$outdir/thermalized_sdf.png", plot_thermal_sdf(filename))
        save("$outdir/chain_coefficients.png", plot_chain_coefficients(filename))
    end
    if Ks == true
        save("$outdir/Ks00.png", plot_Ks(dirdata, 0, 0))
        save("$outdir/Ks01.png", plot_Ks(dirdata, 0, 1))
        if tmax !== Inf
            save("$outdir/Ks00_tmax_$tmax.png", plot_Ks(dirdata, 0, 0; tmax = tmax))
            save("$outdir/Ks01_tmax_$tmax.png", plot_Ks(dirdata, 0, 1; tmax = tmax))
        end
    end
    if chain == true
        animate_chain_occupation(dirdata, outdir)
    end
    if envmodes == true
        animate_envmodes_occupation(dirdata, outdir; xmin = xmin, xmax = xmax)
    end
end
