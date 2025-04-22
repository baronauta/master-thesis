
"""
Save a collection of figures and animations for a given dataset.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `tomostates::Vector{String}`: List of tomography states to plot (default `["Up","Down","Plus","Trans"]`).
- `sdf::Bool`: Whether to save spectral density plots (default `true`).
- `Ks::Bool`: Whether to save K-correlation plots (default `true`).
- `chain::Bool`: Whether to generate the chain occupation animation (default `true`).
- `envmodes::Bool`: Whether to generate the environment modes animation (default `true`).

# Outputs
Saves PNG files and animations to `figs/<dirdataname>`.
"""
function savefigs(
    dirdataname::String;
    tomostates = ["Up", "Down", "Plus", "Trans"],
    sdf = true,
    Ks = true,
    chain = true,
    envmodes = true
)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    filename = dirdata * "/config.json"
    if !isempty(tomostates)
        for state in tomostates
            save("$outdir/$state.png", plot_state(dirdata; state=state))
            save("$outdir/tomostates_trdistance.png", plot_tomo_trdistance(dirdata))
        end
    end
    if sdf
        save("$outdir/sdf.png", plot_sdf(filename))
        save("$outdir/thermalized_sdf.png", plot_thermal_sdf(filename))
        save("$outdir/chain_coefficients.png", plot_chain_coefficients(filename))
    end
    if Ks
        save("$outdir/Ks00.png", plot_Ks(dirdata, 0, 0))
        save("$outdir/Ks01.png", plot_Ks(dirdata, 0, 1))
    end
    if chain
        animate_chain_occupation(dirdata, outdir)
    end
    if envmodes
        animate_envmodes_occupation(dirdata, outdir)
    end
end

"""
Save effective Hamiltonian Ks plots with a specified time cutoff.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `tmax::Float`: Maximum time cutoff to apply to the Ks plots.

# Outputs
Saves `Ks00_tmax_<tmax>.png` and `Ks01_tmax_<tmax>.png` under `figs/<dirdataname>`.
"""
function save_Ks(dirdataname::String, tmax::Float)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    save("$outdir/Ks00_tmax_$tmax.png", plot_Ks(dirdata, 0, 0; tmax = tmax))
    save("$outdir/Ks01_tmax_$tmax.png", plot_Ks(dirdata, 0, 1; tmax = tmax))
end

"""
Generate an animation of environment-mode occupations within specified frequency bounds.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `xmin::Float`: Minimum frequecny bound for the animation.
- `xmax::Float`: Maximum frequency bound for the animation.

# Outputs
Writes the animation files into `figs/<dirdataname>`.
"""
function save_envmodes(dirdataname::String, xmin::Float, xmax::Float)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    animate_envmodes_occupation(dirdata, outdir; xmin = xmin, xmax = xmax)
end