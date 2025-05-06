"""
Save a collection of figures and animations for a given dataset.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `tomostates::Vector{String}`: List of tomography states to plot (default `["Up","Down","Plus","Trans"]`).
- `sdf::Bool`: Whether to save spectral density plots (default `true`).
- `Ks::Bool`: Whether to save effective Hamiltonian plots (default `true`).
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
    tfreqs = true,
    chain = true,
    envmodes = true,
)
    outdir = "figs-unbiased-500-hf/$dirdataname"
    mkpath(outdir)
    dirdata = "data-unbiased-500/$dirdataname"
    filename = dirdata * "/config.json"
    if !isempty(tomostates)
        for state in tomostates
            save("$outdir/tomo$state.png", plot_state(dirdata; state = state))
            save("$outdir/tomo_trdistance.png", plot_tomo_trdistance(dirdata))
        end
    end
    if sdf
        save("$outdir/sdf.png", plot_sdf(filename))
        save("$outdir/sdf_therm.png", plot_thermal_sdf(filename))
        save("$outdir/chain_coefficients.png", plot_chain_coefficients(filename))
    end
    if Ks
        save("$outdir/Ks00.png", plot_Ks(dirdata, 0, 0))
        save("$outdir/Ks01.png", plot_Ks(dirdata, 0, 1))
    end
    if tfreqs
        save("$outdir/trans_freqs.png", plot_transitionfreqs(dirdata))
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
- `tmax::AbstractFloat`: Maximum time cutoff to apply to the Ks plots.

# Outputs
Saves `Ks00_tmax_<tmax>.png` and `Ks01_tmax_<tmax>.png` under `figs/<dirdataname>`.
"""
function save_Ks(dirdataname::String, tmax::AbstractFloat)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    save("$outdir/Ks00_tmax_$tmax.png", plot_Ks(dirdata, 0, 0; tmax = tmax))
    save("$outdir/Ks01_tmax_$tmax.png", plot_Ks(dirdata, 0, 1; tmax = tmax))
end

"""
Save effective transition freqency plot with a specified time cutoff.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `tmax::AbstractFloat`: Maximum time cutoff to apply to the Ks plots.

# Outputs
Saves `tfreqs_tmax_<tmax>.png` under `figs/<dirdataname>`.
"""
function save_tfreqs(dirdataname::String, tmax::AbstractFloat)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    save("$outdir/trans_freqs_tmax_$tmax.png", plot_transitionfreqs(dirdata; tmax = tmax))
end

"""
Generate an animation of environment-mode occupations within specified frequency bounds.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `xmin::AbstractFloat`: Minimum frequecny bound for the animation.
- `xmax::AbstractFloat`: Maximum frequency bound for the animation.

# Outputs
Writes the animation files into `figs/<dirdataname>`.
"""
function save_envmodes(dirdataname::String, xmin::AbstractFloat, xmax::AbstractFloat)
    outdir = "figs/$dirdataname"
    mkpath(outdir)
    dirdata = "data/$dirdataname"
    animate_envmodes_occupation(dirdata, outdir; xmin = xmin, xmax = xmax)
end
