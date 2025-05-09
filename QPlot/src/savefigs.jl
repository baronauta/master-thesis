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
function savefigs(dirdata::String)
    # All figures will be stored in a new directory, e.g. figs/ohmic_lowT
    dirdataname = basename(dirdata)
    outdir = "figs/$dirdataname"
    mkpath(outdir)

    # 
    for state in ["Up", "Down", "Plus", "Trans"]
        save("$outdir/tomo$state.png", plot_state(dirdata; state = state))
        save("$outdir/tomo_trdistance.png", plot_tomo_trdistance(dirdata))
    end

    # save("$outdir/chain_coefficients.png", plot_chain_coefficients(dirdata))

    save("$outdir/Ks00.png", plot_Ks(dirdata, 0, 0))
    save("$outdir/Ks01.png", plot_Ks(dirdata, 0, 1))
    save("$outdir/Ks11.png", plot_Ks(dirdata, 1, 1))
    save("$outdir/trans_freqs.png", plot_transitionfreqs(dirdata))
 
    animate_chain_occupation(dirdata, outdir)
    animate_envmodes_occupation(dirdata, outdir)
end

"""
Save effective Hamiltonian Ks plots with a specified time cutoff.

# Arguments
- `dirdataname::String`: Name of the data directory under `data/`.
- `tmax::AbstractFloat`: Maximum time cutoff to apply to the Ks plots.

# Outputs
Saves `Ks00_tmax_<tmax>.png` and `Ks01_tmax_<tmax>.png` under `figs/<dirdataname>`.
"""
function save_Ks(dirdata::String, tmax::AbstractFloat)
    dirdataname = basename(dirdata)
    mkpath("figs/$dirdataname")
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
function save_tfreqs(dirdata::String, tmax::AbstractFloat)
    dirdataname = basename(dirdata)
    mkpath("figs/$dirdataname")
    save("$outdir/trans_freqs_tmax_$tmax.png", plot_transitionfreqs(dirdata; tmax = tmax))
end
