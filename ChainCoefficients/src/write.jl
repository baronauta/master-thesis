"""
    write_chain_coefficients(dir::AbstractString, freqs::Vector, coups::Vector)

Write frequency and coupling coefficient arrays to CSV files under the given directory.
Returns the paths to the created frequency and coupling files.
"""
function write_chain_coefficients(dir::AbstractString, freqs::Vector, coups::Vector)
    freqs_file = joinpath(dir, "freqs.csv")
    coups_file = joinpath(dir, "coups.csv")
    # Write frequencies
    open(freqs_file, "w") do io
        for x in freqs
            writedlm(io, x, ',')
        end
    end
    # Write couplings
    open(coups_file, "w") do io
        for x in coups
            writedlm(io, x, ',')
        end
    end
    println("Chain coefficients saved in directory: '$dir'")
    return freqs_file, coups_file
end
