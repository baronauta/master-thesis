using DelimitedFiles

# --- Write in .dat files ---

function write_vector(dir::String, filename::String, vec::AbstractVector)
    file = joinpath(dir, filename * ".dat")
    open(file, "w") do io
        writedlm(io, vec, ',')
    end
    return file
end

function write_matrix(dir::String, filename::String, mat::AbstractMatrix)
    file = joinpath(dir, filename * ".dat")
    open(file, "w") do io
        writedlm(io, mat, ',')
    end
    return file
end

function write_3darray(
    dir::String,
    filename::String,
    arr::AbstractArray{<:Any,3},
)
    file = joinpath(dir, filename * ".dat")
    open(file, "w") do io
        for i = 1:size(arr, 1)
            # Flatten arr[i, :, :] in row-major order
            flat_slice = vec(permutedims(arr[i, :, :], (1, 2)))
            # Convert to strings and write as comma-separated values
            println(io, join(string.(flat_slice), ","))
        end
    end
    return file
end

function read_3darray(filepath::String, N::Int)
    lines = readlines(filepath)
    M = length(lines)  # Number of 2D slices (i.e., arr[i, :, :])
    arr = Array{ComplexF64}(undef, M, N, N)

    for (i, line) in enumerate(lines)
        values = parse.(ComplexF64, split(line, ","))
        arr[i, :, :] = reshape(values, (N, N))
    end

    return arr
end


# data = readdlm(joinpath(dir, "measEnv_$(state).dat"), ',', ComplexF64)
# arr_reconstructed = Array{ComplexF64}(undef, size(data, 1), N, N)
# for i in 1:size(data, 1)
#     arr_reconstructed[i, :, :] = reshape(data[i, :], (N, N))
# end
# arr_reconstructed[10,:,:] ≈ cdagc[10,:,:]
# all([cdagc[i,:,:] ≈ arr_reconstructed[i,:,:] for i in 1:size(cdagc, 1)])

# function read_Ks(dir::String)
#     data = readdlm(joinpath(dir, "Ks_matrices.dat"), ',', ComplexF64, '\n')
#     Ks = [
#         reshape(ComplexF64.(row), 2, 2)
#         for row in eachrow(data)
#     ]
#     return Ks
# end
