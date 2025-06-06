function prepareState(
   nchain::Integer,
   oscDim::Integer; 
   state::AbstractString="Up")

    #Environment is always the same
    #NOTE:n-index set consistently with single site system
    env = [Index(oscDim,"Boson,Site,n=$(i+1)") for i in 1:nchain]
    #Set the initial state (vacuum)
    stateenv = ["0" for _ in 1:nchain]

    #Prepare an appo system: initial state according to parameters
    apposys = siteinds("S=1/2",1)
    appostate = [state !==nothing ? state : "Up"]
    res = productMPS(ComplexF64,vcat(apposys,env),vcat(appostate,stateenv))

    return vcat(apposys,env),res
end

"""
    τ is the integration step.
"""
function SpinBoson_evolution_TEBD(
    Gammas,
    Lambdas,
    s;
    ϵ,
    Δ,
    sysenvInt::String,
    ChainLength,
    tau,
    ttotal,
    measStep,
    cutoff = 1E-14,
    minBondDim = 5,
    maxBondDim = 100,
    freqs,
    coups,
    trackBondDim = false,
    trackN = false,
    state = "Up",
    dir = "data"
)

    # Total size
    ntot = ChainLength + 1

    # Define the gates as in the Trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2],
    # F odd gate, G even gate.

    gates = ITensor[]

    h_start =
        0.5 * ϵ * op("Z", s[1]) * op("Id", s[2]) +
        0.5 * Δ * op("X", s[1]) * op("Id", s[2]) +
        0.5 * freqs[1] * op("N", s[2]) * op("Id", s[1]) +
        coups[1] * op(sysenvInt, s[1]) * op("A", s[2]) +
        coups[1] * op("Adag", s[2]) * op(sysenvInt, s[1])

    push!(gates, exp(-im * tau / 2 * h_start)) # first one is always divided by 2

    for j = 2:(ntot-2)

        t = isodd(j) ? tau / 2 : tau

        s1 = s[j]
        s2 = s[j+1]

        hj =
            coups[j] * op("Adag", s1) * op("A", s2) +
            coups[j] * op("A", s1) * op("Adag", s2) +
            0.5 * freqs[j-1] * op("N", s1) * op("Id", s2) +
            0.5 * freqs[j] * op("N", s2) * op("Id", s1)

        Gj = exp(-im * t * hj)
        push!(gates, Gj)

    end

    t = isodd(ntot) ? tau / 2 : tau

    h_end =
        freqs[ntot-1] * op("N", s[ntot]) * op("Id", s[ntot-1]) +
        0.5 * freqs[ntot-2] * op("N", s[ntot-1]) * op("Id", s[ntot]) +
        coups[ntot-1] * op("Adag", s[ntot]) * op("A", s[ntot-1]) +
        coups[ntot-1] * op("Adag", s[ntot-1]) * op("A", s[ntot])

    push!(gates, exp(-im * t * h_end))

    # Open files to save data (streaming approach)
    # - time.dat for time step;
    # - meas_<state>.dat for measurements of X, Y, Z on the TLS;
    # - norm_<state>.dat for norm of the TLS;
    # - meas_N.dat for chain sites occupation;
    # - bondDims.dat for monitoring the bond dimension
    ioTv  = open(joinpath(dir, "time_$state.dat"), "w")
    ioTLS  = open(joinpath(dir, "meas_$state.dat"), "w")
    ioNormCheck = open(joinpath(dir, "norm_$state.dat"), "w")

    if trackN
        ioN = open(joinpath(dir, "meas_N.dat"), "w")
    end

    if trackBondDim
        ioTrackBondDim = open(joinpath(dir, "bondDims_$state.dat"), "w")
    end

    for (step, t) in enumerate(0.0:tau:ttotal)
        println("time: ", t)

        if (step - 1) % measStep == 0
            # time
            println(ioTv, t)
            flush(ioTv)

            # Pauli matrix measurements on TLS
            spinMeasures = Vector{ComplexF64}([])
            appo =
                noprime!(op("X", s[1]) * Gammas[1]) *
                Lambdas[1] *
                dag(Gammas[1] * Lambdas[1])
            push!(spinMeasures, scalar(appo))
            appo =
                noprime!(op("Y", s[1]) * Gammas[1]) *
                Lambdas[1] *
                dag(Gammas[1] * Lambdas[1])
            push!(spinMeasures, scalar(appo))
            appo =
                noprime!(op("Z", s[1]) * Gammas[1]) *
                Lambdas[1] *
                dag(Gammas[1] * Lambdas[1])
            push!(spinMeasures, scalar(appo))
            # Write as X,Y,Z
            writedlm(ioTLS, transpose(spinMeasures), ',')
            flush(ioTLS)

            # Norm measure
            appoNorm = Gammas[1] * Lambdas[1] * dag(Gammas[1] * Lambdas[1])
            writedlm(ioNormCheck, [scalar(appoNorm)], ',')
            flush(ioNormCheck)

            if trackN
                occMeasures = Vector{ComplexF64}([])
                for i = 2:ntot-1
                    appoOcc =
                        Lambdas[i-1] *
                        noprime!(op("N", s[i]) * Gammas[i]) *
                        Lambdas[i] *
                        dag(Lambdas[i-1] * Gammas[i] * Lambdas[i])
                    push!(occMeasures, scalar(appoOcc))
                end
                appoOcc =
                    Lambdas[ntot-1] *
                    noprime!(op("N", s[ntot]) * Gammas[ntot]) *
                    dag(Lambdas[ntot-1] * Gammas[ntot])
                push!(occMeasures, scalar(appoOcc))
                writedlm(ioN, transpose(occMeasures), ',')
                flush(ioN)
            end           

            if trackBondDim
                bondDim = [dim(inds(a)[1]) for a in Lambdas]
                writedlm(ioTrackBondDim, transpose(bondDim), ',')
                flush(ioTrackBondDim)

            end
        end

        Gammas, Lambdas =
            apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim = 5, maxBondDim = maxBondDim)

    end

    close(ioTv)
    close(ioTLS)
    close(ioNormCheck)

    if trackN
        close(ioN)
    end

    if trackBondDim
        close(ioTrackBondDim)
    end

    return Gammas, Lambdas

end


function apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim = 5, maxBondDim)

    #println("first odd")
    Gammas, Lambdas =
        apply_odd(gates, Gammas, Lambdas; cutoff, mindim = 1, maxBondDim = maxBondDim)
    #println("Gammas[1]",Gammas[1])
    #println("Gammas[2]",Gammas[2])
    #println("After first gate:")
    #println(Lambdas)
    #print_alternated_inds(Gammas, Lambdas)
    #pippo = readline()
    #println("the even step")
    Gammas, Lambdas =
        apply_even(gates, Gammas, Lambdas; cutoff, mindim = 1, maxBondDim = maxBondDim)

    #println("second odd")
    #println("Gammas[1]: ", Gammas[1])
    Gammas, Lambdas =
        apply_odd(gates, Gammas, Lambdas; cutoff, mindim = 1, maxBondDim = maxBondDim)
    #print_alternated_inds(Gammas, Lambdas) 
    return Gammas, Lambdas
end

# Application of ODD GATES
function apply_odd(gates, Gammas, Lambdas; cutoff, mindim, maxBondDim)

    #Number of links:
    #we are updating on links, so this is the relevant Number
    N = length(Lambdas)

    #Determine number of threads
    nt = Threads.nthreads()

    #For all variables instantiate a number of copies
    #equal to the number of threads

    #List of variables to pool

    #There are things that are already vectors.
    #The "problem" is that:
    #- sequential approach: we push one element at a time.
    #-multi-threading: we need to record things in the right place. 
    #This place is typically determined by the loop index

    #Try to define a single code that can be executed both
    #sequentially and in multi-threading:

    #########ICE###############
    #BACKUP in TEDB_alt copy.jl
    #########ICE###############

    #Loop index arrays
    LambdasOdd = Array{ITensor}(undef, N)
    #Site tensors are 1 more than link tensors
    GammasOdd = Array{ITensor}(undef, N + 1)



    #print_inds(gates)

    # compute the inverse just for the right Lambdas (even), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension
    # Q: perdo tanta informazione troncando?
    # Q: mindim non è compatibile, ma non capisco perchè
    LambdasInvEven = ITensor[]
    for l = 1:length(Lambdas)
        if (l % 2 == 0)
            #LambdaInv = ITensor(inds(Lambdas[l]))
            LambdaInv = diag_itensor(inds(Lambdas[l]))
            LambdaInv .= inv.(Lambdas[l])
            push!(LambdasInvEven, LambdaInv)
        end
    end

    # println(LambdasInvEven)

    Threads.@threads for i = 1:2:N

        #println(Threads.threadid())

        #println("odd", i)
        #Temp array "local" to  for scope => one for each thread
        SVDLeftIndices = Index[]

        #empty!(SVDLeftIndices)

        #Apply the two site gate

        if i == 1 #first site
            #println("first site")
            #Debug
            # @show Gammas[1]
            # @show Lambdas[1]
            # @show Gammas[2]
            # @show Lambdas[2]
            double_site_psi = Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
            #println("first site completed")
        elseif i == N #last site, if it is an odd one, otherwise it will not be evolved 
            #println("Last site")
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1]
        else
            #println("site",i)
            double_site_psi =
                Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
        end

        #  println("double_site_psi ", double_site_psi)
        #  println("gates[$i] ", gates[i])

        Res = gates[i] * double_site_psi

        #Debug
        # if i==1
        #     @show Res
        # end

        # println("Res ", Res)

        indices = inds(Res)

        #println(indices)

        #find the indices for the SVD
        SVDLeftIndices = find_inds(indices, i)
        # if i==1
        # println("first_site SVDLeftIdx",SVDLeftIndices)
        # end

        #NOTA: se uso mindim=mindim si rompre tutto!!
        A, S, B =
            length(SVDLeftIndices) == 1 ?
            svd(
                Res,
                SVDLeftIndices[1];
                lefttags = "u=$(i)",
                righttags = "v=$(i)",
                cutoff = cutoff,
                use_absolute_cutoff = true,
                maxdim = maxBondDim,
            ) :
            svd(
                Res,
                SVDLeftIndices[1],
                SVDLeftIndices[2];
                lefttags = "u=$(i)",
                righttags = "v=$(i)",
                cutoff = cutoff,
                use_absolute_cutoff = true,
                maxdim = maxBondDim,
            )

        # return to Vidal form
        A = i != 1 ? LambdasInvEven[Int((i - 1) / 2)] * A : A
        B = i != N ? B * LambdasInvEven[Int((i + 1) / 2)] : B

        #All computation done.
        #Now we store the result

        if (i != 1)
            LambdasOdd[i-1] = Lambdas[i-1]
            #push!(LambdasOdd, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        LambdasOdd[i] = S
        GammasOdd[i] = A
        GammasOdd[i+1] = B
        # push!(LambdasOdd, S)
        # push!(GammasOdd, A)
        # push!(GammasOdd, B)

    end


    if (N % 2 == 0) #N is even, so the number of sites is odd
        LambdasOdd[N] = Lambdas[N]
        GammasOdd[N+1] = Gammas[N+1]
        # push!(LambdasOdd, Lambdas[N])
        # push!(GammasOdd, Gammas[N+1])
    end

    #print_alternated_inds(GammasOdd, LambdasOdd)

    # unprime the indices, maybe it can be done in the for loop
    for i in eachindex(Lambdas)
        noprime!(LambdasOdd[i])
        noprime!(GammasOdd[i])
    end

    noprime!(GammasOdd[end])

    return GammasOdd, LambdasOdd

end


# Application of EVEN GATES
function apply_even(gates, Gammas, Lambdas; cutoff, mindim, maxBondDim)


    N = length(Lambdas)

    #println("Even apply")

    LambdasInvOdd = ITensor[]
    LambdasEven = Array{ITensor}(undef, N)
    GammasEven = Array{ITensor}(undef, N + 1)

    # compute the inverse just for the right Lambdas (odd), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension 
    for l = 1:2:length(Lambdas)
        #LambdaInv = ITensor(inds(Lambdas[l]))
        #LambdaInv .= inv.(Lambdas[l])
        LambdaInv = diag_itensor(inds(Lambdas[l]))
        LambdaInv .= (Lambdas[l]) .^ (-1)
        push!(LambdasInvOdd, LambdaInv)
    end

    #println(LambdasInvOdd)

    #there is no first site to evolve, so the gamma is the same
    GammasEven[1] = Gammas[1]
    #push!(GammasEven, Gammas[1])

    Threads.@threads for i = 2:2:N

        #println("Even ",i)


        SVDLeftIndices = Index[]


        #there is no first site to evolve
        if i == N
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1]
        else
            double_site_psi =
                Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
        end

        Res = gates[i] * double_site_psi

        indices = inds(Res)

        #println(indices)

        #find the indices for the SVD
        SVDLeftIndices = find_inds(indices, i)

        #println(SVDLeftIndices)

        A, S, B =
            length(SVDLeftIndices) == 1 ?
            svd(
                Res,
                SVDLeftIndices[1];
                lefttags = "u=$(i)",
                righttags = "v=$(i)",
                cutoff = cutoff,
                use_absolute_cutoff = true,
                maxdim = maxBondDim,
            ) :
            svd(
                Res,
                SVDLeftIndices[1],
                SVDLeftIndices[2];
                lefttags = "u=$(i)",
                righttags = "v=$(i)",
                cutoff = cutoff,
                use_absolute_cutoff = true,
                maxdim = maxBondDim,
            )

        #println(inds(A))
        #println(inds(S))
        #println(inds(B))

        # return to Vidal form
        A = i != 1 ? LambdasInvOdd[Int(i / 2)] * A : A
        B = i != N ? B * LambdasInvOdd[Int((i + 2) / 2)] : B

        #println("gamma", i, ": ", inds(A))
        #println("gamma", (i + 1), ": ", inds(B))

        if (i != 1)
            LambdasEven[i-1] = Lambdas[i-1]
            #push!(LambdasEven, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        LambdasEven[i] = S
        GammasEven[i] = A
        GammasEven[i+1] = B
        # push!(LambdasEven, S)
        # push!(GammasEven, A)
        # push!(GammasEven, B)

    end


    if (N % 2 != 0) #N is odd, so the number of sites is even
        LambdasEven[N] = Lambdas[N]
        GammasEven[N+1] = Gammas[N+1]
        # push!(LambdasEven, Lambdas[N])
        # push!(GammasEven, Gammas[N+1])
    end

    #println("Even apply")
    #print_alternated_inds(GammasEven, LambdasEven)

    # unprime the indices
    for i in eachindex(Lambdas)
        noprime!(LambdasEven[i])
        noprime!(GammasEven[i])
    end

    noprime!(GammasEven[end])

    return GammasEven, LambdasEven

end

# Function to find the right indices (left bond one and site one) for the SVD
function find_inds(indices, i)
    SVDLeftIndices = Index[]
    for idx in indices
        if hastags(idx, "u=$(i-1)")

            push!(SVDLeftIndices, idx)
        end
        if hastags(idx, "n=$i")

            push!(SVDLeftIndices, idx)
        end
    end
    return SVDLeftIndices
end

# Defined to compute the expectation value using ITensors methods 
# left-canonical form
function convert_to_MPS(Gammas, Lambdas, ntot)

    psi = MPS(ntot)
    #println(length(psi))
    psi[1] = Gammas[1] #left normalized (A^dag A = id)

    for i in eachindex(Lambdas)
        psi[i+1] = Lambdas[i] * Gammas[i+1]
    end
    #println(psi)

    return psi
end
