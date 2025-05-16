function convert_to_Vidal(psi)

    #begin from a right canonical form
    orthogonalize!(psi, 1)

    Gammas = ITensor[]
    Lambdas = ITensor[]
    N = length(psi)
    M = psi[1] #initialize the first M as the first site, as it is the orthogonality centre

    for i = 1:N-1

        #DEBUG: println("SVD$i")

        svd_inds = nothing

        svd_inds = uniqueinds(M, psi[i+1])

        A, Lambda, V = svd(M, svd_inds; lefttags = "u=$i", righttags = "v=$i")

        #DEBUG:
        #println("A: ", inds(A))
        #println("Lambda: ", inds(Lambda))
        #println("V: ", inds(V))
        #println(A * Lambda * dag(V) ≈ M)

        # push the Lambda matrix into the Lambdas vector  
        push!(Lambdas, Lambda)

        if i != 1 #do from the second site (NOTE: Lambdas will have N-1 elements) 
            # Assign the same indices (same ids) to the inverse of Lambdas, 
            # so they don't have just the same values but they can be contracted in the right way
            #LambdaInv = randomITensor(inds(Lambdas[i-1]))
            LambdaInv = diag_itensor(inds(Lambdas[i-1]))
            #println("debug ", LambdaInv)
            #println("debug2 ", Lambdas[i-1])
            #LambdaInv .= 0
            #println("debug3 pre assegnazione", LambdaInv)
            #LambdaInv .= Lambdas[i-1] .*2
            LambdaInv .= (Lambdas[i-1]) .^ (-1)
            #println("debug4 post assegnazione", LambdaInv)
            #CheckLambdaInv = prime(LambdaInv, commonind(Lambdas[i-1],LambdaInv))
            #println("debug check ", CheckLambdaInv)
            #ide = delta(noncommoninds(Lambdas[i-1],CheckLambdaInv))
            #println(" ide ", ide)
            #println(Lambdas[i-1]*CheckLambdaInv ≈ ide)
            # define the Gammas as per Vidal's form
            gamma = LambdaInv * A
        else #The first gamma is just the first A, 
            #note it is left normalized A A^dag =id (giusto?) 
            gamma = A
        end
        #push the gamma tensor into the Gammas vector
        push!(Gammas, gamma)

        if i != N - 1 #do unless it is the last site, 
            # in the last site there is no need to define M using the lambda
            M = Lambdas[i] * V * psi[i+1]
            # DEBUG: println(psi0)

        else #in the last site I leave out the Lambda, as it needs to stay out
            M = V * psi[i+1]
            # I don't need to make another SVD, so the M is the last gamma
            # Note that as V and psi[i+1] are both right orthogonal(RO), M is also RO
            push!(Gammas, M)
        end
    end

    return Lambdas, Gammas
end
