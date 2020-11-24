# Support files containing all functions for Track Hypothesis pruning

# M-best pruning
function MBest(Tracks::Array{Track,1},AssMat::AssociationMatrix,Cost::HypothesisCost,TuningParams::TuningParams)

    # M
    M = TuningParams.M

    # Columns
    assmat  = copy(AssMat.TOMHT)
    costmat = copy(Cost.CostMatrix)
    tracks  = copy(Tracks)
    ncol    = size(assmat,2)

    # Output Vars (should I just being doing this by refernce to save time?)
    newassmat = Array{Int64}(undef, 0, ncol)
    newcostmat = Array{Float64}(undef, 0, ncol)
    newtracks  = Array{Track,1}()

    while !isempty(assmat) && size(newassmat,1) < M
        # Find minimum cost hypothesis and the measurment ID for that hypoth
        idx = argmin(costmat[:,end])
        IDs = copy(assmat[idx,:])

        # Write hypothesis to output vars
        newassmat  = vcat(newassmat,assmat[idx,:]')
        newcostmat = vcat(newcostmat,costmat[idx,:]')
        newtracks  = vcat(newtracks,tracks[idx])

        # Remove all incompatible tracks from input vars
        for ii in 1:ncol
            ID  = IDs[ii]
            if ID > 0 # Remove all hypotheses that assign that detection
                rmidx   = findall(x -> x==ID , assmat[:,ii])
                assmat  = assmat[setdiff(1:end, rmidx), :]
                costmat = costmat[setdiff(1:end, rmidx), :]
                tracks  = tracks[setdiff(1:end, rmidx), :][:]
            end
        end
    end

    # Enforce there are at least M compatible tracks
    M = min(M,length(newtracks))

    return newtracks[1:M], AssociationMatrix(newassmat[1:M,:]), HypothesisCost(newcostmat[1:M,:])

end # End MBest

export MBest
