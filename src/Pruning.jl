# Support files containing all functions for Track Hypothesis pruning

# M-best pruning
function MBest(Tracks::Array{Track,1},AssMat::AssociationMatrix,Cost::HypothesisCost;M::Int64=50)

    # Columns
    ncol = size(AssMat.TOMHT,2)

    # Output Vars (should I just being doing this by refernce to save time?)
    newassmat  = zeros(Int64,(M,ncol))
    newcostmat = zeros(Float64,(M,ncol))
    newtracks  = Array{Track,1}(undef,M)

    for ii = 1:M

        # Find minimum cost hypothesis and the measurment ID for that hypoth
        idx = argmin(Cost.CostMatrix[:,end])
        ID  = AssMat.TOMHT[idx,end]

        # Write hypothesis to output vars
        newassmat[ii,:]  = AssMat.TOMHT[idx,:]
        newcostmat[ii,:] = Cost.CostMatrix[idx,:]
        newtracks[ii]    = Tracks[idx]

        # Remove all incompatible tracks from input vars
        if ID > 0 # Remove all hypotheses that assign that detection
            for rmidx in findall(x -> x==ID , AssMat.TOMHT[:,end])
                AssMat.TOMHT    = AssMat.TOMHT[setdiff(:, rmidx), :]
                Cost.CostMatrix = Cost.CostMatrix[setdiff(:, rmidx), :]
                filter!(x -> x ≠ Tracks[rmidx], Tracks)
            end
        else # Only remove the one hypothesis (need to preserve other missed track hypotheses)
            AssMat.TOMHT    = AssMat.TOMHT[setdiff(:, idx), :]
            Cost.CostMatrix = Cost.CostMatrix[setdiff(:, idx), :]
            filter!(x -> x ≠ Tracks[idx], Tracks)
        end
    end

    return newtracks, AssociationMatrix(newassmat), HypothesisCost(newcostmat)

end # End MBest
