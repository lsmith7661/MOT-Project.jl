# Take in image and return detections

function PixelThreshold(Img;kmax::Int=20,T::Int=100, Time::Real=1)
    # Img     - m1xm2 array of intensities, an image
    # k_max   - maximum k number of clusters
    # T       - Threshold
    # mean of noise pdf assumed 0

    # Image size
    𝓂1 = size(Img)[1] #Vertical
    𝓂2 = size(Img)[2] #Horizontal

    # Tracking Window (3σ each diection)
    σ_s = sqrt(var(Img)) # Could try to calculate from known camera properties
    W_m1 = 2*ceil(3*σ_s)
    W_m2 = W_m1;

    # Detections
    pᵢ = findall(x -> x>T,Img)

    # K-means Clustering
    pk = [last.(Tuple.(pᵢ)) first.(Tuple.(pᵢ))]' # Cartesian index results of kmeans are horz x vert - need vert x horz for consistancy
    dist = Euclidean()
    dmat = pairwise(dist, pk, dims=2)
    dsum = zeros(1,kmax)
    C = Array{Any,1}(undef,kmax)
    for k = 2:kmax
        C[k] = kmeans(pk,k,maxiter=200)
        dsum[k] = sum(silhouettes(C[k],dmat)) #dist(C.centers[:,C.assignments[ii]], pk[:,ii])
    end
    k = argmax(dsum)[2]

    # Initialize Output array type
    d = Array{Detection}(undef,k)

    # Center of Each Cluster (in an array)
    k_centers = [(C[k].centers[1,col],C[k].centers[2,col]) for col in 1:size(C[k].centers,2)]

    # Image Heatmap
    #plt1 = heatmap(Img, color = :greys, title="K-means Results")

    # Tracking Window of Each Cluster (vert,horz)
    Wv = []
    Wh = []
    for i in 1:size(k_centers,1)
        push!(Wv,clamp(round(k_centers[i][1]-W_m1/2),1,𝓂1-W_m1))
        push!(Wh,clamp(round(k_centers[i][2]-W_m2/2),1,𝓂2-W_m2))
    end

    # Plot Tracking Windows
    #for i in 1:length(Wx)
    #    plot!(plt1,rectangle(W_m1,W_m2,Wx[i],Wy[i]), linecolor=:green, linestyle=:dash, fillalpha=0)
    #end

    # Display K-Means
    #display(plot!(plt1,k_centers, seriestype=:scatter, marker=:star, markersize=3, legend=false))

    # Loop over each tracking window
    for ii in 1:k
        # Create Subimage from Image
        Imgᵢ = Img[Int(Wv[ii]):Int(Wv[ii])+Int(W_m1),Int(Wh[ii]):Int(Wh[ii])+Int(W_m2)]

        # Extract Signal
        s = 0
        for i in CartesianIndices(Imgᵢ)
            s += Imgᵢ[i]
        end
        #@show s

        # Signal Uncertainty
        𝓂   = W_m1*W_m2
        σn   = σ_s # this σ_s is "signal diffraction width" which is different than the signal uncertainty
        σ_s² = (𝓂*σn^2)
        #@show σ_s²

        # Extract Centroid
        ∑ζI = 0
        ∑ηI = 0
        ∑I  = 0
        for i in CartesianIndices(Imgᵢ)
            ∑I  += Imgᵢ[i]
            ∑ζI += i[2]*Imgᵢ[i] # cartesian indicies are (horz x vert)
            ∑ηI += i[1]*Imgᵢ[i] # cartesian indicies are (horz x vert)
        end
        ζ̂c = ∑ζI/∑I + Wv[ii] # vert
        η̂c = ∑ηI/∑I + Wh[ii] # horz
        #@show ζ̂c,η̂c

        # Centroid Uncertainty
        ∑ζ  = 0
        ∑η  = 0
        for i in CartesianIndices(Imgᵢ)
            ∑ζ  += (i[2]-(W_m1/2)-.5)^2 # 0.5 accounts for index starting at 1
            ∑η  += (i[1]-(W_m2/2)-.5)^2 # 0.5 accounts for index starting at 1
        end
        σ₁²   = (𝓂*σn^2)/(s^2 + 𝓂*σn^2)
        σ₂_ζ² = (∑ζ*σn^2)/(s^2 + 𝓂*σn^2)
        σ₂_η² = (∑η*σn^2)/(s^2 + 𝓂*σn^2)
        σ₃²   = (1/(2*sqrt(π)))*(1/(12*𝓂))
        σc_ζ² = (ζ̂c*σ₁²+σ₂_ζ²+σ₃²)
        σc_η² = (η̂c*σ₁²+σ₂_η²+σ₃²)
        #@show σc_ζ²,σc_η²

        #display(heatmap(Imgᵢ, color=:grays, title=string("Tracking Window #",ii)))

        d[ii] = Detection(Time,SVector(ζ̂c,η̂c),SVector(σc_ζ²,σc_η²),SVector(s),SVector(σ_s²),SVector(η̂c-W_m2/2, ζ̂c+W_m1/2, η̂c+W_m2/2, ζ̂c-W_m1/2))
    end

    return d, dsmum
end


function HierarchicalAgglomerative(Img;dmax::Int=50,T::Int=100, Time::Real=1)
    # Img     - m1xm2 array of intensities, an image
    # dmax    - maximum distance between two clusters to merge (condition to end loop)
    # T       - Threshold
    # mean of noise pdf assumed 0

    # Image size
    𝓂1 = size(Img)[1] #Vertical
    𝓂2 = size(Img)[2] #Horizontal

    # Tracking Window (3σ each diection)
    σ_s = sqrt(var(Img)) # Could try to calculate from known camera properties
    W_m1 = 10#2*ceil(3*σ_s) # Hard code window size???
    W_m2 = W_m1;

    # Detections using thresholding
    pᵢ = findall(x -> x>T,Img)

    # HA Clustering
    pk = [last.(Tuple.(pᵢ)) first.(Tuple.(pᵢ))]' # Cartesian index results of kmeans are horz x vert - need vert x horz for consistancy
    numclust = size(pk,2)
    dist = Euclidean()
    centroids = Float64.(pk)
    assign = [1:numclust;]
    Csizes = ones(size(assign))
    dmat = pairwise(dist, centroids, dims=2)
    dmat[dmat.==0] .= Inf # Avoid distance with itself
    while minimum(dmat) < dmax

        # Find Pair of clusters closest together
        indx = argmin(dmat)
        ii = indx[1]
        jj = indx[2]

        # Ensure ii<jj
        if ii > jj
            t = ii
            ii = jj
            jj = t
        end

        # Merge jj into ii - assign all points in jj to cluster ii
        assign[assign.==jj] .= assign[ii]

        # Calculate new centroid of cluster ii and set centroid of jj to Inf to ignore
        centroids[:,ii] = mean(pk[:,assign.==ii],dims=2)
        centroids[:,jj] .= Inf

        # Update Sizes of clusters ii and jj
        Csizes[ii] = Csizes[ii] + Csizes[jj]
        Csizes[jj] = 0

        # Update dmat with new distance
        pairwise!(dmat, dist, centroids, dims=2)
        dmat[dmat.==0] .= Inf
        dmat[isnan.(dmat)] .= Inf

        # Update numclust
        numclust -= 1
    end

    # Center of Each Cluster (in an array)
    k = numclust
    c1=centroids[1,:]
    c2=centroids[2,:]
    filter!(x->x != Inf,c1)
    filter!(x->x != Inf,c2)
    k_centers = [(c2[col],c1[col]) for col in 1:size(c1,1)]

    # Initialize Output array type
    d = Array{Detection}(undef,k)

    # Tracking Window of Each Cluster (vert,horz)
    Wv = []
    Wh = []
    for i in 1:size(k_centers,1)
        push!(Wv,clamp(round(k_centers[i][1]-W_m1/2),1,𝓂1-W_m1))
        push!(Wh,clamp(round(k_centers[i][2]-W_m2/2),1,𝓂2-W_m2))
    end

    # Loop over each tracking window
    for ii in 1:k
        # Create Subimage from Image
        Imgᵢ = Img[Int(Wv[ii]):Int(Wv[ii])+Int(W_m1),Int(Wh[ii]):Int(Wh[ii])+Int(W_m2)]

        # Extract Signal
        s = 0
        for i in CartesianIndices(Imgᵢ)
            s += Imgᵢ[i]
        end
        #@show s

        # Signal Uncertainty
        𝓂   = W_m1*W_m2
        σn   = σ_s # this σ_s is "signal diffraction width" which is different than the signal uncertainty
        σ_s² = (𝓂*σn^2)
        #@show σ_s²

        # Extract Centroid
        ∑ζI = 0
        ∑ηI = 0
        ∑I  = 0
        for i in CartesianIndices(Imgᵢ)
            ∑I  += Imgᵢ[i]
            ∑ζI += i[2]*Imgᵢ[i] # cartesian indicies are (horz x vert)
            ∑ηI += i[1]*Imgᵢ[i] # cartesian indicies are (horz x vert)
        end
        ζ̂c = clamp(∑ζI/∑I + Wv[ii],1,𝓂1) # vert
        η̂c = clamp(∑ηI/∑I + Wh[ii],1,𝓂2) # horz
        #@show ζ̂c,η̂c

        # Centroid Uncertainty
        ∑ζ  = 0
        ∑η  = 0
        for i in CartesianIndices(Imgᵢ)
            ∑ζ  += (i[2]-(W_m1/2)-.5)^2 # 0.5 accounts for index starting at 1
            ∑η  += (i[1]-(W_m2/2)-.5)^2 # 0.5 accounts for index starting at 1
        end
        σ₁²   = (𝓂*σn^2)/(s^2 + 𝓂*σn^2)
        σ₂_ζ² = (∑ζ*σn^2)/(s^2 + 𝓂*σn^2)
        σ₂_η² = (∑η*σn^2)/(s^2 + 𝓂*σn^2)
        σ₃²   = (1/(2*sqrt(π)))*(1/(12*𝓂))
        σc_ζ² = (ζ̂c*σ₁²+σ₂_ζ²+σ₃²)
        σc_η² = (η̂c*σ₁²+σ₂_η²+σ₃²)
        #@show σc_ζ²,σc_η²

        #display(heatmap(Imgᵢ, color=:grays, title=string("Tracking Window #",ii)))

        d[ii] = Detection(Time, SVector(ζ̂c,η̂c),SVector(σc_ζ²,σc_η²),SVector(s),SVector(σ_s²),SVector(η̂c-W_m2/2, ζ̂c+W_m1/2, η̂c+W_m2/2, ζ̂c-W_m1/2))
    end

    return d
end

export PixelThreshold, HierarchicalAgglomerative
