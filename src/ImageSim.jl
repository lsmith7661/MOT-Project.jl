function ImageSim(target,scope,n_img)
    # target  - 1xn array of SNRs for targets
    #            - SNR
    #            - velocity magnitude (pixels/image)
    # scope   - Telescope and Camera parameters
    #            - f (m)
    #            - D (m)
    #            - σ (pixels) sensor noise
    #            - pixel pitch (m)
    #            - λ (m)
    #            - 𝓂 (m1,m2)
    # n_img   - number of images to create in sequence

    # Define output dictionary
    ImgSeq = Dict()

    # Extract Target Parameters
    SNR = target.SNR
    vel = target.vel

    # Extract Telescope Parameters
    f  = scope.f
    D  = scope.D
    σ  = scope.σ
    p  = scope.p
    λ  = scope.λ
    m1 = scope.𝓂[1]  # Image size x
    m2 = scope.𝓂[2]  # Image size y
    σ_s = scope.σ_s

    # Derived Values
    n    = length(SNR)   # number of objects in image
    𝓂   = m1*m2         # total number of pixels
    μx   = rand(1:m1,n)  # random centroid of each object, x
    μy   = rand(1:m2,n)  # random centroid of each object, y
    W_m1 = 2*ceil(3*σ_s) # reasonable size of the tracking window
    W_m2 = W_m1          # reasonable size of the tracking window
    W𝓂  = W_m1*W_m2     # tracking window pixels
    s    = [snr*σ*√W𝓂 for snr in SNR]
    # TODO: Add loop for each μx and μy for case with multiple objects in single image
    vdir = ([m1/2 m2/2] - [μx μy])/norm([m1/2 m2/2] - [μx μy]) # Velocity direction towards center of frame
    V    = round.(vel.*vdir)

    # Noise - Normal Distribution
    n_pdf = Normal(0, σ)

    # Loop over Image Sequence
    for seq_i in 1:n_img

        # Image
        Img = Float64[] # Start with an empty array

        # Add Sensor Noise to Image
        for pixel in 1:𝓂
            push!(Img,rand(n_pdf))
        end

        # Reshape the Image into the correct size
        Img = reshape(Img, (m1,m2))

        # Create Signal Distributions - Multivariate Normal Distribution
        s_pdf = []
        if μx[1] > 0 || μx[1] <= m1
            if μy[1] > 0 || μy[1] <= m2
                for i in 1:n
                    push!(s_pdf, MvNormal([μx[i];μy[i]], σ_s^2*Matrix(I, 2, 2)))
                end
            end
        end

        # Photons - Estimate where each photon lands
        sᵢ = []
        if μx[1] > 0 || μx[1] <= m1
            if μy[1] > 0 || μy[1] <= m2
                for i in 1:n
                    push!(sᵢ,Int.(round.(clamp.(rand(s_pdf[i],Int(s[i])),1,m1))))
                    [Img[sᵢ[i][1,j],sᵢ[i][2,j]]+=1 for j in 1:size(sᵢ[i],2)]
                end
            end
        end

        #display(heatmap(Img, color = :greys, title = string("Velocity = ",V)))

        # Object truth information
        Objects = Dict(i => (μx[i],μy[i],s[i],s_pdf[i], V) for i in 1:n)

        # Update Object Location
        [μx[ii] += V[1] for ii in 1:n]
        [μy[ii] += V[2] for ii in 1:n]

        # Dictionary to return all images and truth locations
        push!(ImgSeq, seq_i => (Img , Objects))
    end

    return sort(OrderedDict(ImgSeq))
end
