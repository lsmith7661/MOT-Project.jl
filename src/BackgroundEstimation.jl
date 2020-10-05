# Background Estiamtion

function SigmaClipping(img,Wpsf::Int=10,n::Int=5)
    # Image
    # Estimated psf size for smoothing window
    # Number of interations for clipping

    Im = img
    B = []
    #display(heatmap(Im, title="Raw Image"))

    # Noise Estimation
    σn = sqrt(var(img))

    # Smoothing Filters
    Wl = ones(20*Wpsf,20*Wpsf)
    Ws = ones(10*Wpsf,10*Wpsf);

    for i in 1:n
        # Smooth Image - background estimation
        if isodd(i)
            B = imfilter(Im,reflect(centered(Wl)),"replicate")./sum(Wl[:])
        else
            B = imfilter(Im,reflect(centered(Ws)),"replicate")./sum(Ws[:])
        end
        #display(heatmap(B.+(2*σn)))

        # Image Clipping
        Im = min.(Im,B.+(2*σn))
        #display(heatmap(Im, title=string("Sigma Clipping - Iteration",i)))
    end

    B = imfilter(Im,reflect(centered(Wl)),"replicate")./sum(Wl[:])
    #display(heatmap(B, title = "Background Estimation Image"))

    # Background Estimation Error
    E = zeros(size(img))
    kernal_Wpsf = ones(Wpsf,Wpsf)
    Bigkernal_Wpsf = ones(5*Wpsf,5*Wpsf)
    localmean = imfilter((img-B),reflect(centered(kernal_Wpsf)),"replicate")./sum(kernal_Wpsf[:])
    Biglocalmean = imfilter((img-B),reflect(centered(Bigkernal_Wpsf)),"replicate")./sum(Bigkernal_Wpsf[:])
    for pix in cat(findall(x-> x>σn,localmean), findall(x-> x>10*σn,Biglocalmean), dims=1)
        E[pix] += 1
        E[pix] = min(E[pix],1)
    end

    # Noise Image
    N = 1 .- E

    # Local Mean Residuals
    NI = imfilter((N.*(img-B)),reflect(centered(kernal_Wpsf)),"replicate")
    Nf = imfilter((N),reflect(centered(kernal_Wpsf)),"replicate")
    μr = NI./(Nf.+1e-9)

    Ib = img-B
    σb = sqrt(var(μr))

    return σb, B, Ib
end
