# Support files containing all functions for Kalman Filtering and MHTs

function InitKF(TrackSeq::Array{Detection,1})
    # TrackSeq - Array of Detections

    # Z matrix
    Δt = 1                       # 1 Frame
    TL = zeros(2,2)              # Top Left
    TR = Matrix(I, 2, 2)         # Top Right
    BL = (-1/Δt)*Matrix(I, 2, 2) # Bottom Left
    BR = (1/Δt)*Matrix(I, 2, 2)  # Bottom Right
    Z₀ = [TL TR; BL BR]          # Make Z Matrix

    # Centroid Measurments
    y_prev = [TrackSeq[1].Centroid[1]; TrackSeq[1].Centroid[2]]
    y₀     = [TrackSeq[2].Centroid[1]; TrackSeq[2].Centroid[2]]

    # Initialized State x₀
    x₀ = Z₀*[y_prev ; y₀]

    # Centroid Uncertainties
    R_prev = [TrackSeq[1].CentroidUncertainty[1]; TrackSeq[1].CentroidUncertainty[2]].*Matrix(I, 2, 2)
    R₀     = [TrackSeq[2].CentroidUncertainty[1]; TrackSeq[2].CentroidUncertainty[2]].*Matrix(I, 2, 2)
    R      = [R_prev TL; TL R₀]

    # Initialize Covariance P₀
    P₀ = Z₀*R*Z₀'

    return x₀, P₀
end
