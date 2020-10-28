# Support files containing all functions for Kalman Filtering and MHTs

function InitKF(d1::Detection,d2::Detection)
    # d1 - first detection of object
    # d2 - second detection of object (assumed to be associated)

    # Z matrix
    Δt = 1                       # 1 Frame
    TL = zeros(2,2)              # Top Left
    TR = Matrix(I, 2, 2)         # Top Right
    BL = (-1/Δt)*Matrix(I, 2, 2) # Bottom Left
    BR = (1/Δt)*Matrix(I, 2, 2)  # Bottom Right
    Z₀ = [TL TR; BL BR]          # Make Z Matrix

    # Centroid Measurments
    y_prev = [d1.Centroid[1]; d1.Centroid[2]]
    y₀     = [d2.Centroid[1]; d2.Centroid[2]]

    # Initialized State x₀
    x₀ = Z₀*[y_prev ; y₀]

    # Centroid Uncertainties
    R_prev = [d1.CentroidUncertainty[1]; d1.CentroidUncertainty[2]].*Matrix(I, 2, 2)
    R₀     = [d2.CentroidUncertainty[1]; d2.CentroidUncertainty[2]].*Matrix(I, 2, 2)
    R      = [R_prev TL; TL R₀]

    # Initialize Covariance P₀
    P₀ = Z₀*R*Z₀'

    # Initialize residuals r₀
    r₀ = zeros(2,1)

    # Create Track Object
    T = Track([d1,d2],SMatrix{4,1}(x₀[:]), SMatrix{4,4}(P₀[:]), SMatrix{2,1}(r₀[:]))

    return TrackHistory([T])
end

# Init Kalman Filter with KF method
function KF(d1::Detection,d2::Detection)
    return InitKF(d1,d2)
end

# Linear Kalman Filter for Image Space
function KF!(T::TrackHistory,d::Detection)

    # Time
    tk  = d.Time
    tk₋ = T.Tracks[end].Measurements[end].Time

    # State Transition Matix
    Φ(tk,tk₋) = [Matrix(I, 2, 2) (tk-tk₋)*Matrix(I, 2, 2); zeros(2,2) Matrix(I, 2, 2)]

    # Process Noise
    function Qk(tk,tk₋)#,W)
        Δt = tk - tk₋
        σw² = 0.01#var(W)
        Q  = σw²*Matrix(I, 2, 2)
        return [(Q*Δt^3)/3 (Q*Δt^2)/2; (Q*Δt^2)/2 Q*Δt]
    end

    # Discrete-time Equivalent of B
    Γ(tk,tk₋) = [(exp(tk-tk₋)-1)*Matrix(I, 2, 2) ; (tk-tk₋)*Matrix(I, 2, 2)]

    # Read Next Observation
    yk = [d.Centroid[1]; d.Centroid[2]]
    Rk = [d.CentroidUncertainty[1]; d.CentroidUncertainty[2]].*Matrix(I, 2, 2)

    # Time Update
    xk₋⁺ = T.Tracks[end].State
    Pk₋⁺ = T.Tracks[end].Covariance
    uk₋  = zeros(2,1) # No input
    #W    = ModTrack[tk][2][1][5]
    xk⁻  = Φ(tk,tk₋)*xk₋⁺ + Γ(tk,tk-1)*uk₋
    Pk⁻  = Φ(tk,tk₋)*Pk₋⁺*(Φ(tk,tk-1))' + Qk(tk,tk-1)#,W)

    # Residuals, Kalman Gain
    H   = [Matrix(I, 2, 2) zeros(2,2)]
    rk⁻ = yk - H*xk⁻
    Kk  = Pk⁻*H'*inv(H*Pk⁻*H' + Rk) # Changed from Notes

    # Measurement Update
    xk⁺ = xk⁻ + Kk*rk⁻ # Changed from Notes
    rk⁺ = yk - H*xk⁺
    Pk⁺ = (Matrix(I, 4, 4) - Kk*H)*Pk⁻*(Matrix(I, 4,4)-Kk*H)' + Kk*Rk*Kk'

    # Measurement History
    detections = [T.Tracks[end].Measurements;d]

    # Update Tracks
    push!(T.Tracks,Track(detections,SMatrix{4,1}(xk⁺[:]), SMatrix{4,4}(Pk⁺[:]), SMatrix{2,1}(rk⁺[:])))

    return T
end

export KF, KF!
