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


#= Linear Kalman Filter for Image Space
function KF(TrackSeq::Array{Detection,1})

    #Initiate KF
    k = 1
    x₀, P₀ = InitKF(TrackSeq)
    r₀ = zeros(2,1)
    KF = OrderedDict(k => [x₀, P₀, r₀])

    # State Transition Matix
    Φ(tk,tk₋) = [Matrix(I, 2, 2) (tk-tk₋)*Matrix(I, 2, 2); zeros(2,2) Matrix(I, 2, 2)]

    # Process Noise
    function Qk(tk,tk₋,W)
        Δt = tk - tk₋
        σw² = 0.01#var(W)
        Q  = σw²*Matrix(I, 2, 2)
        return [(Q*Δt^3)/3 (Q*Δt^2)/2; (Q*Δt^2)/2 Q*Δt]
    end

    # Discrete-time Equivalent of B
    Γ(tk,tk₋) = [(exp(tk-tk₋)-1)*Matrix(I, 2, 2) ; (tk-tk₋)*Matrix(I, 2, 2)]

    # Remove First Track
    ModTrack = filter(x -> x.first != 1, TrackSeq)

    # Iterate for each next image
    for (key,value) in ModTrack

        # Read Next Observation
        tk = key
        yk = [ModTrack[tk][2][1][3][1]; ModTrack[tk][2][1][3][2]]
        Rk = [ModTrack[tk][2][1][4][1]; ModTrack[tk][2][1][4][2]].*Matrix(I, 2, 2)

        # Time Update
        xk₋⁺ = KF[key-1][1]
        Pk₋⁺ = KF[key-1][2]
        uk₋  = zeros(2,1) # No input
        W    = ModTrack[tk][2][1][5]
        xk⁻  = Φ(tk,tk-1)*xk₋⁺ + Γ(tk,tk-1)*uk₋
        Pk⁻  = Φ(tk,tk-1)*Pk₋⁺*(Φ(tk,tk-1))' + Qk(tk,tk-1,W)

        # Residuals, Kalman Gain
        H   = [Matrix(I, 2, 2) zeros(2,2)]
        rk⁻ = yk - H*xk⁻
        Kk  = Pk⁻*H'*inv(H*Pk⁻*H' + Rk) # Changed from Notes

        # Measurement Update
        xk⁺ = xk⁻ + Kk*rk⁻ # Changed from Notes
        rk⁺ = yk - H*xk⁺
        Pk⁺ = (Matrix(I, 4, 4) - Kk*H)*Pk⁻*(Matrix(I, 4,4)-Kk*H)' + Kk*Rk*Kk'

        # Update KF Dictionary
        push!(KF, key => [xk⁺,Pk⁺,rk⁺])

    end

    return KF

end

export KF
=#
