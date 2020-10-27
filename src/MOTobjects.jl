struct telescope
    f::Float64 # Focal Length, meters
    D::Float64 # Aperture, meters
    σ::Float64 # Sensor Noise, pixels
    p::Float64 # pixel pitech, m
    λ::Float64 # wavelength, m
    𝓂::AbstractArray{Int,1} # Image Size (m1 x m2), pixels
    σ_s::Float64 # Std Dev of signal derived from telescope parameters
    telescope(f,D,σ,p,λ,𝓂) = new(f,D,σ,p,λ,𝓂,0.45*λ*f/D/p)
end

struct target
    SNR::Float64 # Signal to Noise Ratio
    vel::Float64 # Target Velocity Magnitude, pixels/image
end

struct Detection
    Time::Float64
    Centroid::SVector{2,Float64}
    CentroidUncertainty::SVector{2,Float64}
    Signal::SVector{1,Float64}
    SignalUncertainty::SVector{1,Float64}
    TrackingWindow::SVector{4,Float64}
end

struct DataPipe{T1, T2, T3}
    BackgroundMethod::T1
    ClusterMethod::T2
    Threshold::T3
end

struct Track
    Detections::Array{Array{Detection}}
    State::Array{SVector{Float64,1}}
    Covariance::Array{SMatrix{2,2,Float64}}
    Residuals::Array{SVector{Float64,1}}
end

export telescope, target, Detection, DataPipe, Track
