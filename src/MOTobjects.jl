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

# --------
# A Detection is the result of extracting the point sources from an SSA image
struct Detection
    ID::Int
    Time::Float64
    Pd::Float64
    Centroid::SVector{2,Float64}
    CentroidUncertainty::SVector{2,Float64}
    Signal::SVector{1,Float64}
    SignalUncertainty::SVector{1,Float64}
    TrackingWindow::SVector{4,Float64}
end
# --------

#=struct DataPipe{T1, T2, T3}
    BackgroundMethod::T1
    ClusterMethod::T2
    Threshold::T3
end =#

# --------
# Object to hold information about a track (result of Kalman Filter)
struct Track
    Measurements::Array{Detection,1}
    State::Array{SMatrix{4,1,Float64},1}
    Covariance::Array{SMatrix{4,4,Float64},1}
    Residuals::Array{SMatrix{2,1,Float64},1}
end
# --------

# --------
# TrackHistory holds an array of Tracks to keep the history
struct TrackHistory
    Tracks::Array{Track,1}
end
# --------

# --------
# An AssociationMatrix is an TxS matrix of integers where each row, T
# is a unique track and each Column, S, is a time history containing
# the index of the Detection that is associated to that track at time, ts.
struct AssociationMatrix
    TOMHT::Array{Int,2}
end
# --------

# --------
# Clusters holds a single Track and all Detections within some measurement
# gating distance away (ex: mahalanobis distance)
struct Cluster
    Track::Track
    Detections::Array{Detection,1}
end
# --------

# --------
# HypothesisCost is the same size as AssociationMatrix and has a running value
# of the cost function for each Hypothesis. Col(1) = Cost(1), Col(2) = Cost(1) +
# Cost(2), Col(3) = Cost(2) + Cost(3) = Cost(1) + Cost(2) + Cost(3),... etc.
struct HypothesisCost
    CostMatrix::Array{Float64,2}
end
# --------

# --------
# Tuning object to hold all of the parameters I dont want to recompile to change
struct TuningParams
    ProcessNoise::Float64
    M::Int64
    ScanningWindow::Int64
    Pfa::Float64
    Inflator::Float64
    TrackingWindow::Int64
end
# --------
# --------
# Tuning object to hold all of the parameters I dont want to recompile to change
struct Frame
    Raw::Array{Float64,2}
    Background::Array{Float64,2}
    Detections::Vector{Detection}
    Tracks::Vector{Track}
    Associations::AssociationMatrix
    Cost::HypothesisCost
end
# --------

# --------
# Catalog - filled with all tracks that leave the scene and persist at the end of the processing
struct Catalog
    EndFrames::Vector{Int64} # Frame number where track ends
    Tracks::Vector{Track}    # Track
end
# --------

export telescope, target, Detection, Track, TrackHistory, AssociationMatrix, Cluster, HypothesisCost, TuningParams, Frame, Catalog
