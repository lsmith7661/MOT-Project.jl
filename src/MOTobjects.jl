struct telescope
    f::Float64 # Focal Length, meters
    D::Float64 # Aperture, meters
    σ::Float64 # Sensor Noise, pixels
    p::Float64 # pixel pitech, m
    λ::Float64 # wavelength, m
    𝓂::AbstractArray{Int,1} # Image Size (m1 x m2), pixels
end

struct target
    SNR::Float64 # Signal to Noise Ratio
    vel::Float64 # Target Velocity Magnitude, pixels/image
end
