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
end # end InitKF

# Init Kalman Filter with KF method
function KF(d1::Detection,d2::Detection)
    return InitKF(d1,d2)
end # end KF

function TimeUpdate(x::Track, tk::Float64, T::TuningParams)

    σw² = T.ProcessNoise

    tk₋ = x.Measurements[end].Time

    # State Transition Matix
    Φ(tk,tk₋) = [Matrix(I, 2, 2) (tk-tk₋)*Matrix(I, 2, 2); zeros(2,2) Matrix(I, 2, 2)]

    # Discrete-time Equivalent of B
    Γ(tk,tk₋) = [(exp(tk-tk₋)-1)*Matrix(I, 2, 2) ; (tk-tk₋)*Matrix(I, 2, 2)]

    # Process Noise
    function Qk(tk,tk₋,σw²)#,W)
        Δt = tk - tk₋
        #σw² = 0.1 #var(W)
        Q  = σw²*Matrix(I, 2, 2)
        return [(Q*Δt^3)/3 (Q*Δt^2)/2; (Q*Δt^2)/2 Q*Δt]
    end

    xk₋⁺ = x.State[end] # Adding [end] here probs breaks KF!
    Pk₋⁺ = x.Covariance[end]
    uk₋  = zeros(2,1) # No input
    #W    = ModTrack[tk][2][1][5]
    xk⁻  = Φ(tk,tk₋)*xk₋⁺ + Γ(tk,tk-1)*uk₋
    Pk⁻  = Φ(tk,tk₋)*Pk₋⁺*(Φ(tk,tk-1))' + Qk(tk,tk-1,σw²)#,W)

    return SMatrix{4,1}(xk⁻[:]), SMatrix{4,4}(Pk⁻[:])
end # End TimeUpdate

function MeasurementUpdate(xk⁻::SMatrix{4,1}, Pk⁻::SMatrix{4,4}, yk::SMatrix{2,1}, Rk::SMatrix{2,2})
    # Residuals, Kalman Gain
    H   = [Matrix(I, 2, 2) zeros(2,2)]
    rk⁻ = yk - H*xk⁻
    Kk  = Pk⁻*H'*inv(H*Pk⁻*H' + Rk) # Changed from Notes

    # Measurement Update
    xk⁺ = xk⁻ + Kk*rk⁻ # Changed from Notes
    rk⁺ = yk - H*xk⁺
    Pk⁺ = (Matrix(I, 4, 4) - Kk*H)*Pk⁻*(Matrix(I, 4,4)-Kk*H)' + Kk*Rk*Kk'

    return SMatrix{4,1}(xk⁺[:]),SMatrix{4,4}(Pk⁺[:]),SMatrix{2,1}(rk⁺[:])

end # End MeasurementUpdate

# Linear Kalman Filter for Image Space # CURRENTLY UNTESTED AFTER CHANGES
function KF!(T::TrackHistory,d::Detection,TP::TuningParams)

    # Time Update
    xk⁻, Pk⁻ = TimeUpdate(T.Tracks[end], d.Time, TP)

    # Read Next Observation
    yk = SMatrix{2,1}(d.Centroid[1], d.Centroid[2])
    R = [d.CentroidUncertainty[1]; d.CentroidUncertainty[2]].*Matrix(I, 2, 2)
    Rk = SMatrix{2,2}(R[:])

    # Measurement Update
    xk⁺,Pk⁺, rk⁺ = MeasurementUpdate(xk⁻, Pk⁻, yk, Rk)

    # Measurement History
    detections = [T.Tracks[end].Measurements;d]

    # Update Tracks
    push!(T.Tracks,Track(detections,SMatrix{4,1}(xk⁺[:]), SMatrix{4,4}(Pk⁺[:]), SMatrix{2,1}(rk⁺[:])))

    return T
end # end KF!

# Make a Dummy Detection (for missed detections)
function dummydetect(Time::Float64;Pd::Float64=1.)
    return Detection(0, Time, Pd, SVector(0.,0.),SVector(0.,0.),SVector(0.),SVector(0.),SVector(0., 0., 0., 0.))
end

# Birth Model
function BirthModel(d::Detection;n_dummys::Int64=0)
    # Starlink Birth Model - velocity
    tanFOV = (11/50)    # Camera half height is 11mm and 60mm focal length
    r = 6378.1+500      # Km radius of starlink
    v = sqrt(398600/r)  # km/s velocity (circular orbit)
    x = r*tanFOV        # half FOV in km
    pixperkm = 1080/x   # Pixels per KM
    vpix = v*pixperkm   # velocity in pixels/sec
    vpixframe = vpix/30 # velcity pixels per frame (30Hz framerate)

    state = vcat(d.Centroid[:], 0.0, 0.0)
    covar = Matrix(1.0I, 4, 4)
    covar[1,1] = d.CentroidUncertainty[1]
    covar[2,2] = d.CentroidUncertainty[2]
    covar[3,3] = sqrt(5) # This is hard coded for now
    covar[4,4] = sqrt(5) # This is hard coded for now
    r = [0. 0.]

    # Add a bunch of dummy detects for prior missed frames when birthing new track
    detectarray = [d]
    states = [SMatrix{4,1}(state[:])]
    covars = [SMatrix{4,4}(covar[:])]
    resids = [SMatrix{2,1}(r[:])]
    if n_dummys>0
        for ii = 1:n_dummys
            newdummy = dummydetect(d.Time-ii)
            detectarray = vcat(newdummy,detectarray)
            states = vcat([SMatrix{4,1}([0.;0.;0.;0.])],states)
            covars = vcat([SMatrix{4,4}([0. 0. 0. 0.; 0. 0. 0. 0.; 0. 0. 0. 0.; 0. 0. 0. 0.])],covars)
            resids = vcat([SMatrix{2,1}([0.;0.])],resids)
        end
    end

    return Track(detectarray,states,covars,resids)
end

# TrackCost without a PrevCost input calculates the cost of the current
# time step using all the available history in the track and assignments
function TrackCost(track::Track, assignments::Array{Int64,1})

    S = length(assignments)
    Cost = 0.
    λν = 0.2 # I dont know what this should be
    λc = 1080*1920*1e-5 # Hard code for now but should pass this in
    λex = λν + λc
    for ii in 1:S
        #Pd = clamp(track.Measurements[end].Pd,0.001,.999)
        Pd = .999
        Σ = track.Covariance[ii][1:2,1:2]
        r = track.Residuals[ii]
        if assignments[ii] == 0
            δ = 0
            Term1 = (δ-1)*log(1-Pd) +100
            Term2 = 0
            Term3 = 0
        else
            δ = 1
            Term1 = (δ-1)*log(1-Pd)
            Term2 = -δ*log(Pd/(λex*sqrt(norm(2*π*Σ))))
            Term3 = δ*(1/2*r'*inv(Σ)*r)[end] # Do this in the IF becuase zeros'*inv(zeros)*zeros fails for dummy detections
        end
        c = Term1 + Term2 + Term3
        Cost += c
    end
    return Cost
end # End TrackCost

# TrackCost with a PrevCost input calculates updates the cost assuming the
# previous cost is a running sum of all previous cost
function TrackCost(track::Track, assignments::Array{Int64,1}, PrevCost::Float64)

    λν = 0. # I dont know what this should be but since im cheating on gamma anyways its probably fine for now
    λc = 1080*1920*1e-5 # Hard code for now but should pass this in
    λex = λν + λc
    #Pd = clamp(track.Measurements[end].Pd,0.001,.999)
    Pd = .999
    Σ = track.Covariance[end][1:2,1:2]
    r = track.Residuals[end]
    if assignments[end] == 0
        δ = 0
        Term1 = (δ-1)*log(1-Pd) + 100
        Term2 = 0
        Term3 = 0
    else
        δ = 1
        Term1 = (δ-1)*log(1-Pd)
        Term2 = -δ*log(Pd/(λex*sqrt(norm(2*π*Σ))))
        Term3 = δ*(1/2*r'*inv(Σ)*r)[end] # Do this in the IF becuase zeros'*inv(zeros)*zeros fails for dummy detections
    end
    #show δ, Term1, Term2, Term3, assignments
    c = Term1 + Term2 + Term3
    Cost = c + PrevCost

    return Cost
end # End TrackCost

# Init TOMHT on first frame with only Detections
function TOMHT(Detections::Array{Detection,1})
    # Number of Detections
    nd = length(Detections)

    # Create a dummy track for each detection and assign it to an AssociationMatrix
    tracks = Array{Track,1}(undef,nd)
    AssMat = Array{Int,2}(undef,nd,1)
    cost = zeros(Float64,size(AssMat))
    for ii in 1:nd
        tracks[ii] = BirthModel(Detections[ii])
        AssMat[ii] = Detections[ii].ID
        cost[ii]   = TrackCost(tracks[ii], [AssMat[ii]], 0.0)
    end

    return tracks, AssociationMatrix(AssMat), HypothesisCost(cost), Catalog(Array{Int64,1}(undef,0),Array{Track,1}(undef,0))
end # end TOMHT

# Track Oriented Multiple Hypothesis Tracker (TOMHT)
function TOMHT!(Detections::Array{Detection,1},Tracks::Array{Track,1},AssMat::AssociationMatrix,Cost::HypothesisCost, catalog::Catalog, T::TuningParams, FrameNum::Int64)

    # Check AssMat and Tracks are the same size
    assmat  = copy(AssMat.TOMHT)
    CostMat = copy(Cost.CostMatrix)
    nTracks = length(Tracks)
    nAssMat = size(assmat,1)
    ndetect = length(Detections)

    # Add a new empty column to AssMat and CostMat that needs to be filled
    TrackHypoth = copy(Tracks) # Dont want to return a reference - or do i?
    newcol      = Array{Int,2}(undef,nAssMat,1)
    newassmat   = hcat(assmat,newcol)
    newcostmat  = hcat(CostMat,newcol)
    newcat      = copy(catalog.Tracks)
    frameends = copy(catalog.EndFrames)

    # Create Hypotheses for each track-detection pair in each cluster
    clusters = Array{Cluster,1}(undef,nTracks)
    rowoffset = 0 # Row Offset Value (add to this when a new row (hypothesis) is inserted for a track)
    for ii in 1:nTracks

        # Time Update to get predicted state at current frame
        tk = Detections[1].Time
        xk⁻, Pk⁻ = TimeUpdate(Tracks[ii], tk, T)

        # Measurement Gating - make clusters of detections near each track (in measurment space)
        H = [Matrix(I, 2, 2) zeros(2,2)]
        P = H*Pk⁻*H'
        MD = zeros(Float64,size(Detections))
        for ii in 1:length(Detections)
            d = Detections[ii]
            R = [d.CentroidUncertainty[1]; d.CentroidUncertainty[2]].*Matrix(I, 2, 2)
            S = P + R
            rm = reshape(d.Centroid,2,1) - H*xk⁻
            MD[ii] = (rm'*inv(S)*rm)[end] # Mahalanobis Distance
        end
        boolindx = vec(MD.<5^2) # 2 sigma mahalanobis distance
        clusters[ii] = Cluster(Tracks[ii], Detections[boolindx])

        # Missed Detection Hypothesis - will always be the first hypothesis
        dummy = dummydetect(Detections[1].Time)#,Pd=Tracks[ii].Measurements[end].Pd)
        TrackHypoth[ii+rowoffset] = Track(vcat(Tracks[ii].Measurements,dummy),vcat(Tracks[ii].State,[SMatrix{4,1}(xk⁻[:])]), vcat(Tracks[ii].Covariance,[SMatrix{4,4}(Pk⁻[:])]), vcat(Tracks[ii].Residuals,[SMatrix{2,1}([0.;0.])]))
        newassmat[ii+rowoffset,end] = dummy.ID # Missed Detection
        newcostmat[ii+rowoffset,end] = TrackCost(TrackHypoth[ii+rowoffset], newassmat[ii+rowoffset,:], CostMat[ii,end])

        # Assign each detection in cluster as possible association hypothesis
        if !isempty(clusters[ii].Detections)
            for d in clusters[ii].Detections
                # Assignment
                newrow    = copy(newassmat[ii+rowoffset,:]) # Make a copy of the track (row)
                newrow[end] = d.ID                          # Assign detection to the track
                frontass  = newassmat[1:ii+rowoffset,:]
                backass   = newassmat[ii+rowoffset+1:end,:]
                newassmat = vcat(frontass, newrow', backass)

                # Track Hypothesis Creation
                xk⁻ = TrackHypoth[ii+rowoffset].State[end] # Remember current state has been time updated
                Pk⁻ = TrackHypoth[ii+rowoffset].Covariance[end]
                yk  = SMatrix{2,1}(d.Centroid[1], d.Centroid[2])
                R   = [d.CentroidUncertainty[1]; d.CentroidUncertainty[2]].*Matrix(I, 2, 2)
                Rk  = SMatrix{2,2}(R[:])
                xk⁺,Pk⁺,rk⁺ = MeasurementUpdate(xk⁻, Pk⁻, yk, Rk)
                newTrack    = Track(vcat(Tracks[ii].Measurements[1:end],d),vcat(Tracks[ii].State[1:end],[SMatrix{4,1}(xk⁺[:])]), vcat(Tracks[ii].Covariance[1:end],[SMatrix{4,4}(Pk⁺[:])]), vcat(Tracks[ii].Residuals[1:end],[SMatrix{2,1}(rk⁺)]))
                fronttrack  = TrackHypoth[1:ii+rowoffset]
                backtrack   = TrackHypoth[ii+rowoffset+1:end]
                TrackHypoth = vcat(fronttrack, newTrack, backtrack)

                # Cost Calculation
                newcostrow = copy(newcostmat[ii+rowoffset,:])
                newcostrow[end] = TrackCost(TrackHypoth[ii+rowoffset], newrow, CostMat[ii,end])
                frontcost  = newcostmat[1:ii+rowoffset,:]
                backcost   = newcostmat[ii+rowoffset+1:end,:]
                newcostmat = vcat(frontcost, newcostrow', backcost)

                # Increment Counter
                rowoffset += 1
            end
        end
    end

    # Lastly, detections could be a new track
    n_prevframe = T.ScanningWindow
    for ii in 1:ndetect
        TrackHypoth = vcat(TrackHypoth, BirthModel(Detections[ii],n_dummys=n_prevframe))
        newrow = zeros(Int64,size(newassmat[1,:]))
        newrow[end] = Detections[ii].ID
        newassmat = vcat(newassmat,newrow')
        newcostrow = zeros(Float64,size(newcostmat[1,:]))
        TrackMultiplyer = max(1.,FrameNum-T.ScanningWindow)
        newcostrow[1] = TrackMultiplyer*TrackCost(TrackHypoth[end], [newrow[1]], 0.0)
        for jj in 2:length(newcostrow)-1
            newcostrow[jj] = newcostrow[jj-1] + newcostrow[1]
        end
        newcostrow[end] = TrackCost(TrackHypoth[end], [newrow[end]], newcostrow[end-1])
        newcostmat = vcat(newcostmat,newcostrow')
    end

    # Only save info for 5 frames
    if size(newassmat,2) > T.ScanningWindow
        # Cut First Column
        newassmat  = newassmat[:,2:end]
        newcostmat = newcostmat[:,2:end]
        # A rows of all dummy tracks means track is gone (ex: [1 0 0 0 0 0] -> [0 0 0 0 0] -> Write to catalog)
        rmidx       = findall(x -> x==0 , sum(newassmat,dims=2)[:])
        newassmat   = newassmat[setdiff(1:end, rmidx), :]
        newcostmat  = newcostmat[setdiff(1:end, rmidx), :]
        # Put Cut Tracks into the catalog
        newcat      = vcat(catalog.Tracks, TrackHypoth[rmidx])
        frames      = Vector{Int64}(undef,length(rmidx))
        frames     .= FrameNum
        frameends   = vcat(catalog.EndFrames,frames)
        TrackHypoth = TrackHypoth[setdiff(1:end, rmidx), :][:]
    end

return TrackHypoth, AssociationMatrix(newassmat), HypothesisCost(newcostmat), Catalog(frameends,newcat)

end # end TOMHT!

export KF, KF!, TOMHT, TOMHT!
