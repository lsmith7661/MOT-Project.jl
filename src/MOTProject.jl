module MOTProject

# Remember to ]dev MOTProject to "add" it while developing
# ]activate MOTProject - add pkgs to module
# Base.compilecache(Base.PkgId(MOTProject)) - recompile without restarting julia

# Using pkgs for the whole module
using StaticArrays
using VideoIO, Plots, ImageView
using Distributions, LinearAlgebra
using DataStructures, ImageFiltering
using Clustering, Distances

# Include my support files with custome functions and objects and stuff
include("MOTobjects.jl")
include("ImageSim.jl")
include("Preprocess.jl")
include("BackgroundEstimation.jl")
include("ClusterMethods.jl")

# Export specific function (dont have to call MOTProject.myfun to use)
export
    Detection
    ImageSim
    DataPipe, SigmaClipping, PixelThreshold, HierarchicalAgglomerative
    Preprocess

#=
Make a DataPipe which contains an background estimation method, clustering method,
and threshold value

DataPipe{BackgroundMethod, ClusterMethod, Threshold}
=#

#function DataPipe(Img::AbstractMatrix, DataPipe::DataPipe)
#    σb, B, Ib = DataPipe.BackgroundMethod(Img)
    #ClusterResult = DataPipe.ClusterMethod(Ib,k,DataPipe.Threshold)
#end


end # end module
