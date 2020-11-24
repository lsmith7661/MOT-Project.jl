module MOTProject

# Remember to ]dev MOTProject to "add" it while developing
# ]activate MOTProject - add pkgs to module
# Base.compilecache(Base.PkgId(MOTProject)) - recompile without restarting julia

# Using pkgs for the whole module
using StaticArrays
using VideoIO, ImageView
using Plots: Gray, @layout
using Distributions, LinearAlgebra
using DataStructures, ImageFiltering, SpecialFunctions
using Clustering
using Distances

# Include my support files with custome functions and objects and stuff
include("MOTobjects.jl")
include("ImageSim.jl")
include("Preprocess.jl")
include("BackgroundEstimation.jl")
include("ClusterMethods.jl")
include("Filters.jl")
include("Pruning.jl")

# Export specific function (dont have to call MOTProject.myfun to use)
export
    Detection,
    ImageSim,
    SigmaClipping, PixelThreshold, HierarchicalAgglomerative,
    Preprocess,
    KF, KF!, TOMHT, TOMHT!,
    MBest

end # end module
