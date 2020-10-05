module MOTProject

# Remember to ]dev MOTProject to "add" it while developing
# ]activate MOTProject - add pkgs to module

# Using pkgs for the whole module
using Distributions
using LinearAlgebra
using DataStructures
using ImageFiltering

# Include my support files with custome functions and objects and stuff
include("MOTobjects.jl")
include("ImageSim.jl")
include("BackgroundEstimation.jl")

# Export specific function (dont have to call MOTProject.myfun to use)
export
    ImageSim
    SigmaClipping

end
