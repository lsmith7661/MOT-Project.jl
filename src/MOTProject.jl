module MOTProject

# Remember to ]dev MOTProject to "add" it while developing

# Using pkgs for the whole module
using Distributions
using LinearAlgebra
using DataStructures

# Include my support files with custome functions and objects and stuff
include("MOTobjects.jl")
include("ImageSim.jl")

# Export specific function (dont have to call MOTProject.myfun to use)
export ImageSim

end
