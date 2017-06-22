#module BigRidge  # Think later about organizing this into a module
#using StatsBase
using Match
#using Plots

# The types exported
#export MyOptions, Prob, Method, Output
# The exported functions
#export solveRidge
#export plot_outputs_Plots

type MyOptions
    sketchsize::Int
    sketch::AbstractString   #Type of sketch
    tol::Float64
    restol::Float64
    maxiter::Int
    skip_error_calculation::Int   #How many iterations should be skipped before computing error
    max_time::Float64
    printiters::Bool 
    exacterror::Bool
    AUX::Array{Float64} # auxiliary options used for random stuff
end


type Prob
    A::Array{Float64}
    b::Array{Float64}
    xsol::Array{Float64}
    n::Int
    lambda::Float64
    name::AbstractString
end

type SketchMethod
    flopsperiter::Int
    name::AbstractString
    stepmethod::Function
    bootmethod::Function
    SA::Array{Float64}
    SAS::Array{Float64}
    Sb::Array{Float64}
    ind::Array{Int64}
    sigs::Array{Int64}
end

type CGMethod
    flopsperiter::Int
    name::AbstractString
    stepmethod::Function
    bootmethod::Function
    p::Array{Float64}
    r::Array{Float64}
    rr::Float64
end



type Output
    iterations::Int
    flopsperiter::Int
    times::Array{Float64}
    errors::Array{Float64}
    residuals::Array{Float64}
    name::AbstractString
    fail::AbstractString
end
#Including method wrappers 
include("solveRidge.jl")
include("boot_method.jl")
#Including test and problem generating functions
include("dataLoad.jl")
include("defineProb.jl")
include("get_randProb.jl") # generates a random Gaussian matrix
#Including iterative methods
include("grad.jl")
include("CD.jl")
include("CG.jl")
include("Hadamard.jl")
include("rademacher.jl")
include("countmin.jl")
#Including utilities, plotting, data analysis
include("plot_outputs_Plots.jl")
#include("call_all_methods.jl")


#end
