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
    n::Int
    m::Int
    sketchsize::Int
    sketch::AbstractString   #Type of sketch
    tol::Float64
    restol::Float64
    maxiter::Int
    skip_error_calculation::Int   #How many iterations should be skipped before computing error
    max_time::Float64
    printiters::Bool 
    exacterror::Bool
end


type Prob
    A::Array{Float64}
    b::Array{Float64}
    xsol::Array{Float64}
    name::AbstractString
end

type Method
    flopsperiter::Int
    name::AbstractString
    stepmethod::Function
    bootmethod::Function
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

#Including iterative methods for calculating pseudoinverse

#Including utilities, plotting, data analysis
include("plot_outputs_Plots.jl")


#end