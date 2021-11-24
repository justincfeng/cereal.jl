using Test

using LinearAlgebra, Combinatorics, Statistics
using OrdinaryDiffEq, ForwardDiff, DiffResults

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

include("testfunctions.jl")
include("../cereal.jl")

    @time @testset "cereal tests:" begin 
            include("cereal_test.jl") 
    end

nothing
