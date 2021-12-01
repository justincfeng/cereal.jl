using Test, LinearAlgebra, Combinatorics, Statistics

include("../src/type.jl")
include("../src/cereal.jl")

    @time @testset "cereal tests:" begin 
            include("minkowski_test.jl")
            include("levicivita_test.jl")
            include("alg_test.jl")
    end

nothing
