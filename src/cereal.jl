#-----------------------------------------------------------------------
module cereal      # cereal module
#-----------------------------------------------------------------------

using LinearAlgebra
using Combinatorics
using Statistics

include("type.jl")
include("minkowski.jl")
include("levicivita.jl")
include("multi.jl")
include("eval.jl")

include("alg/RTC21.jl")
include("alg/CFM10.jl")
include("alg/FHC21.jl")

#-----------------------------------------------------------------------
"""
# Locator function

    locatorfunc( N::Int , Method::String )

This function selects the locator function to use. The first argument
`N` determines the minimum number of emission points that the locator
function will use. `N=4` or `N≥5` is recommended. The second argument
`Method` selects the formula or algorithm used to compute the
intersection point ``X_c``. The available methods are `"FHC21"`,
`"CFM10"`, and for `N≥5`, `"RTC21"`.

By default, the locator function assumes `N=5` and `Method="RTC21"`:

    julia> locatorfunc() == locatorfunc(5,"RTC21")
    true

"""
function locatorfunc( N::Int=5 , Method::String="RTC21" )
    if N < 4 
        print("No methods available.")
        return X -> zeros(Float64,4)
    elseif  N == 4
        if Method == "RTC21"
            print("RTC21 requires 5 emission points.")
            return X -> zeros(Float64,4)
        elseif Method == "CFM10"
            return locator4CFM10
        elseif Method == "FHC21"
            return locator4FHC21
        else 
            print("Unrecognized method.")
            return X -> zeros(Float64,4)
        end
    elseif  N == 5 && Method == "RTC21"
        return locator5RTC21
    elseif  N > 5 && Method == "RTC21"
        return (X,q=1e-14) -> mlocator(X,locator5RTC21,5,false,q)
    elseif  N >= 5 && Method != "RTC21"
        if  Method == "CFM10"
            return (X,q=1e-14) -> mlocator(X,locator4CFM10,4,true,q)
        elseif  Method == "FHC21"
            return (X,q=1e-14) -> mlocator(X,locator4FHC21,4,true,q)
        else
            print("Unrecognized method.")
            return (X,q=1e-14) -> zeros(Float64,4)
        end
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # End scope of module cereal
#-----------------------------------------------------------------------