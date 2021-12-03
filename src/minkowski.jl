#-----------------------------------------------------------------------
#   Minkowski product, norm and matrix
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
Minkowski product

    η( V1::RealVec , V2::RealVec )

The function `η` takes two vectors ``V_1`` and ``V_2`` of equal length, and
computes their Minkowski product ``η(V_1,V_2)``, with the assumption
that the first element of the vectors corresponds to the time component.

"""
function η( V1::RealVec , V2::RealVec )
    l1 = length(V1)
    l2 = length(V2)

    met = -V1[1]*V2[1]

    if l1==l2
        for i=2:l1
            met += V1[i]*V2[i]
        end
        return met
    else
        print("Vectors are of a different dimension.")
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
Minkowski norm

    mnorm( V ) 

The function `mnorm` computes the Minkowski norm, which is mathematically
equivalent to the absolute value of the square root of the Minkowski
product ``√|η(V_1,V_2)|``. However, this function computes the result
according to the formula used in the hypot function.

"""
function mnorm( V )
    l = length(V)

    t = abs(V[1])
    x = norm(V[2:l])

    if t > x
        return abs(t - ( x/t )*( x / ( 1 + √( 1 - (x/t)^2 ) ) ))
    elseif x > t
        return abs(x - ( t/x )*( t / ( 1 + √( 1 - (t/x)^2 ) ) ))
    else
        return zero( typeof(V[1]) )
    end
end      #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
Minkowski components

    ημν( tpfl::DataType , dim::Int )

The function `ημν` constructs the components of the Minkowski metric. The
first argument `tpfl` specifies the floating point datatype (typically
Float64 or Double64 if one uses the DoubleFloats package) and the second
argument `dim` specifies the dimension. The default values are 
`tpfl=Float64` and `dim=4`:

    julia> ημν() == ημν( Float64 , 4 )
        true

"""
function ημν( tpfl::DataType=Float64 , dim::Int=4 )   # Minkowski metric
    gη = Matrix(one(tpfl)*(I(dim)))

    gη[1,1] = -gη[1,1]

    return gη
end     #---------------------------------------------------------------