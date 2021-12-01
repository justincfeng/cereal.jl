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

#-----------------------------------------------------------------------
module ceval    # Evaluation module for cereal
#-----------------------------------------------------------------------

using LinearAlgebra

include("type.jl")
include("minkowski.jl")

#-----------------------------------------------------------------------
"""
# Vector generator

    vgenerator( tpfl::DataType )

This function generates a random 3-vector of unit length. Returns a 
three component vector.

"""
function vgenerator( tpfl::DataType=Float64 ) 
    b = false
    v = zeros(tpfl,3)

    while b == false    # Generate random points in a unit box
        vx = 2*(rand(tpfl)-1/2)
        vy = 2*(rand(tpfl)-1/2)
        vz = 2*(rand(tpfl)-1/2)
        vl = norm([vx;vy;vz])
        if vl>1         # Cut out points outside unit sphere
            b = false
        else
            b = true
        end
        v = tpfl[ vx ; vy ; vz ]/vl     # Normalization
    end
    return v
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Null vector generator

    nullgen( tpfl::DataType )

This function generates a random past directed null vector. Returns a
four component vector.

"""
function nullgen( tpfl::DataType=Float64 )   
    v = vgenerator(tpfl)
    return [ -norm(v) ; v[1] ; v[2] ; v[3] ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Intersection and emission point generator

    pgen( tpfl::DataType , N::Int )

This function generates a point `Xc` and `N` random emission points (in
a ``4×```N` matrix `X`) on the past null cone of `Xc`. Returns a tuple
`(X,Xc)`.

"""
function pgen( tpfl::DataType=Float64 , N::Int=4 )     
    v = (one(tpfl) + rand(tpfl))*vgenerator(tpfl)
    Xc = [ one(tpfl) + rand(tpfl) ; v[1] ; v[2] ; v[3] ]   # Generate Xc
    X = zeros(tpfl,4,N)     # Create container
    for i=1:N
        λ = (one(tpfl)+rand(tpfl))      # Affine parameter
        k = nullgen(tpfl)   # Null vector
        X[:,i] = Xc + λ*k   # Emission point
    end
    return (X,Xc)           # Return emission points and target point Xc
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Restricted intersection and emission point generator

    xgen( xc::Real , r1::Real , r2::Real , N::Int )

This function generates a point `Xc=[0;xc;0;0]` and `N` random emission
points (in a ``4×```N` matrix `X`) on the past null cone of `Xc` at a
radius `r` such that `r1<r<r2`. Returns a tuple `(X,Xc)`.

"""
function xgen( xc::Real , r1::Real , r2::Real=r1 , N::Int=4 )
    tpfl = typeof(xc)
	X = zeros(tpfl,4,N)	
	for i=1:N
		if r1 != r2
			r = (r2-r1)*rand(tpfl) + r1
		else
			r = r1
		end
		v = vgenerator(tpfl)
        x = r*v[1]
		y = r*v[2]
		z = r*v[3]
		X[1,i] = - norm(tpfl[(x - xc);y;z])
		X[2,i] = x
		X[3,i] = y
		X[4,i] = z
	end	#end for
	return (X,tpfl[0;xc;0;0])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Comparison function

    comp( q::Real , P::RealVec , Xc::RealVec )

This function compares the vectors `P` and `Xc` by taking the L1 norm of
the difference and comparing it with the L1 norm of `Xc`. If the ratio
`δ=|P-Xc|/|Xc|` is small compared to the threshold parameter `q`, the
function returns a tuple `(true,δ)`. The function returns `(false,δ)`
otherwise.

"""
function comp( q::Real , P::RealVec , Xc::RealVec )
    # Checks if location points P are recovered
    tpfl=typeof(q)

    δ = abs( norm(tpfl.(P-Xc),1) / norm(tpfl.(Xc),1) )

    if δ < abs(q)
        return ( true  , δ )
    else
        return ( false , δ )
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   Evaluation function
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Evaluation function

    main( locator::Function , N::Number , q::Real , k::Number ,     
          counter::Bool , usexgen::Bool )

This function tests the user specified `locator` function for `N`
stochastically generated test cases. The results produced by the
`locator` are compared (by way of the `comp` function) to the generated
intersection points up to a threshold value of `q`. The variable `k` is
the number of emission points to generate for each test case, and the
variable `usexgen` replaces the function `pgen` with `xgen` for test
case generation.

Examples:

    ceval.main(cereal.locatorfunc(4,"CFM10"),1e5,1e-6,4)

    ceval.main(cereal.locatorfunc(4,"FHC21"),1e5,1e-6,4)

    ceval.main(cereal.locatorfunc(5,"RTC21"),1e5,1e-9,5)

    ceval.main(cereal.locatorfunc(6,"RTC21"),1e5,5e-13,6)

"""
function main( locator::Function , N::Number , q::Real , k::Number=5 ,  
               usexgen::Bool=false )
    # Main test function--q determines floating point datatype
    tpfl    = typeof(q)
    lb      = false
    mi      = zero(Int)
    P       = zeros(tpfl,4)

    print(Int(N)," cases"," \n")
    print(k," emission points"," \n")

    for i=1:Int(N)
        if usexgen
            xc  = tpfl(2 + rand(tpfl))
	        Xp  = xgen(xc,tpfl(1))
        else
            Xp  = pgen(tpfl,k)
        end
        print("\r$i")
        P   = locator(Xp[1])
        if length(P) == 4
            Sd  = comp(q,P,Xp[2])
            if Sd[1] != true
                print("\n",Xp[1],"\n",Xp[2],"\n",Sd[2],"\n")
                lb  = true
                mi  += 1
            end
        elseif length(P) == 2 && typeof(P[1]) == typeof(P[2])
            Sda  = comp(q,P[1],Xp[2])
            Sdb  = comp(q,P[2],Xp[2])
            if Sda[1] != true && Sdb[1] != true
                print("\n",Xp[1],"\n",Xp[2],"\n",Sda[2],"\n",Sdb[2],
                      "\n")
                lb = true
                mi += 1
            end
        else
            print("\n","Check locator function"," \n")
            lb = true
            mi += 1        
        end
    end
    if lb
        print("\rTest ended with ",mi," failed cases out of "
                ,Int(N)," \n")
    else
        print("\rTest ended with zero failed cases out of "
                ,Int(N)," \n")
    end
end  # End full

#-----------------------------------------------------------------------
end     # End scope of module ceval
#-----------------------------------------------------------------------
