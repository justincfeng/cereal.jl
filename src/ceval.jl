
# INCOMPATIBLE WITH CURRENT CODE

#-----------------------------------------------------------------------
module ceval    # Evaluation module for cereal
#-----------------------------------------------------------------------

using LinearAlgebra

include("type.jl")
include("minkowski.jl")

#-----------------------------------------------------------------------
#   Point generators
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function vgenerator( tpfl::DataType )    # Generates random 3-vector
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
end  # End vgenerator

#-----------------------------------------------------------------------
function nullgen( tpfl::DataType )   
    # Generates random past-directed null vector
    v = vgenerator(tpfl)

    V = [ -norm(v) ; v[1] ; v[2] ; v[3] ]

    return V
end  # End nvgenerator

#-----------------------------------------------------------------------
function pgen( tpfl::DataType , N::Int=4 )     
    # Generates N random emission points on past null cone of Xc

    v = (one(tpfl) + rand(tpfl))*vgenerator(tpfl)
    
    Xc = [ one(tpfl) + rand(tpfl) ; v[1] ; v[2] ; v[3] ]   # Generate Xc

    X = zeros(tpfl,4,N)     # Create container

    for i=1:N
        λ = (one(tpfl)+rand(tpfl))      # Affine parameter
        k = nullgen(tpfl)   # Null vector
        X[:,i] = Xc + λ*k   # Emission point
    end

    return (X,Xc)           # Return emission points and target point Xc
end  # End pointgenerator

#-----------------------------------------------------------------------
function xgen( xc::Real , r1::Real , r2::Real=r1 , N::Int=4 )
    # Generates N random emission points on past null cone of [0;xc;0;0]
    # at radius r s.t. r1<r<r2
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

    v = zeros(tpfl,3,3)
    # Tetrahedron volume calculator 
    v[:,1] = X[2:4,2] - X[2:4,1]
    v[:,2] = X[2:4,3] - X[2:4,1]
    v[:,3] = X[2:4,4] - X[2:4,1]
    vol = abs(det([ v[:,1]  v[:,2]  v[:,3] ]))/6

	return (X,tpfl[0;xc;0;0],vol)
end	# End xgen

#-----------------------------------------------------------------------
#   Comparison function
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function compdirect( q::Real , P::RealVec , Xc::RealVec )
    # Checks if location points P are recovered
    tpfl=typeof(q)

    δ = abs( norm(tpfl.(P-Xc),1) / norm(tpfl.(Xc),1) )

    if δ < abs(q)
        return ( true  , δ )
    else
        return ( false , δ )
    end
end  # End compdirect

#-----------------------------------------------------------------------
function nrat( V::RealVec )        # Norm ratio
    return abs(η(V,V)/dot(V,V))
end  # End nrat

#-----------------------------------------------------------------------
function compnull( q::Real , P1::RealVec , P2::RealVec , X::RealMtx , 
                   Xc::RealVec )
    # Checks if separation vectors are null
    # This is not used; nrat is too sensitive to machine precision errs.
    tpfl=typeof(q)

    VT1 = [ nrat( P1 - X[:,1] ) ; nrat( P1 - X[:,2] ) ;
            nrat( P1 - X[:,3] ) ; nrat( P1 - X[:,4] ) ]

    VT2 = [ nrat( P2 - X[:,1] ) ; nrat( P2 - X[:,2] ) ;
            nrat( P2 - X[:,3] ) ; nrat( P2 - X[:,4] ) ]

    A1  = ( VT1[1] + VT1[2] + VT1[3] + VT1[4] ) / 4
    A2  = ( VT2[1] + VT2[2] + VT2[3] + VT2[4] ) / 4

    Δ   = Delta(X)
    χ   = Δ[2]

    if A1 < q && A2 < q
        return ( true  , X , A1 , A2 , Xc , P1 , P2 , Δ[1] ,
                η(χ,χ) )
    else
        return ( false , X , A1 , A2 , Xc , P1 , P2 , Δ[1] ,
                η(χ,χ) )
    end
end  # End compnull

#-----------------------------------------------------------------------
#   Test functions
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function full( iters::Number , q::Real , erm::Bool=true , 
               counter::Bool=true )
    # Main test function--q determines floating point datatype
    tpfl=typeof(q)
    lb=false
    mi = zero(Int64)

    for i=1:Int64(iters)
        B = true
        if counter
            print("\r$i")
        end
        Xp   = pgen(tpfl,4)
        P    = cereal.slocator(tpfl.(Xp[1]),erm)
        Sda  = compdirect(q,P[1],Xp[2])
        Sdb  = compdirect(q,P[2],Xp[2])
        Δ    = Delta( tpfl.(Xp[1]) )
        # Sn   = compnull(q,P[1],P[2],Xp[1],Xp[2])
        if Sda[1] != true && Sdb[1] != true # || Sn[1] != true
            print("\n",Xp[1],"\n",Xp[2],"\n",P[1],"\n",P[2],"\n",
                  Δ[1]," ",Δ[2]," ",Δ[3]," ",Δ[4]," ","\n",
                  sort([Sda[2],Sdb[2]])[1],"\n")
            lb = true
            mi += 1
        end
    end
    if lb
        print("\rTest ended with ",mi," failed cases out of "
              ,iters," \n")
    else
        print("\rTest ended with zero failed cases out of ",iters," \n")
    end
end  # End full

#-----------------------------------------------------------------------
function xtest( iters::Number , q::Real , erm::Bool=true , 
                counter::Bool=true )
    tpfl=typeof(q)
    lb=false
    mi = zero(Int64)

    for i=1:Int64(iters)
        B = true
        if counter
            print("\r$i")
        end
	    xc  = tpfl(2 + rand(tpfl))
	    X   = xgen(xc,tpfl(1))
        P   = cereal.slocator(tpfl.(X[1]),erm)
        Sda = compdirect(q,P[1],X[2])
        Sdb = compdirect(q,P[2],X[2])
        Δ   = Delta( X[1] )
        # Sn   = compnull(q,P[1],P[2],Xp[1],Xp[2])
        if Sda[1] != true && Sdb[1] != true # || Sn[1] != true
            print("\n",X[1],"\n",X[2],"\n",P[1],"\n",P[2],"\n",
                  Δ[1]," ",Δ[2]," ",Δ[3]," ",Δ[4]," ","\n",
                  X[3]," ",sort([Sda[2],Sdb[2]])[1],"\n")
            lb = true
            mi += 1
        end
    end
    if lb
        print("\rTest ended with ",mi," failed cases out of "
              ,iters," \n")
    else
        print("\rTest ended with zero failed cases out of ",iters," \n")
    end
end # End xtest

#-----------------------------------------------------------------------
function fullmulti( iters::Number , q::Real ,
                     N::Number=5 , erm::Bool=true , counter::Bool=true )
    # Main test function--q determines floating point datatype
    tpfl    = typeof(q)
    lb      = false
    mi      = zero(Int64)
    P       = zeros(tpfl,4)

    for i=1:Int64(iters)
        if counter
            print("\r$i")
        end
        Xp  = pgen(tpfl,N)
        P   = cereal.mlocator(tpfl.(Xp[1]),q,erm)
        Sd  = compdirect(q,P[1],Xp[2])
        if Sd[1] != true
            print("\n",Xp[1],"\n",Xp[2],"\n",Sd[2],"\n")
            lb  = true
            mi  += 1
        end
    end
    if lb
        print("\rTest ended with ",mi," failed cases out of "
                ,iters," \n")
    else
        print("\rTest ended with zero failed cases out of "
                ,iters," \n")
    end
end  # End fullmulti

#-----------------------------------------------------------------------
end     # End scope of module ceval
#-----------------------------------------------------------------------
