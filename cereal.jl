#-----------------------------------------------------------------------
module cereal      # cereal module
#-----------------------------------------------------------------------

using LinearAlgebra
using Combinatorics
using Statistics

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

#-----------------------------------------------------------------------
#   Minkowski product and norm
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function η( V1::RealVec , V2::RealVec )  # Computes Minkowski product
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
end  # End η

#-----------------------------------------------------------------------
function mnorm( V )  # Computes Minkowski norm
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
end  # End mnorm

#-----------------------------------------------------------------------
#   Levi-Civita and Hodge
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function ϵ( a , b , c , d )  # Levi-Civita
    if typeof.((a,b,c,d)) == (Int,Int,Int,Int)
        return levicivita([a,b,c,d])
    elseif typeof(a) <: RealVec && typeof(b) <: RealVec &&
        typeof(c) <: RealVec && typeof(d) <: RealVec
        v = zero( typeof(a[1]) )
        for i=1:4, j=1:4, k=1:4, l=1:4
            v += levicivita([i,j,k,l])*a[i]*b[j]*c[k]*d[l]
        end
        return v
    else
        return 0
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function multivec( X::RealMtx , k::Int , vecvec::Bool )  
    # Multi vector function
    tpfl = typeof(X[1,1])
    d    = size(X)[1]
    np   = size(X)[2]
    Z    = [zeros(tpfl,d) for _ in 1:np]

    for i=1:np
        Z[i] = X[:,i]
    end

    w   = collect(combinations(Z,k))

    if vecvec
        return w
    else
        l   = length(w)
        z   = [zeros(tpfl,d,k) for _ in 1:l]
        for i=1:l, j=1:k 
                z[i][:,j] = w[i][j]
        end

        return z
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function HodgeV( U::RealVec , V::RealVec , W::RealVec )
    # Hodge of three vectors
    tpfl=typeof(U[1]) 

    v = zeros(tpfl,4)
    for i=1:4, j=1:4, k=1:4, l=1:4
        v[i] += ϵ(i,j,k,l)*U[j]*V[k]*W[l]
    end

    v[1] = -v[1]    # Raise index

    return v
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Hodge2( U::RealVec , V::RealVec )
    tpfl=typeof(U[1]) 

    W = zeros(tpfl,4,4)
    for i=1:4
    for j=1:4
        Z  = zero(tpfl)
        for k=1:4
        for l=1:4
            Z += -cereal.ϵ(i,j,k,l)*U[k]*V[l]
        end
        end
        W[i,j] = Z
    end
    end

    return W
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   Transformation functions
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function LTM( NVa::RealVec )     
    # Constructs Lorentz transformation matrix
    tpfl=typeof(NVa[1])
    normv = mnorm(NVa)  # Normalization factor
    NV    = NVa/normv             # Normalize Vector
    l     = one(tpfl)           # l has a value of one
    γ     = NV[1]               # γ is the Lorentz factor
    if γ<0
        NV = -NV
        γ  = -γ
    end
        δ    = (γ - l)             # A useful quantity
        norms = norm(NV[2:4])
    if γ != l || norms != zero(tpfl) || δ != zero(tpfl)
        nf = l/norms    # Normalization factor
        β = abs(norms/γ) # β is v/c

        vx = nf*NV[2]           # x-component of unit vector
        vy = nf*NV[3]           # y-component of unit vector
        vz = nf*NV[4]           # z-component of unit vector
        return [  γ        -γ*β*vx       -γ*β*vy       -γ*β*vz      ;
                 -γ*β*vx   l + δ*(vx^2)  δ*vx*vy       δ*vx*vz      ;
                 -γ*β*vy   δ*vy*vx       l + δ*(vy^2)  δ*vy*vz      ;
                 -γ*β*vz   δ*vz*vx       δ*vz*vy       l + δ*(vz^2) ]
    else
        return Matrix(l*I(4))
    end
end  # End LTM

#-----------------------------------------------------------------------
function MRz( v::RealVec )      # Rotate to z-adapted frame
    tpfl=typeof(v[1])
    l = one(tpfl)
    o = zero(tpfl)
	if length(v)==3             # 3×3 rotation matrix
		rv=norm(v[1:3])
		Rv=norm(v[1:2])
        if Rv > 0 && rv > 0
		    Cosϑ=v[3]/rv
		    Sinϑ=Rv/rv
		    Cosφ=v[1]/Rv
		    Sinφ=v[2]/Rv
		    return [ Cosϑ*Cosφ  Cosϑ*Sinφ   -Sinϑ   ;
                     -Sinφ      Cosφ        o       ;
                     Sinϑ*Cosφ  Sinϑ*Sinφ   Cosϑ    ]
        else
            return I(3)        
        end
    elseif length(v)==4         # 4×4 rotation matrix
        rv=norm(v[2:4])
		Rv=norm(v[2:3])
        if Rv > 0 && rv > 0
            Cosϑ=v[4]/rv
            Sinϑ=Rv/rv
		    Cosφ=v[2]/Rv
		    Sinφ=v[3]/Rv
            return [ one(tpfl)  o           o           o       ;
                        o       Cosϑ*Cosφ   Cosϑ*Sinφ   -Sinϑ   ;
                        o       -Sinφ       Cosφ        o       ;
                        o       Sinϑ*Cosφ   Sinϑ*Sinφ   Cosϑ    ]
        else
            return I(4)
        end
    else
        return I(4)
	end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Lrot( X::RealMtx , Λ::RealMtx )  # Spacetime rotation operator
    tpfl = typeof(X[1,1])

    Xprime = zeros(tpfl,4,4)

    for i=1:4
        Xprime[:,i] = Λ*X[:,i]
    end

    return Xprime
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function NormflipS( Vsl::RealVec )    # Flip vectors
    tpfl = typeof(Vsl[1])

    normVsl = mnorm(Vsl)
    normsq = η(Vsl,Vsl)

    if Vsl[1] < 0
        Vsl = -Vsl
    end

    norms = norm(Vsl[2:4])

    if normsq > 0
        NVsl = Vsl/normVsl
        nlt  = NVsl[1]
        nls  = norm(NVsl[2:4])
        nvsl = Vsl[2:4]/norms
        NVs  = zeros(tpfl,4)
        NVs[1]   = nls
        NVs[2:4] = nlt * nvsl
        return NVs
    else
        return Vsl/normVsl
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   Adapted frame intersection point finders
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function IPFinderS( Y::RealMtx )   
    # Finds intersection of light cones in adapted frame
    tpfl=typeof(Y[1,1])

    x = zeros(tpfl,3,4)     # Create containers
    v = zeros(tpfl,3,3)

    x[1,:] = Y[2,:]         # Spatial points for corners of tetrahedron
    x[2,:] = Y[3,:]
    x[3,:] = Y[4,:]

    v[:,1] = x[:,2] - x[:,1]    # "Frame" vectors centered on x[:,1]
    v[:,2] = x[:,3] - x[:,1]
    v[:,3] = x[:,4] - x[:,1]
    
    B = [ x[:,2]⋅x[:,2] - x[:,1]⋅x[:,1] ; x[:,3]⋅x[:,3] - x[:,1]⋅x[:,1];
          x[:,4]⋅x[:,4] - x[:,1]⋅x[:,1] ] / 2    # Constructing B vector

    A = transpose([ v[:,1]  v[:,2]  v[:,3] ])    # Constructing A matrix

    xc = inv(A)*B           # Compute circumcenter
    rc = (norm(xc-x[:,1]) + norm(xc-x[:,2]) + # Distance to circumcenter
          norm(xc-x[:,3]) + norm(xc-x[:,4])) / 4
    tc = rc + Y[1,1]        # Compute time coordinate
    return [ tc ; xc[1] ; xc[2] ; xc[3] ]
end  # End IPfinder

#-----------------------------------------------------------------------
function IPFinderT( Y::RealMtx )   
    # Finds intersection of light cones in adapted frame
    # Spacelike subconfiguration plane
    tpfl=typeof(Y[1,1])

    #   Defining variables
    z  = ( Y[4,1] + Y[4,2] + Y[4,3] + Y[4,4] )/4
    ( t1 , t2 , t3 ) = ( Y[1,1] , Y[1,2] , Y[1,3] )
    ( x1 , x2 , x3 ) = ( Y[2,1] , Y[2,2] , Y[2,3] )
    ( y1 , y2 , y3 ) = ( Y[3,1] , Y[3,2] , Y[3,3] )
    ( t4 , x4 , y4 ) = ( Y[1,4] , Y[2,4] , Y[3,4] )

    #   Denominator for time coordinate:
    dT = -(t3*x2*y1) + t4*x2*y1 + t2*x3*y1 - t4*x3*y1 - t2*x4*y1 + 
        t3*x4*y1 + t3*x1*y2 - t4*x1*y2 - t1*x3*y2 + t4*x3*y2 + 
        t1*x4*y2 - t3*x4*y2 - t2*x1*y3 + t4*x1*y3 + t1*x2*y3 - 
        t4*x2*y3 - t1*x4*y3 + t2*x4*y3 + t2*x1*y4 - t3*x1*y4 - 
        t1*x2*y4 + t3*x2*y4 + t1*x3*y4 - t2*x3*y4

    #   Denominator for spatial coordinates:
    dS = t3*x2*y1 - t4*x2*y1 - t2*x3*y1 + t4*x3*y1 + t2*x4*y1 - 
         t3*x4*y1 - t3*x1*y2 + t4*x1*y2 + t1*x3*y2 - t4*x3*y2 - 
         t1*x4*y2 + t3*x4*y2 + t2*x1*y3 - t4*x1*y3 - t1*x2*y3 + 
         t4*x2*y3 + t1*x4*y3 - t2*x4*y3 - t2*x1*y4 + t3*x1*y4 + 
         t1*x2*y4 - t3*x2*y4 - t1*x3*y4 + t2*x3*y4
    
    #   t coordinate:
    tc  = ( -(t3^2*x2*y1) + t4^2*x2*y1 + t2^2*x3*y1 - t4^2*x3*y1 - 
            x2^2*x3*y1 + x2*x3^2*y1 - t2^2*x4*y1 + t3^2*x4*y1 + 
            x2^2*x4*y1 - x3^2*x4*y1 - x2*x4^2*y1 + x3*x4^2*y1 + 
            t3^2*x1*y2 - t4^2*x1*y2 - t1^2*x3*y2 + t4^2*x3*y2 + 
            x1^2*x3*y2 - x1*x3^2*y2 + t1^2*x4*y2 - t3^2*x4*y2 - 
            x1^2*x4*y2 + x3^2*x4*y2 + x1*x4^2*y2 - x3*x4^2*y2 + 
            x3*y1^2*y2 - x4*y1^2*y2 - x3*y1*y2^2 + x4*y1*y2^2 - 
            t2^2*x1*y3 + t4^2*x1*y3 + t1^2*x2*y3 - t4^2*x2*y3 - 
            x1^2*x2*y3 + x1*x2^2*y3 - t1^2*x4*y3 + t2^2*x4*y3 + 
            x1^2*x4*y3 - x2^2*x4*y3 - x1*x4^2*y3 + x2*x4^2*y3 - 
            x2*y1^2*y3 + x4*y1^2*y3 + x1*y2^2*y3 - x4*y2^2*y3 + 
            x2*y1*y3^2 - x4*y1*y3^2 - x1*y2*y3^2 + x4*y2*y3^2 + 
            t2^2*x1*y4 - t3^2*x1*y4 - t1^2*x2*y4 + t3^2*x2*y4 + 
            x1^2*x2*y4 - x1*x2^2*y4 + t1^2*x3*y4 - t2^2*x3*y4 - 
            x1^2*x3*y4 + x2^2*x3*y4 + x1*x3^2*y4 - x2*x3^2*y4 + 
            x2*y1^2*y4 - x3*y1^2*y4 - x1*y2^2*y4 + x3*y2^2*y4 + 
            x1*y3^2*y4 - x2*y3^2*y4 - x2*y1*y4^2 + x3*y1*y4^2 + 
            x1*y2*y4^2 - x3*y2*y4^2 - x1*y3*y4^2 + x2*y3*y4^2 ) /
          (2*dT)
    
    #   x coordinate:
    xc = ( -(t2^2*t3*y1) + t2*t3^2*y1 + t2^2*t4*y1 - t3^2*t4*y1 - 
            t2*t4^2*y1 + t3*t4^2*y1 + t3*x2^2*y1 - t4*x2^2*y1 - 
            t2*x3^2*y1 + t4*x3^2*y1 + t2*x4^2*y1 - t3*x4^2*y1 + 
            t1^2*t3*y2 - t1*t3^2*y2 - t1^2*t4*y2 + t3^2*t4*y2 + 
            t1*t4^2*y2 - t3*t4^2*y2 - t3*x1^2*y2 + t4*x1^2*y2 + 
            t1*x3^2*y2 - t4*x3^2*y2 - t1*x4^2*y2 + t3*x4^2*y2 - 
            t3*y1^2*y2 + t4*y1^2*y2 + t3*y1*y2^2 - t4*y1*y2^2 - 
            t1^2*t2*y3 + t1*t2^2*y3 + t1^2*t4*y3 - t2^2*t4*y3 - 
            t1*t4^2*y3 + t2*t4^2*y3 + t2*x1^2*y3 - t4*x1^2*y3 - 
            t1*x2^2*y3 + t4*x2^2*y3 + t1*x4^2*y3 - t2*x4^2*y3 + 
            t2*y1^2*y3 - t4*y1^2*y3 - t1*y2^2*y3 + t4*y2^2*y3 - 
            t2*y1*y3^2 + t4*y1*y3^2 + t1*y2*y3^2 - t4*y2*y3^2 + 
            t1^2*t2*y4 - t1*t2^2*y4 - t1^2*t3*y4 + t2^2*t3*y4 + 
            t1*t3^2*y4 - t2*t3^2*y4 - t2*x1^2*y4 + t3*x1^2*y4 + 
            t1*x2^2*y4 - t3*x2^2*y4 - t1*x3^2*y4 + t2*x3^2*y4 - 
            t2*y1^2*y4 + t3*y1^2*y4 + t1*y2^2*y4 - t3*y2^2*y4 - 
            t1*y3^2*y4 + t2*y3^2*y4 + t2*y1*y4^2 - t3*y1*y4^2 - 
            t1*y2*y4^2 + t3*y2*y4^2 + t1*y3*y4^2 - t2*y3*y4^2) /
          (2*dS)
    
    #   y coordinate:
    yc = (  t2^2*t3*x1 - t2*t3^2*x1 - t2^2*t4*x1 + t3^2*t4*x1 + 
            t2*t4^2*x1 - t3*t4^2*x1 - t1^2*t3*x2 + t1*t3^2*x2 + 
            t1^2*t4*x2 - t3^2*t4*x2 - t1*t4^2*x2 + t3*t4^2*x2 + 
            t3*x1^2*x2 - t4*x1^2*x2 - t3*x1*x2^2 + t4*x1*x2^2 + 
            t1^2*t2*x3 - t1*t2^2*x3 - t1^2*t4*x3 + t2^2*t4*x3 + 
            t1*t4^2*x3 - t2*t4^2*x3 - t2*x1^2*x3 + t4*x1^2*x3 + 
            t1*x2^2*x3 - t4*x2^2*x3 + t2*x1*x3^2 - t4*x1*x3^2 - 
            t1*x2*x3^2 + t4*x2*x3^2 - t1^2*t2*x4 + t1*t2^2*x4 + 
            t1^2*t3*x4 - t2^2*t3*x4 - t1*t3^2*x4 + t2*t3^2*x4 + 
            t2*x1^2*x4 - t3*x1^2*x4 - t1*x2^2*x4 + t3*x2^2*x4 + 
            t1*x3^2*x4 - t2*x3^2*x4 - t2*x1*x4^2 + t3*x1*x4^2 + 
            t1*x2*x4^2 - t3*x2*x4^2 - t1*x3*x4^2 + t2*x3*x4^2 + 
            t3*x2*y1^2 - t4*x2*y1^2 - t2*x3*y1^2 + t4*x3*y1^2 + 
            t2*x4*y1^2 - t3*x4*y1^2 - t3*x1*y2^2 + t4*x1*y2^2 + 
            t1*x3*y2^2 - t4*x3*y2^2 - t1*x4*y2^2 + t3*x4*y2^2 + 
            t2*x1*y3^2 - t4*x1*y3^2 - t1*x2*y3^2 + t4*x2*y3^2 + 
            t1*x4*y3^2 - t2*x4*y3^2 - t2*x1*y4^2 + t3*x1*y4^2 + 
            t1*x2*y4^2 - t3*x2*y4^2 - t1*x3*y4^2 + t2*x3*y4^2 ) /
          (2*dS)

    #   Hyperboloid distance calculation
    R  = (  mnorm( [t1-tc;x1-xc;y1-yc] ) + 
            mnorm( [t2-tc;x2-xc;y2-yc] ) +
            mnorm( [t3-tc;x3-xc;y3-yc] ) +
            mnorm( [t4-tc;x4-xc;y4-yc] ) ) / 4

    #   z coordinate
    Δz = R

    return  (   [ tc ; xc ; yc ; z + Δz ] 
              , [ tc ; xc ; yc ; z - Δz ] )
end  # End IPfinder

#-----------------------------------------------------------------------
#   Locator functions
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function slocator( X::RealMtx , erm::Bool=true )
    #   Computes location from a single set of four emission points
    tpfl = typeof(X[1,1])
    XF = tpfl.(X)

    #   Containers
        nv  = zeros(tpfl,4)
        Xc  = zeros(tpfl,4)

        XA = zeros(tpfl,4,4)
        XB = zeros(tpfl,4,4)

        Mz  = zeros(tpfl,4,4)
        My  = zeros(tpfl,4,4)
        Λ   = zeros(tpfl,4,4)

        nv = HodgeV(X[:,1] - X[:,4],X[:,2] - X[:,4],X[:,3] - X[:,4])

    if η(nv,nv) < 0         # Timelike normal vector n
        #---------------------------------------------------------------
        #   Computation for spacelike configuration plane
        #---------------------------------------------------------------
        Λ  = LTM(nv)
        Xc = inv(Λ)*IPFinderS( Lrot(XF,Λ) )
        return (Xc,Xc)
    elseif η(nv,nv) > 0     # Spacelike normal vector n
        #---------------------------------------------------------------
        #   Computation for timelike configuration plane
        #---------------------------------------------------------------
        Λn  = LTM(NormflipS(nv))
        Mz  = MRz(nv)
        Λ   = Mz*Λn
        XT  = IPFinderT( Lrot(XF,Λ) )
        return (inv(Λ)*XT[1],inv(Λ)*XT[2])
    elseif η(nv,nv) == 0     # Null normal vector n
        if erm
		    print("nv norm zero.")
        end
        return (Xc,Xc)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function mlocator( X::RealMtx , q::Real=1e-14 , erm::Bool=true )
    #   Calculates location for more than four emission points
    tpfl = typeof(q*X[1,1])

    l = size(X)

    XF = tpfl.(X)

    if l[2] == 4
        return slocator(XF,erm)
    elseif l[2] > 4
        W   = multivec(X,4,false)

        k = length(W)
        Xa  = [(zeros(tpfl,4),zeros(tpfl,4)) for _ =1:k ]
        VV  = (zeros(tpfl,4),zeros(tpfl,4))
        w   = [zeros(tpfl,4) for _ =1:2*k ]

        for a=1:k
            (w[a],w[k+a]) = slocator( W[a] , erm )
        end

        #   The following picks out points closely clustered together
            k = length(w)
            ΔVm   = [zeros(tpfl,4) for _ =1:k ]
            ΔVp   = [zeros(tpfl,4) for _ =1:k ]
            δV    = zeros(tpfl,k)

            #   Sort points according to their norm
            wn     = norm.(w)
            iws    = sortperm(wn)
            ws     = w[iws]
        
            #   Compute differences between points
            j = 0
            for i=2:k-1
                ΔVm[i] = (ws[i] - ws[i-1])
                ΔVp[i] = (ws[i+1] - ws[i])
                δV[i]  = sort( [dot(ΔVm[i],ΔVm[i]),
                                dot(ΔVp[i],ΔVp[i]) ] )[1]
        
                if dot(ΔVm[i],ΔVm[i]) < q || dot(ΔVp[i],ΔVp[i]) < q
                    j += 1
                end
        
                if i == 2 && dot(ΔVm[i],ΔVm[i]) < q
                    j += 1    
                elseif i == k-1 && dot(ΔVp[i],ΔVp[i]) < q
                    j += 1
                end
            end
            δV[1] = dot(ΔVm[2],ΔVm[2])
            δV[k] = dot(ΔVp[k-1],ΔVp[k-1])

            #   Return index vector for sorting of differences
            is    = sortperm(δV)

            if j==0
                j=1
            end

            #   The following should return vectors closely clustered
            WS = ws[is][1:j]

        #   Re-sort location points according to smallest Minkowski norm
            δW = zeros(tpfl,j)
            for a=1:j
                for i=1:4
                    δW[a] += abs(η(X[:,i]-WS[a],X[:,i]-WS[a]))
                end
            end
            iW = sortperm(δW)

            return (WS[iW][1],zeros(tpfl,4))
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # End scope of module cereal
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
module ceval    # Evaluation module for cereal
#-----------------------------------------------------------------------

using LinearAlgebra
import ..cereal         # Importing functions from cereal

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

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
    return abs(cereal.η(V,V)/dot(V,V))
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

    Δ   = cereal.Delta(X)
    χ   = Δ[2]

    if A1 < q && A2 < q
        return ( true  , X , A1 , A2 , Xc , P1 , P2 , Δ[1] ,
                cereal.η(χ,χ) )
    else
        return ( false , X , A1 , A2 , Xc , P1 , P2 , Δ[1] ,
                cereal.η(χ,χ) )
    end
end  # End compnull

#-----------------------------------------------------------------------
function Delta( X::RealMtx )
    tpfl=typeof(X[1,1])
    E = zeros(tpfl,4,3)
    for i=1:3
        E[:,i] = X[:,i] - X[:,4]
    end

    χ = cereal.HodgeV(E[:,1],E[:,2],E[:,3])

    ξ = ones(tpfl,4)

    if χ[1] != 0.0
        ξ[1] = ( 1.0 - (χ[2]*ξ[2] + χ[3]*ξ[3] + χ[4]*ξ[4]) ) / χ[1]
    else
        ξ = ξ/(χ[2]*ξ[2] + χ[3]*ξ[3] + χ[4]*ξ[4])
    end

    h23 = cereal.Hodge2( E[:,2] , E[:,3] )
    h31 = cereal.Hodge2( E[:,3] , E[:,1] )
    h12 = cereal.Hodge2( E[:,1] , E[:,2] )

    Ω1 = cereal.η(E[:,1],E[:,1])/2
    Ω2 = cereal.η(E[:,2],E[:,2])/2
    Ω3 = cereal.η(E[:,3],E[:,3])/2

    H  = Ω1*h23 + Ω2*h31 + Ω3*h12

    ys = zeros(tpfl,4)

    for i=1:4
    Z  = zero(tpfl)
    for j=1:4
        Z += ξ[j]*H[j,i]
    end
    ys[i] = Z
    end

    ys = ys / cereal.η(ξ,χ)

    ys[1] = -ys[1]

    Δ = cereal.η(ys,χ)^2 - cereal.η(ys,ys)*cereal.η(χ,χ)

    den = abs(cereal.η(ys,ys))

    return (Δ/den,abs(cereal.η(ys,χ))/den,cereal.η(χ,χ)/den,
            cereal.η(χ,χ)*den/(abs(cereal.η(ys,χ))^2))
end     #---------------------------------------------------------------

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
