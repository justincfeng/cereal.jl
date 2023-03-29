#-----------------------------------------------------------------------
#   The functions presented here implement the relativistic location 
#   algorithm of the authors
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Transformation functions
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   LORENTZ TRANSFORMATION MATRIX CONSTRUCTOR
#-----------------------------------------------------------------------
"""
    LTM( NT::RealVec )

This function takes a timelike normal vector `NT` and constructs a 
Lorentz transformation matrix so that the transformed vector `NT'` has 
the form `NT'=[ N0 ; 0 ; 0 ; 0 ]`.

"""
function LTM( NT::RealVec )             # Lorentz transformation matrix
    tpfl=typeof(NT[1])
    normv = mnorm(NT)           # Normalization factor
    NV    = NT/normv            # Normalize Vector
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
#   ROTATION MATRIX CONSTRUCTOR
#-----------------------------------------------------------------------
"""
    MRz( v::RealVec )

This function takes a vector `v` and constructs a rotation matrix so
that the transformed vector `v'` has the form `v'=[ 0 ; ... ; 0 ; vz ]`.
The three-dimensional and four-dimensional cases are considered.

"""
function MRz( v::RealVec )                   # Rotate to z-adapted frame
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
            return Matrix(l*I(3))      
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
            return Matrix(l*I(4))
        end
    else
        return Matrix(l*I(4))
	end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   SPACETIME ROTATION OPERATOR
#-----------------------------------------------------------------------
"""
    Lrot( X::RealMtx , Λ::RealMtx )

This function applies the transformation transformation matrix `Λ` to
four emission points, which form the columns of a ``4×4`` matrix `X`. The
function returns a ``4×4`` matrix of transformed points `X'`.

"""
function Lrot( X::RealMtx , Λ::RealMtx )  # Spacetime rotation operator
    tpfl = typeof(X[1,1])

    Xprime = zeros(tpfl,4,4)

    for i=1:4
        Xprime[:,i] = Λ*X[:,i]
    end

    return Xprime
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   FLIP SPACELIKE VECTORS ACROSS A NULL DIRECTION
#-----------------------------------------------------------------------
"""
    NormflipS( Vsl::RealVec )

This function takes a spacelike vector `Vsl` and "flips" it across a
null direction to form a timelike vector. The null direction is chosen
to be the one tangent to the plane spanned by `Vsl` and the time axis
in the given coordinate system.

This function is used to construct a timelike vector from a spacelike
configuration vector (defined to be the vector normal to the plane 
spanned by the emission points).

"""
function NormflipS( Vsl::RealVec )
    tpfl = typeof(Vsl[1])

    normVsl = mnorm(Vsl)
    normsq = η(Vsl,Vsl)

    if Vsl[1] < 0                   # Ensures vector is future pointing
        Vsl = -Vsl
    end

    norms = norm(Vsl[2:4])          # Compute spatial norm

    if normsq > 0
        NVsl = Vsl/normVsl          # Normalization of spacetime vector
        nvsl = Vsl[2:4]/norms       # Normalized spatial part
        nlt  = NVsl[1]              # Temporal projection
        nls  = norm(NVsl[2:4])      # Spatial projection

        NVs  = zeros(tpfl,4)        # Construct and return flipped vec
        NVs[1]   = nls
        NVs[2:4] = nlt * nvsl
        return NVs
    else
        return Vsl/normVsl
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   INTERSECTION POINT FINDER (SPACELIKE CONFIGURATION HYPERPLANE)
#-----------------------------------------------------------------------
"""
# Intersection point finder (Spacelike configuration hyperplane)

    IPFinderS( Y::RealMtx ) 

Given a ``4×4`` matrix `Y` of emission points in an adapted frame 
(defined such that the emission points have the same time coordinate),
this function computes the coordinates of the intersection point `` X_c `` 
for the future pointing light cones of the emission points.
    
The procedure here is to find the circumcenter of a sphere which 
circumscribes the emission points. The spatial coordinates of the 
circumcenter yields the spatial coordinates of the intersection point
`` X_c `` and the time coordinate is given by the common time coordinate 
of the emission points plus the radius of the circumsphere.

"""
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
#   INTERSECTION POINT FINDER (TIMELIKE CONFIGURATION HYPERPLANE)
#-----------------------------------------------------------------------
"""
    IPFinderT( Y::RealMtx ) 

Given a ``4×4`` matrix `Y` of emission points in an adapted frame 
(defined such that the emission points have the same z coordinate),
this function computes the coordinates of the two intersection points 
`` X_c `` for the future pointing light cones of the emission points.
    
The procedure in this case is to find the vertex of the hyperboloid 
which passes through the emission points. The coordinates of the vertex 
on the `` (t,x,y) `` hyperplane yields the `` (t_c,x_c,y_c) `` coordinates of 
the  intersection points, and the Minkowski distance `` R `` between the 
vertex and the hyperboloid yields the `` z_c = z_e ± R `` coordinate of the 
intersection points `` X_c ``, where `` z_e `` is the common `` z `` coordinate
for the emission points in the adapted frame.

"""
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

function IPFinderTlin( Y::RealMtx )   
    # Finds intersection of light cones in adapted frame
    # Spacelike subconfiguration plane
    tpfl=typeof(Y[1,1])

    #   Defining variables
    z  = ( Y[4,1] + Y[4,2] + Y[4,3] + Y[4,4] )/4
    ( t1 , t2 , t3 ) = ( Y[1,1] , Y[1,2] , Y[1,3] )
    ( x1 , x2 , x3 ) = ( Y[2,1] , Y[2,2] , Y[2,3] )
    ( y1 , y2 , y3 ) = ( Y[3,1] , Y[3,2] , Y[3,3] )
    ( t4 , x4 , y4 ) = ( Y[1,4] , Y[2,4] , Y[3,4] )

    #   3x3 RTC-style matrix
    C  = [  t1-t2   x2-x1   y2-y1   ;
            t2-t3   x3-x2   y3-y2   ;
            t3-t4   x4-x3   y4-y3   ]

    B  = [  x2^2 - x1^2 + y2^2 - y1^2 + t1^2 - t2^2  ;
            x3^2 - x2^2 + y3^2 - y2^2 + t2^2 - t3^2  ;
            x4^2 - x3^2 + y4^2 - y3^2 + t3^2 - t4^2  ]

    #   Solving for the vertex
    Xc = inv(2 .* C) * B
    (tc , xc , yc) = (Xc[1], Xc[2], Xc[3])

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
#   FOUR POINT LOCATOR FUNCTION (FHC22)
#-----------------------------------------------------------------------
"""
    locator4FHC22( X::RealMtx )

This function implements the four point relativistic location algorithm
of the authors. It outputs a pair of location points.

"""
function locator4FHC22( X::RealMtx )
    #   Computes location from a single set of four emission points
    tpfl = typeof(X[1,1])
    XF = tpfl.(X)

    #   Containers
        nv  = zeros(tpfl,4)
        Xc  = zeros(tpfl,4)

        XT  = (zeros(tpfl,4,4), zeros(tpfl,4,4))

        Mz  = zeros(tpfl,4,4)
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
		print("nv norm zero.")
        return (Xc,Xc)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   FOUR POINT LOCATOR FUNCTION (FHC23 - development version)
#-----------------------------------------------------------------------
"""
    locator4FHC23( X::RealMtx )

This function implements the four point relativistic location algorithm
of the authors with modified procedure for timelike configuration 
plane. It outputs a pair of location points.

"""
function locator4FHC23( X::RealMtx )
    #   Computes location from a single set of four emission points
    tpfl = typeof(X[1,1])
    XF = tpfl.(X)

    #   Containers
        nv  = zeros(tpfl,4)
        Xc  = zeros(tpfl,4)

        XT  = (zeros(tpfl,4,4), zeros(tpfl,4,4))

        Mz  = zeros(tpfl,4,4)
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
        XT  = IPFinderTlin( Lrot(XF,Λ) )
        return (inv(Λ)*XT[1],inv(Λ)*XT[2])
    elseif η(nv,nv) == 0     # Null normal vector n
		print("nv norm zero.")
        return (Xc,Xc)
    end
end     #---------------------------------------------------------------

