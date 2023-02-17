#Experimental tool using the Broyden-based locator exclusively in flat spacetime

#-----------------------------------------------------------------------
#
#   BROYDEN ALGORITHM FROM SQUIRREL.JL
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    JiSMU( ΔF::RealVec , Δx::RealVec , Ji::RealMtx )

The function `JiSMU` implements the Sherman-Morrison update formula, 
returning an updated value of the Jacobian matrix `Ji`, given the 
respective differences `ΔF` and `Δx` for the function ``F(x)`` and its 
argument ``x``.

"""
function JiSMU( ΔF::RealVec , Δx::RealVec , Ji::RealMtx )
    ΔxTJi = transpose(Δx)*Ji
    return Ji + ( ( Δx - Ji*ΔF )/( ΔxTJi*ΔF ) )*( ΔxTJi )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    bsolve( F::Function , J::RealMtx , f0::RealVec , x0::RealVec ,
            nb::Int=24 )

The function `bsolve` implements the Broyden algorithm; in particular,
it finds the roots of the function `F(x)`, given an initial Jacobian 
matrix `J` and the initial guesses `f0` and `x0`. The parameter `nb` 
specifies the maximum number of iterations.

"""
function bsolve( F::Function , J::RealMtx , f0::RealVec , x0::RealVec ,
                 nb::Int=24 )
    tpfl = typeof(x0[1])

    Ji = inv(J)
    F0 = f0

    p = dot(F0,F0)

    ( b , k ) = ( true , 1 )

    if nb > 0
        Δx = - Ji*f0
        xi = x0 + Δx
        xs = xi
        if nb > 1
            Fi = F(xi)
            ΔF = Fi - F0
            F0 = Fi
            p  = norm(F0)
            ps = p
            i  = 2
            while i <= nb && b
                Ji = JiSMU( ΔF , Δx , Ji )
                Δx = - Ji*F0
                xi = xi + Δx
                Fi = F(xi)
                ( ΔF , F0 ) = ( Fi - F0 , Fi )
                pc = norm(F0)
                if pc<p && pc<ps
                    ( p , ps , xs ) = ( pc , pc , xi )
                elseif pc>p && k>=4
                    b = false
                else
                    p  = pc
                    k += 1
                end
                i += 1
            end
        end
    end

    return xs
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
# UNMODIFIED COMPONENT FUNCTIONS FROM SQUIRREL.JL LOCATOR
#
#-----------------------------------------------------------------------

function VidF( Zi::RealMtx )
    return vcat( Zi[6:8,1] , Zi[6:8,2] , Zi[6:8,3] , Zi[6:8,4] )
end     #---------------------------------------------------------------

function zFc( Zf::RealMtx )
    return vcat( Zf[1:4,1] - Zf[1:4,2] , Zf[1:4,1] - Zf[1:4,3] ,
                 Zf[1:4,1] - Zf[1:4,4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
# COMPONENT FUNCTIONS OF SQUIRREL.JL LOCATOR
# MODIFIED FOR STRAIGHT-LINE GEODESICS OF FLAT SPACETIME
#
#-----------------------------------------------------------------------

function V34null( V3::RealVec )
    tpfl = typeof(V3[1])
    return [ norm(V3) ; V3[1] ; V3[2] ; V3[3] ]
end     #end V34null

function gejacstr( Xi::RealVec , Vi::RealVec )
    v = Vi[2:4]
    result = DiffResults.JacobianResult(vcat(Xi,Vi),v)
    result = ForwardDiff.jacobian!(result,var->vcat(Xi+V34null(var),V34null(var))
                                    , v )
    return ( DiffResults.value(result)
             , DiffResults.jacobian(result)[1:4,:] )
end     #end gejacstr

function zFstr( Vid::RealVec , Zi::RealMtx )
    # Want to find roots of this wrt initial velocity
    tpfl=typeof(Vid[1])

    Xf = [zeros(tpfl,4) for _ in 1:4]

    for i=1:4
        ( k1 , k2 ) = ( 1 + 3*(i-1) , 3 + 3*(i-1) )
        Xf[i] = Zi[1:4,i] + V34null( Vid[k1:k2] ) 
    end

    return vcat( Xf[1] - Xf[2] , Xf[1] - Xf[3] , Xf[1] - Xf[4] )
end     #end zFstr

function geocJstr( Zi::RealMtx )
    tpfl=typeof(Zi[1,1])
    v   = zeros(tpfl,3)
    J   = zeros(tpfl,12,12)
    dXV = [zeros(tpfl,4,3) for _ in 1:4]
    Zf  = copy(Zi)

    for i=1:4
        res = gejacstr( Zi[1:4,i] , Zi[5:8,i] )
        ( Zf[:,i] , dXV[i] ) = ( res[1] , res[2] )
    end

    J[1:4,1:3]  = dXV[1]
    J[5:8,1:3]  = dXV[1]
    J[9:12,1:3] = dXV[1]

    J[1:4,4:6]    = - dXV[2]
    J[5:8,7:9]    = - dXV[3]
    J[9:12,10:12] = - dXV[4]

    return (Zf,J)
end     #end geocJstr

function idfstr( Zi::RealMtx , nb::Int )
    tpfl=typeof(Zi[1,1])

    Zf = copy(Zi)
    V0 = VidF(Zi)
    res = geocJstr( Zi )

    ( F0 , J ) = ( zFc(res[1]) , res[2] )

    V = bsolve( v->zFstr(v,Zi) , J , F0 , V0 , nb )

    Zf[5:8,1] = V34null(  V[1:3]  )
    Zf[5:8,2] = V34null(  V[4:6]  )
    Zf[5:8,3] = V34null(  V[7:9]  )
    Zf[5:8,4] = V34null( V[10:12] )

    return Zf
end     #end idfstr

#-----------------------------------------------------------------------
#
# FOUR-POINT LOCATOR FUNCTION USING BROYDEN ALGORITHM IN FLAT SPACETIME
#
#-----------------------------------------------------------------------
function locator4flbr( X::RealMtx , Xc::RealVec , nb::Int=24 , idv::Bool=false , 
                   V::RealMtx=zeros(Float64,4,4) )
    tpfl=typeof(X[1,1])

    Zi  = zeros(tpfl,8,4)
    Zf = copy(Zi)

    if idv 
        for i=1:4
            Zi[:,i] = vcat( X[:,i] , V[:,i] )
        end
    else
        for i=1:4
            Zi[:,i] = vcat( X[:,i] , Xc - X[:,i] )
        end
    end

    Zi = idfstr( Zi , nb )
    for i=1:4
        Zf[:,i] = vcat(Zi[1:4,i] + Zi[5:8,i], Zi[5:8,i])
    end
    return (Zf[1:4,1] + Zf[1:4,2] + Zf[1:4,3] + Zf[1:4,4])/4 
end     #end locator4flbr
