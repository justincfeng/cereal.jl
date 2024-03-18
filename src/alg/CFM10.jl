#-----------------------------------------------------------------------
#   The functions presented here implement the relativistic location
#   formula given by Coll, Ferrando, and Morales-Lladosa in 
#   Coll et al., Class.Quant.Grav. 27 (2010) 065013
#   and
#   Coll et al., Phys. Rev. D 86, 084036 (2012).
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   FRAME FUNCTION
#-----------------------------------------------------------------------
"""
    Frame( X::RealMtx )

This function constructs a matrix of frame vectors `` e_i `` formed by the
differences between the emission points `` X_I ``:

`` e_i = X_i - X_4 ``,

where the index `` i `` runs from ``1`` through ``3``.

"""
function Frame( X::RealMtx ) # Constructs spatial frame
    tpfl=typeof(X[1,1])
    s = 1
    E = zeros(tpfl,4,3)
    for i=1:3
        E[:,i] = X[:,i] - X[:,4]
    end
    return E
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   CONFIGURATION VECTOR CONSTRUCTOR
#-----------------------------------------------------------------------
"""
    ConfVec( E::RealMtx )

This function constructs the configuration vector `` χ ``, which is 
normal to the hyperplane spanned by the three vectors `` e_1 ``, `` e_2 ``,
`` e_3 ``.

"""
function ConfVec( E::RealMtx )
    tpfl=typeof(E[1,1])

    return -HodgeV( E[:,1] , E[:,2] , E[:,3] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   AUXILIARY VECTOR CONSTRUCTOR
#-----------------------------------------------------------------------
"""
    AuxVec( χ::RealVec )

This function constructs the vector `` ξ ``, which is arbitrary except
for the requirement that it be transverse to `` χ ``. Here, the vector 
`` ξ `` is constructed so that `` η(ξ,χ) = 1``.

"""
function AuxVec( χ::RealVec )
    tpfl=typeof(χ[1])

    ξ = ones(tpfl,4)

    if χ[1] != zero(tpfl)
        ξ[1] = ( 1.0 + (χ[2]*ξ[2] + χ[3]*ξ[3] + χ[4]*ξ[4]) ) / χ[1]
    else
        ξ = -ξ/(χ[2]*ξ[2] + χ[3]*ξ[3] + χ[4]*ξ[4])
    end

    return -ξ
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   y* CONSTRUCTOR
#-----------------------------------------------------------------------
"""
    ystar( E::RealMtx , χ::RealVec , ξ::RealVec )

This function computes the vector `` y_* ``:

`` y_* = (η(ξ,χ))^{-1} i_ξ H ``,

where `` H `` is a bivector given by (`h(⋅,⋅)=-ϵ(0,0,⋅,⋅)`):

`` H = Ω_1 h(e_2,e_3) + Ω_2 h(e_3,e_1) + Ω_3 h(e_1,e_2) ``,

and `` Ω_i := η(e_i,e_i) ``.

"""
function ystar( E::RealMtx , χ::RealVec , ξ::RealVec )
    tpfl=typeof(E[1,1])

    h23 = -ϵ( 0 , 0 , E[:,2] , E[:,3] )
    h31 = -ϵ( 0 , 0 , E[:,3] , E[:,1] )
    h12 = -ϵ( 0 , 0 , E[:,1] , E[:,2] )

    Ω1 = η(E[:,1],E[:,1])/2
    Ω2 = η(E[:,2],E[:,2])/2
    Ω3 = η(E[:,3],E[:,3])/2

    H  = Ω1*h23 + Ω2*h31 + Ω3*h12

    ys = zeros(tpfl,4)

    for i=1:4
        Z  = zero(tpfl)
        for j=1:4
            Z += ξ[j]*H[j,i]
        end
        ys[i] = Z
    end

    ys = ys / η(ξ,χ)

    ys[1] = -ys[1]

    return ys
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   DELTA
#-----------------------------------------------------------------------
"""
    Delta( χ::RealVec , ys::RealVec )

This function takes two vectors and inputs, and computes the following
scalar quantity

`` Δ := η(y_*,χ)^2 - η(y_*,y_*) η(χ,χ) ``.

"""
function Delta( χ::RealVec , ys::RealVec )
    return η(ys,χ)^2 - η(ys,ys)*η(χ,χ)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   FOUR POINT LOCATOR FUNCTION (CFM10)
#-----------------------------------------------------------------------
"""
    locator4CFM10( X::RealMtx )

The function presented here implements the four point relativistic
location formula given by Coll, Ferrando, and Morales-Lladosa. It
outputs a pair of location points.

"""
function locator4CFM10( X::RealMtx )     # Locator function
    tpfl=typeof(X[1,1])

    E   = Frame( X )
    χ   = ConfVec( E )
    ys  = ystar( E , χ , AuxVec(χ) )

    Δ   = Delta( χ , ys )
    yχ  = η(ys,χ)

    d1  = ( η(ys,χ) + sqrt( abs(Δ) ) )
    d2  = ( η(ys,χ) - sqrt( abs(Δ) ) )

    if d1 != 0
        φ1 = η(ys,ys) / d1
        x1 = X[:,4] + ys - φ1*χ
    else
        x1 = zeros(tpfl,4)
    end

    if d2 != 0
        φ2 = η(ys,ys) / d2
        x2 = X[:,4] + ys - φ2*χ
    else
        x2 = zeros(tpfl,4)
    end

    return ( x1 , x2 )
end     #---------------------------------------------------------------
