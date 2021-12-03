#-----------------------------------------------------------------------
#   Levi-Civita
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
Levi-Civita symbol and tensor

    ϵ( a , b , c , d )

The function `ϵ` computes the Levi-Civita symbol and implements the
corresponding tensor. The arguments `( a , b , c , d )` can either be
integers or vectors; in the case of integers, it returns the components
of the Levi-Civita symbol, and in the case of vectors, it returns the
contraction of the Levi-Civita symbol with the vectors. 

This function also returns Hodge duals of vector tensor products. For
instance, the following

    ϵ( 0 , 0 , U , V )

yields the rank-2 quantity ``ϵ_{i,j,k,l} U^i V^i``. Make sure all `0`
valued arguments come first.

"""
function ϵ( a , b , c , d )              # Levi-Civita symbol and tensor
    if a == 0 && typeof(b) <: RealVec && typeof(c) <: RealVec && 
        typeof(d) <: RealVec
        v = zeros( typeof(d[1]) , 4 )
        for i=1:4, j=1:4, k=1:4, l=1:4
            v[i] += levicivita([i,j,k,l])*b[j]*c[k]*d[l]
        end
        return v
    elseif a == 0 && b == 0 && typeof(c) <: RealVec && 
        typeof(d) <: RealVec
        v = zeros( typeof(d[1]) , 4 , 4 )
        for i=1:4, j=1:4, k=1:4, l=1:4
            v[i,j] += levicivita([i,j,k,l])*c[k]*d[l]
        end
        return v
    elseif a == 0 && b == 0 && c == 0 && typeof(d) <: RealVec
        v = zeros( typeof(d[1]) , 4 , 4 , 4 )
        for i=1:4, j=1:4, k=1:4, l=1:4
            v[i,j,k] += levicivita([i,j,k,l])*d[l]
        end
        return v
    elseif typeof.((a,b,c,d)) == (Int,Int,Int,Int) && a != 0 || b != 0 ||
        c != 0 || d != 0
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
"""
Hodge vector function

    HodgeV( U::RealVec , V::RealVec , W::RealVec )

The function `HodgeV` computes `ϵ(0,U,V,W)`` and effectively raises the
index of the resulting (dual) vector.

"""
function HodgeV( U::RealVec , V::RealVec , W::RealVec )
    v = ϵ(0,U,V,W)
    v[1] = -v[1]        # Raise index
    return v
end     #---------------------------------------------------------------
