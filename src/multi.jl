#-----------------------------------------------------------------------
"""
# Multiple emission set vector generator

    multivec( X::RealMtx , k::Int )

This function takes a ``m×n`` matrix `X`, and for `k```<n``,
constructs a vector of ``m×```k` matrices constructed from all choices
of `k` columns from the ``m×n`` matrix `X`.

"""     #---------------------------------------------------------------
function multivec( X::RealMtx , k::Int )  
    tpfl = typeof(X[1,1])
    d    = size(X)[1]
    np   = size(X)[2]
    Z    = [zeros(tpfl,d) for _ in 1:np]

    for i=1:np
        Z[i] = X[:,i]
    end

    w   = collect(combinations(Z,k))

    l   = length(w)
    z   = [zeros(tpfl,d,k) for _ in 1:l]
    for i=1:l, j=1:k
        z[i][:,j] = w[i][j]
    end

    return z
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Multiple emission set locator

    mlocator( X::RealMtx , locator::Function , Nbase::Int , dual::Bool, q::Real )

This computes the intersection point ``X_c`` for a large number of
emission points. Given an ``m×n`` matrix `X` of ``n`` emission
points, and a function `locator` designed to work with `Nbase```<n``
emission points, this function is designed to minimize errors and solve
the bifurcation problem when needed.

For `locator` functions designed to work with ``N = 4`` emission
points, the bifurcation problem is solved by way of a rudimentary
clustering algorithm. The clustering algorithm sorts points according to
their norms and identifies closely clustered points according to the
differences between the neighbors of the sorted points, and the tolerance
parameter `q`.

In all cases, the errors are minimized by sorting points according to
their Minkowski norms and selecting the point with the smallest norm.

"""     #---------------------------------------------------------------
function mlocator( X::RealMtx , locator::Function , Nbase::Int=5 , 
                   dual::Bool=false , q::Real=1e-14 )
    tpfl = typeof(X[1,1])

    l = size(X)

    if Nbase < 4
        print("Need at least 4 base points.")
        return zeros(tpfl,4)
    elseif dual || Nbase == 4
        if l[2] == 4
            return locator(X)[1]    # This only returns one point.
        elseif l[2] > 4
            W   = multivec(X,Nbase)
            k   = length(W)
            w   = [zeros(tpfl,4) for _ =1:2*k ]
            for a=1:k
                if Nbase==4
                    (w[a],w[k+a]) = locator( W[a] )
                else
                    w[a] = locator( W[a] )
                end
            end
    
            #   The following picks out closely clustered points
            k     = length(w)
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
    
            #   The following returns vectors closely clustered
            WS = ws[is][1:j]
    
            #   Re-sort according to smallest Minkowski norm
            δW = zeros(tpfl,j)
            for a=1:j
                for i=1:4
                    δW[a] += abs(η(X[:,i]-WS[a],X[:,i]-WS[a]))
                end
            end
            iW = sortperm(δW)
    
            return WS[iW][1]
        end
    elseif dual==false || Nbase > 4
        if l[2] <= Nbase
            return locator(X)
        elseif l[2] > Nbase
            W   = multivec(X,Nbase)
            k   = length(W)

            w   = [zeros(tpfl,4) for _ =1:2*k ]

            for a=1:k
                w[a] = locator(W[a])
            end

            #   Re-sort according to smallest Minkowski norm
            δW = zeros(tpfl,k)
            for a=1:k, i=1:4
                    δW[a] += abs(η(X[:,i]-w[a],X[:,i]-w[a]))
            end
            iW = sortperm(δW)

            return w[iW][1]
        end
    else
        print("Invalid arguments.")
        return zeros(tpfl,4)
    end
end     #---------------------------------------------------------------