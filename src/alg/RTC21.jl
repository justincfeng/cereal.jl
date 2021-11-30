#-----------------------------------------------------------------------
#   RTC21
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
# Five point locator function (RTC21)

    locator5RTC21( X::RealMtx )

The function presented here implements the five point relativistic 
location formula given by Ruggiero, Tartaglia, and Casalino in 
Ruggiero et al., arxiv:2111.13423 (2021).

"""     #---------------------------------------------------------------
function locator5RTC21( X::RealMtx )   
    tpfl=typeof(X[1,1])
    l = size(X)

    if l[2] >= 5 || l[1] != 4
        (t1,t2,t3,t4,t5) = (X[1,1],X[1,2],X[1,3],X[1,4],X[1,5])
        (x1,x2,x3,x4,x5) = (X[2,1],X[2,2],X[2,3],X[2,4],X[2,5])
        (y1,y2,y3,y4,y5) = (X[3,1],X[3,2],X[3,3],X[3,4],X[3,5])
        (z1,z2,z3,z4,z5) = (X[4,1],X[4,2],X[4,3],X[4,4],X[4,5])

        C  = [  t1-t2   x2-x1   y2-y1   z2-z1    ;
                t2-t3   x3-x2   y3-y2   z3-z2    ;
                t3-t4   x4-x3   y4-y3   z4-z3    ;
                t4-t5   x5-x4   y5-y4   z5-z4    ]

        B  = [  x2^2 - x1^2 + y2^2 - y1^2 + z2^2 - z1^2 + t1^2 - t2^2  ;
                x3^2 - x2^2 + y3^2 - y2^2 + z3^2 - z2^2 + t2^2 - t3^2  ;
                x4^2 - x3^2 + y4^2 - y3^2 + z4^2 - z3^2 + t3^2 - t4^2  ;
                x5^2 - x4^2 + y5^2 - y4^2 + z5^2 - z4^2 + t4^2 - t5^2  ]

        return inv(2 .* C) * B
    else
        return zeros(tpfl,4)
    end
end      #---------------------------------------------------------------
