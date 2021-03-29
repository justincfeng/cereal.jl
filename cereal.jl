module cereal      # cereal module

using LinearAlgebra      # The only external library needed here

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

function η( V1::RealVec , V2::RealVec )        # Computes Minkowski inner product
    nv1 = length(V1)
    nv2 = length(V2)

    met = -V1[1]*V2[1]

    if nv1==nv2
        for i=2:nv1
            met += V1[i]*V2[i]
        end
        return met
    else
        print("Vectors are of a different dimension.")
    end
end  # End η

function Frame( X::RealMtx )        # Constructs spatial frame
    tpfl=typeof(X[1,1])       # Extracting floating point datatype
    s = 1       # Integer for constructing frame
    b = true    # Boolean variable - stops routine if not spacelike

    E = zeros(tpfl,4,3)

    for i=1:4
        for j=i+1:4
            V = X[:,i] - X[:,j]       # Computing difference between points
            if η( V , V ) > zero(tpfl)
                if s<=3
                    E[:,s]=V          # Frame vector
                end
                b == true
                s += 1
            else
                b == false
            end
        end
    end

    if b
        return E
    else
        println("Points not spacelike separated")
    end
end  # End Frame

function NormVecF( E::RealMtx )        # Constructs timelike normal to frame
    tpfl=typeof(E[1,1])
        NVr = [ E[2,3]*E[3,2]*E[4,1] - E[2,2]*E[3,3]*E[4,1] -
                E[2,3]*E[3,1]*E[4,2] + E[2,1]*E[3,3]*E[4,2] +
                E[2,2]*E[3,1]*E[4,3] - E[2,1]*E[3,2]*E[4,3]   ;
                #
                E[1,3]*E[3,2]*E[4,1] - E[1,2]*E[3,3]*E[4,1] -
                E[1,3]*E[3,1]*E[4,2] + E[1,1]*E[3,3]*E[4,2] +
                E[1,2]*E[3,1]*E[4,3] - E[1,1]*E[3,2]*E[4,3]   ;
                #
                E[1,2]*E[2,3]*E[4,1] - E[1,3]*E[2,2]*E[4,1] +
                E[1,3]*E[2,1]*E[4,2] - E[1,1]*E[2,3]*E[4,2] -
                E[1,2]*E[2,1]*E[4,3] + E[1,1]*E[2,2]*E[4,3]   ;
                #
                E[1,3]*E[2,2]*E[3,1] - E[1,2]*E[2,3]*E[3,1] -
                E[1,3]*E[2,1]*E[3,2] + E[1,1]*E[2,3]*E[3,2] +
                E[1,2]*E[2,1]*E[3,3] - E[1,1]*E[2,2]*E[3,3]   ]
        norm = η( NVr , NVr )
        if norm < zero(tpfl)
            if NVr[1] < zero(tpfl)
                NV = -NVr / sqrt(abs(norm))
            elseif NVr[1] > zero(tpfl)
                NV = NVr / sqrt(abs(norm))
            else
                println("Temporal component of normal vector too small.")
            end
        else
            println("Normal vector not timelike.")
        end
end  # End NormVecF

function LTM( NV::RealVec )        # Constructs Lorentz transformation matrix
    tpfl=typeof(NV[1])
        γ = NV[1]               # γ is the Lorentz factor
        δ = (γ - one(tpfl))     # A useful quantity
        NVs = NV[2]^2 + NV[3]^2 + NV[4]^2           # Square of spatial part
        if NVs > 0              # Normalization factor check
            nf = 1/sqrt(NVs)    # Normalization factor
            β = sqrt(NVs/(γ*γ)) # β is v/c
        else
            nf = zero(tpfl)     # Set normalization factor to zero
            β = zero(tpfl)      # β is v/c
        end
        vx = nf*NV[2]     # x-component of unit vector
        vy = nf*NV[3]     # y-component of unit vector
        vz = nf*NV[4]     # z-component of unit vector
    return [  γ        -γ*β*vx               -γ*β*vy               -γ*β*vz      ;
             -γ*β*vx   one(tpfl) + δ*(vx^2)  δ*vx*vy               δ*vx*vz      ;
             -γ*β*vy   δ*vy*vx               one(tpfl) + δ*(vy^2)  δ*vy*vz      ;
             -γ*β*vz   δ*vz*vx               δ*vz*vy               one(tpfl) + δ*(vz^2) ]
end  # End LTM

function IPfinder( Y::RealMtx )    # Finds intersection of light cones in adapted frame
  tpfl=typeof(Y[1,1])

      X = zeros(tpfl,3,4)
      V = zeros(tpfl,3,3)

      X[1,:] = Y[2,:]    # Spatial points for corners of tetrahedron
      X[2,:] = Y[3,:]
      X[3,:] = Y[4,:]

      V[:,1] = X[:,2] - X[:,1]    # "Frame" vectors centered on point X[:,1]
      V[:,2] = X[:,3] - X[:,1]
      V[:,3] = X[:,4] - X[:,1]

      B = [ X[:,2]⋅X[:,2] - X[:,1]⋅X[:,1] ; X[:,3]⋅X[:,3] - X[:,1]⋅X[:,1] ;
            X[:,4]⋅X[:,4] - X[:,1]⋅X[:,1] ] / 2    # Constructing B vector

      A = transpose([ V[:,1]  V[:,2]  V[:,3] ])    # Constructing A matrix

      xc = inv(A)*B    # Compute circumcenter

      rc = sqrt( (xc-X[:,1])⋅(xc-X[:,1]) )    # Compute distance to circumcenter

      tc = rc + Y[1,1]    # Compute time coordinate

  return [ tc ; xc[1] ; xc[2] ; xc[3] ]
end  # End IPfinder

function locator( X::RealMtx )    # Special relativistic locator code
  tpfl=typeof(X[1,1])
        E = Frame(X)            # Constructing spatial frame
        NV = NormVecF(E)        # Constructing normal vector to frame
        Λ = LTM(NV)             # Lorentz transformation matrix
        Y = zeros(tpfl,4,4)      # Create container
        for i=1:4
            Y[:,i] =  Λ * X[:,i]    # Lorentz transformation
        end

        XP = IPfinder(Y)        # Intersection point in adapted frame

        ΛR = inv(Λ)           # Lorentz transformation matrix
        Z = zeros(tpfl,4)    # Create container
        for i=1:4
           Z = ΛR * XP       # Lorentz transformation
        end
    return Z
end  # End locator

end  # End scope of module cereal


module cerealtest    # Test module for cereal

using LinearAlgebra      # The only external library needed here
import ..cereal       # Importing functions from cereal

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

function vgenerator( tpfl::DataType )    # Generates random 3-component vector
    θ = pi*rand(tpfl)
    ϕ = 2*pi*rand(tpfl)
    V = [   sin(θ)*cos(ϕ)   ;
            sin(θ)*sin(ϕ)   ;
            cos(θ)          ]
    return V
end  # End vgenerator

function nvgenerator( tpfl::DataType )   # Generates random 4-velocity
    v = vgenerator(tpfl)

    v = rand(tpfl)*rand(tpfl)*v

    V = [   one(tpfl)  ;
              v[1]     ;
              v[2]     ;
              v[3]     ]

    fl = one(tpfl) / sqrt(abs(cereal.η(V,V)))

    NV = fl*V

    return NV
end  # End nvgenerator

function pointgenerator( tpfl::DataType )     # Generates four random points W all with same time coordinate
    V1 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
    V2 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
    V3 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
    V4 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))

    Vt = rand(tpfl)

    W = [   Vt     Vt     Vt     Vt     ;
            V1[1]  V2[1]  V3[1]  V4[1]  ;
            V1[2]  V2[2]  V3[2]  V4[2]  ;
            V1[3]  V2[3]  V3[3]  V4[3]  ]

    return W
end  # End pointgenerator

function LT( W::RealMtx, NV::RealVec )     # Lorentz transforms points W to frame where vector NV is in t direction
    tpfl=typeof(W[1,1])

    Λ = cereal.LTM(NV)

    β = rand(tpfl)/(2*one(tpfl))
    γ = one(tpfl)/sqrt( one(tpfl) - β^2 )
    δ = ( γ - one(tpfl) )             # A useful quantity

    vx = β*NV[1]        # x-component of unit vector
    vy = β*NV[2]        # y-component of unit vector
    vz = β*NV[3]        # z-component of unit vector

    Λ = [  γ        -γ*β*vx        -γ*β*vy        -γ*β*vz       ;
           -γ*β*vx   one(tpfl) + δ*(vx^2)  δ*vx*vy        δ*vx*vz      ;
           -γ*β*vy   δ*vy*vx        one(tpfl) + δ*(vy^2)  δ*vy*vz      ;
           -γ*β*vz   δ*vz*vx        δ*vz*vy        one(tpfl) + δ*(vz^2) ]

    Z = zeros(tpfl,4,4)

    for i=1:4
        Z[:,i] = Λ * W[:,i] # Lorentz transformation
    end

    return Z
end  # End LT

function epgen( tpfl::DataType )     # Randomly generates a set of four emission points
    W  = pointgenerator(tpfl)
    NV = nvgenerator(tpfl)
    return LT(W,NV)
end  # End epgen

function single( q::Real , P::RealVec, X::RealMtx )     # Checks if separation vectors are null
    tpfl=typeof(q)

    Y = P

    dX1 = Y - X[:,1]
    dX2 = Y - X[:,2]
    dX3 = Y - X[:,3]
    dX4 = Y - X[:,4]

    n1 = cereal.η(dX1,dX1)/dot(dX1,dX1)
    n2 = cereal.η(dX2,dX2)/dot(dX2,dX2)
    n3 = cereal.η(dX3,dX3)/dot(dX3,dX3)
    n4 = cereal.η(dX4,dX4)/dot(dX4,dX4)

    VT = [  cereal.η(dX1,dX1)/dot(dX1,dX1)  ;
            cereal.η(dX2,dX2)/dot(dX2,dX2)  ;
            cereal.η(dX3,dX3)/dot(dX3,dX3)  ;
            cereal.η(dX4,dX4)/dot(dX4,dX4)  ]

    av = ( VT[1] + VT[2] + VT[3] + VT[4] ) / (4*one(tpfl))

    if av < q
        return true
    else
        return X
    end
end  # End ctesti

function full( iters::Number, q::Real )     # Main test function--q determines floating point datatype
    tpfl=typeof(q)
    lb=false
    mi = zero(Int64)
    for i=1:Int64(iters)
        a = epgen(tpfl)
        b = single(q,cereal.locator(a),a)
        if b != true
            print(b,"\n \n")
            lb = true
            mi += 1
        end
    end
    if lb
        print("Test ended with ",mi," failed cases \n")
    else
        print("Test ended with zero failed cases \n")
    end
end  # End full

end  # End scope of module cerealtest
