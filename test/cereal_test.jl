#-----------------------------------------------------------------------
#   MINKOWSKI PRODUCT TEST
#-----------------------------------------------------------------------

tv  = [1.;0.;0.;0.]
xv  = [0.;1.;0.;0.]
yv  = [0.;0.;1.;0.]
zv  = [0.;0.;0.;1.]

@test cereal.η(tv,tv) ≈ -1
@test cereal.η(xv,xv) ≈ 1
@test cereal.η(yv,yv) ≈ 1
@test cereal.η(zv,zv) ≈ 1

@test cereal.η(tv,xv) ≈ 0
@test cereal.η(tv,yv) ≈ 0
@test cereal.η(tv,zv) ≈ 0
@test cereal.η(xv,yv) ≈ 0
@test cereal.η(xv,zv) ≈ 0
@test cereal.η(yv,zv) ≈ 0

k1  = π^2*[1.;1.;0.;0.]
k2  = π^2*[1.;0.;1.;0.]
k3  = π^2*[1.;0.;0.;1.]

kn1 = π^2*[1.;-1.;0.;0.]
kn2 = π^2*[1.;0.;-1.;0.]
kn3 = π^2*[1.;0.;0.;-1.]

@test cereal.η(k1,k1) ≈ 0
@test cereal.η(k2,k2) ≈ 0
@test cereal.η(k3,k3) ≈ 0

@test cereal.η(kn1,kn1) ≈ 0
@test cereal.η(kn2,kn2) ≈ 0
@test cereal.η(kn3,kn3) ≈ 0

#-----------------------------------------------------------------------
#   MINKOWSKI NORM TEST
#-----------------------------------------------------------------------

@test cereal.mnorm( tv ) ≈ 1
@test cereal.mnorm( xv ) ≈ 1
@test cereal.mnorm( yv ) ≈ 1
@test cereal.mnorm( zv ) ≈ 1

@test cereal.mnorm( k1 ) ≈ 0
@test cereal.mnorm( k2 ) ≈ 0
@test cereal.mnorm( k3 ) ≈ 0

@test cereal.mnorm( kn1 ) ≈ 0
@test cereal.mnorm( kn2 ) ≈ 0
@test cereal.mnorm( kn3 ) ≈ 0

#-----------------------------------------------------------------------
#   LEVI-CIVITA TESTS
#-----------------------------------------------------------------------

@test cereal.ϵ(1,2,3,4) ≈ 1
@test cereal.ϵ(4,3,2,1) ≈ 1
@test cereal.ϵ(1,4,2,3) ≈ 1
@test cereal.ϵ(3,1,2,4) ≈ 1

@test cereal.ϵ(1,4,3,2) ≈ -1
@test cereal.ϵ(2,1,3,4) ≈ -1
@test cereal.ϵ(1,3,2,4) ≈ -1
@test cereal.ϵ(3,2,1,4) ≈ -1

@test cereal.ϵ(1,1,3,4) ≈ 0
@test cereal.ϵ(1,2,1,2) ≈ 0
@test cereal.ϵ(4,2,3,4) ≈ 0
@test cereal.ϵ(1,3,4,4) ≈ 0

@test cereal.ϵ(tv,xv,yv,zv) ≈ 1
@test cereal.ϵ(zv,yv,xv,tv) ≈ 1
@test cereal.ϵ(tv,zv,xv,yv) ≈ 1
@test cereal.ϵ(yv,tv,xv,zv) ≈ 1

@test cereal.ϵ(tv,zv,yv,xv) ≈ -1
@test cereal.ϵ(xv,tv,yv,zv) ≈ -1
@test cereal.ϵ(tv,yv,xv,zv) ≈ -1
@test cereal.ϵ(yv,xv,tv,zv) ≈ -1

@test cereal.ϵ(tv,tv,yv,zv) ≈ 0
@test cereal.ϵ(tv,xv,tv,xv) ≈ 0
@test cereal.ϵ(zv,xv,yv,zv) ≈ 0
@test cereal.ϵ(tv,yv,zv,zv) ≈ 0

#-----------------------------------------------------------------------
#   MULTIVEC
#-----------------------------------------------------------------------

k   = 4
X   = rand(Float64,4,6)

w1  = cereal.multivec(X,k,true)
w2  = cereal.multivec(X,k,false)

@test length(w1) == binomial(6,4)
@test length(w2) == binomial(6,4)

@test length(w1[1]) == 4
@test size(w2[1]) == (4,4)

#-----------------------------------------------------------------------
#   HODGE
#-----------------------------------------------------------------------

U   = rand(Float64,4)
V   = rand(Float64,4)
W   = rand(Float64,4)

Z0  = cereal.HodgeV(xv,yv,zv)
Z   = cereal.HodgeV(U,V,W)
Y1  = rand(Float64)*U + rand(Float64)*V + rand(Float64)*W

ω0  = cereal.Hodge2(yv,zv)
ω   = cereal.Hodge2(V,W)
Y2  = rand(Float64)*V + rand(Float64)*W
Y3  = rand(Float64)*V + rand(Float64)*W

Qv1 = 0
Qv0 = 0

for i=1:4, j=1:4
    global Qv1  += ω[i,j]*Y2[i]*Y3[j]
    global Qv0  += ω0[i,j]*tv[i]*xv[j]
end

@test cereal.η(Y1,Z)    ≈ 0   atol=1e-13
@test cereal.η(tv,Z0)   ≈ 1   atol=1e-13
@test Qv1     ≈ 0     atol=1e-13
@test Qv0    ≈ -1    atol=1e-13

#-----------------------------------------------------------------------
#   LORENTZ TRANSFORMATION TEST
#-----------------------------------------------------------------------

NV  = [1,0,0,0] + rand(4)/20
NV  = NV/cereal.mnorm(NV)
Λ   = cereal.LTM(NV)

@test Λ*NV ≈ [1,0,0,0] atol=1e-13

#-----------------------------------------------------------------------
#   Z ROTATION TEST
#-----------------------------------------------------------------------

mv    = rand(Float64,4)
mv[1] = zero(Float64)
mv    = mv/cereal.mnorm(mv)
Mz    = cereal.MRz(mv)

@test Mz*mv ≈ [0,0,0,1] atol=1e-13

#-----------------------------------------------------------------------
#   LROT TEST
#-----------------------------------------------------------------------

X       = rand(Float64,4,4)

Xprime  = cereal.Lrot(X,Λ)

Xp0 = zeros(Float64,4,4)

for i=1:4
    Xp0[:,i] = Xprime[:,i] - Λ*X[:,i]
end

@test Xp0 ≈ zeros(Float64,4,4)  atol=1e-13

#-----------------------------------------------------------------------
#   NORMFLIPS TEST
#-----------------------------------------------------------------------

vx  = rand(Float64)
vy  = rand(Float64)
vz  = rand(Float64)

vn  = norm([vx;vy;vz])
vx  = vx/vn
vy  = vy/vn
vz  = vz/vn

vt  = rand(Float64)
o   = zero(Float64)

Vsl = [ vt ; vx ; vy ; vz ]

Vfl = [ one(Float64) ; vt*vx ; vt*vy ; vt*vz ]/cereal.mnorm(Vsl)

@test cereal.NormflipS( Vsl ) ≈ Vfl  atol=1e-13

#-----------------------------------------------------------------------
#   IPFINDER TEST
#-----------------------------------------------------------------------

Xc = rand(Float64,4)

ne = 4
ID4S = idg(ne,Xc)

@test cereal.IPFinderS(ID4S[1]) ≈ Xc atol=1e-13

ID4T = idg(ne,Xc,false)

if      norm(cereal.IPFinderT(ID4T[1])[1]-Xc) < 
        norm(cereal.IPFinderT(ID4T[1])[2]-Xc)
            @test cereal.IPFinderT(ID4T[1])[1] ≈ Xc atol=1e-6
elseif  norm(cereal.IPFinderT(ID4T[1])[2]-Xc) < 
        norm(cereal.IPFinderT(ID4T[1])[1]-Xc)
            @test cereal.IPFinderT(ID4T[1])[2] ≈ Xc atol=1e-6
end

#-----------------------------------------------------------------------
#   LOCATOR TEST
#-----------------------------------------------------------------------

q   = 1e-6
Xp  = ceval.pgen(Float64)
P   = cereal.slocator(Float64.(Xp[1]),true)
Sda = ceval.compdirect(q,P[1],Xp[2])
Sdb = ceval.compdirect(q,P[2],Xp[2])

@test Sda[1] == true || Sdb[1] == true

Xp  = ceval.pgen(Float64)
P   = cereal.mlocator(Float64.(Xp[1]),q,true)
Sd  = ceval.compdirect(q,P[1],Xp[2])

@test Sd[1] == true
