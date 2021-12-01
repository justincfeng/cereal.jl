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
#   MINKOWSKI COMPONENTS TEST
#-----------------------------------------------------------------------

@test det(cereal.ημν())    ≈ -1
@test tr(cereal.ημν())     ≈ 2
