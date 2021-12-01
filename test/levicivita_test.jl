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

tv  = [1.;0.;0.;0.]
xv  = [0.;1.;0.;0.]
yv  = [0.;0.;1.;0.]
zv  = [0.;0.;0.;1.]

@test cereal.ϵ(0,0,yv,xv)[1,4] ≈ -1
@test cereal.ϵ(0,0,yv,xv)[4,1] ≈ 1

@test cereal.ϵ(0,0,tv,yv)[2,4] ≈ -1
@test cereal.ϵ(0,0,tv,yv)[4,2] ≈ 1
@test cereal.ϵ(0,0,tv,xv)[4,3] ≈ -1
@test cereal.ϵ(0,0,tv,xv)[3,4] ≈ 1
@test cereal.ϵ(0,0,xv,yv)[1,4] ≈ 1
@test cereal.ϵ(0,0,xv,yv)[4,1] ≈ -1
@test cereal.ϵ(0,0,xv,zv)[1,3] ≈ -1
@test cereal.ϵ(0,0,xv,zv)[3,1] ≈ 1
@test cereal.ϵ(0,0,yv,zv)[1,2] ≈ 1
@test cereal.ϵ(0,0,yv,zv)[2,1] ≈ -1

@test cereal.ϵ(0,0,tv,tv) ≈ zeros(Float64,4,4)
@test cereal.ϵ(0,0,xv,xv) ≈ zeros(Float64,4,4)
@test cereal.ϵ(0,0,yv,yv) ≈ zeros(Float64,4,4)
@test cereal.ϵ(0,0,zv,zv) ≈ zeros(Float64,4,4)

@test cereal.ϵ(0,xv,yv,zv) ≈  tv
@test cereal.ϵ(0,tv,yv,zv) ≈ -xv
@test cereal.ϵ(0,xv,tv,zv) ≈ -yv
@test cereal.ϵ(0,xv,yv,tv) ≈ -zv
@test cereal.ϵ(0,tv,zv,yv) ≈ xv
@test cereal.ϵ(0,zv,tv,xv) ≈ yv
@test cereal.ϵ(0,yv,xv,tv) ≈ zv

@test cereal.ϵ(0,tv,tv,xv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,tv,tv,yv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,tv,tv,zv) ≈ zeros(Float64,4)

@test cereal.ϵ(0,tv,xv,xv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,yv,xv,xv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,zv,xv,xv) ≈ zeros(Float64,4)

@test cereal.ϵ(0,yv,tv,yv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,yv,xv,yv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,yv,zv,yv) ≈ zeros(Float64,4)

@test cereal.ϵ(0,zv,zv,zv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,zv,zv,zv) ≈ zeros(Float64,4)
@test cereal.ϵ(0,zv,zv,zv) ≈ zeros(Float64,4)

#-----------------------------------------------------------------------
#   HODGE
#-----------------------------------------------------------------------

@test cereal.HodgeV(xv,yv,zv) ≈ -tv
@test cereal.HodgeV(tv,yv,zv) ≈ -xv
@test cereal.HodgeV(xv,tv,zv) ≈ -yv
@test cereal.HodgeV(xv,yv,tv) ≈ -zv
@test cereal.HodgeV(tv,zv,yv) ≈ xv
@test cereal.HodgeV(zv,tv,xv) ≈ yv
@test cereal.HodgeV(yv,xv,tv) ≈ zv

@test cereal.HodgeV(tv,tv,xv) ≈ zeros(Float64,4)
@test cereal.HodgeV(tv,tv,yv) ≈ zeros(Float64,4)
@test cereal.HodgeV(tv,tv,zv) ≈ zeros(Float64,4)

@test cereal.HodgeV(tv,xv,xv) ≈ zeros(Float64,4)
@test cereal.HodgeV(yv,xv,xv) ≈ zeros(Float64,4)
@test cereal.HodgeV(zv,xv,xv) ≈ zeros(Float64,4)

@test cereal.HodgeV(yv,tv,yv) ≈ zeros(Float64,4)
@test cereal.HodgeV(yv,xv,yv) ≈ zeros(Float64,4)
@test cereal.HodgeV(yv,zv,yv) ≈ zeros(Float64,4)

@test cereal.HodgeV(zv,zv,zv) ≈ zeros(Float64,4)
@test cereal.HodgeV(zv,zv,zv) ≈ zeros(Float64,4)
@test cereal.HodgeV(zv,zv,zv) ≈ zeros(Float64,4)
