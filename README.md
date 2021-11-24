# cereal.jl
This is a simple, Special Relativistic Location (SRL/cereal) tool for use in relativistic positioning systems. It computes the intersection of future pointing light cones from at least four distinct emission points which lie on a spacelike or timelike hyperplane.

## Algorithm

The algorithm is conceptually simple for a spacelike hyperplane:

  1. Find and perform a Lorentz transformation such that the four points *X<sub>i</sub>* (labeled with *i* ∈ {1,2,3,4}) have the same time coordinate in the transformed frame. In the transformed frame, the four points form a tetrahedron in space at some instant in time.

  2. Assuming the tetrahedron has finite volume, find the circumcenter and circumradius. The circumcenter provides the spatial coordinates of the intersection point in the new frame, and the circumradius provides the time coordinate (via time of flight).

  3. Transform back to obtain the coordinates of the intersection point in the original frame.

The computed point *P* should satisfy the four constraints *dX<sub>i</sub>*<sup>2</sup>:=*η<sub>μν</sub>*(*P<sup>μ</sup>-X<sup>μ</sup><sub>i</sub>*)(*P<sup>ν</sup>-X<sup>ν</sup><sub>i</sub>*)=0, where *μ*, *ν* are spacetime indices and *η<sub>μν</sub>* is the Minkowski metric. 

The case of a timelike hyperplane is conceptually similar, but more intricate:

  1. Find and perform a Lorentz transformation such that the four points *X<sub>i</sub>* (labeled with *i* ∈ {1,2,3,4}) lie on a plane *P<sub>z</sub>* specified by a *z* coordinate in the transformed frame. In the transformed frame, the four points lie on an elliptic hyperboloid.

  2. Find the vertex of the hyperboloid, and the Minkowski distance *R* from the vertex to a point on the hyperboloid. There are two intersection points. The *z* coordinate for the intersection points is a distance *R* in the direction normal to the *P<sub>z</sub>* plane.

  3. Transform back to obtain the coordinates of the intersection point in the original frame.

There are two intersection points in the case of a timelike hyperplane; this is known as the bifurcation problem (Coll et al., Phys. Rev. D 86, 084036 (2012)). The addition of a fifth emission point will in most cases permit the selection of a single emission point.

## Example and tests

### Four emission points

To try out the code, open the Julia REPL in the directory containing ```cereal.jl``` , and run:

```julia
include("cereal.jl")
```
In the REPL, you can change directory using the command ```cd("[directory]")```, with ```[directory]``` being the directory containing ```cereal.jl```.

The following generates an 4x4 array array ```X```, the column vectors of which are the coordinates for each emission point, and the intersection point ```Xc```

```julia
(X,Xc) = ceval.pgen(Float64)
```
For higher precision calculations, one can use [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl):

```julia
using DoubleFloats
(X,Xc) = ceval.pgen(Double64)
```

Alternatively, ```X``` may be typed explicitly (by default, the elements will be interpreted as double precision [Float64] values):

```julia
X = [   0.444567   0.382671   0.422239   0.448663   ;
       -0.594917   0.338193  -0.255980  -0.662920   ;
       -0.287310   0.289673   0.345037  -0.216705   ;
       -0.693082  -0.274554  -0.381769   0.595151   ]
```

To compute the intersection points (the location code returns two), run

```julia
P = cereal.slocator(X)
```

To check the result for each point ```P[1]``` and ```P[2]```, run the function:

```julia
q = 1e-7
ceval.compdirect(q,P[1],Xc)
ceval.compdirect(q,P[2],Xc)
```

The first element of the output for the test function ```cerealtest.single(q,P,Xc)``` is ```true``` if the L1 norm of ```P-Xc``` (rescaled according to the L1 norm of ```Xc```) is smaller than a threshold value ```q```, and ```false``` otherwise. THe second element is the quantity being compared ```norm(P-Xc,1)/norm(Xc,1)```.

For a more complete test, run the test function ```cerealtest.full(n,q)```, which performs the test described above for ```n``` randomly generated test cases. The following is an example which generates a million test cases with a threshold of ```1e-6```:

```julia
n = Int64(1e6)
q = 1e-6
ceval.full(n,q)
```

The default precision is Float64 (double precision), and one typically has 1-3 failed cases with the specified threshold. One may increase the precision by using [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl) instead:

```julia
using DoubleFloats
n = Int64(1e6)
q = Double64(1e-18)
ceval.full(n,q)
```

With the increased precision, the cereal code passes this test with a much lower threshold of ```1e-18```.

### Five or more emission points

One may generate more five emission points with the following:

```julia
using DoubleFloats
N = 5
(X,Xc) = ceval.pgen(Double64,N)
```

The function ```mlocator``` computes a single intersection point in the first element:

```julia
q = Double64(1e-20)
p = cereal.mlocator(X,q)[1]
ceval.compdirect(q,p,Xc)
```

For the full test, one may perform the following test:
```julia
using DoubleFloats
N = 5
n = Int64(1e6)
q = Double64(1e-20)
ceval.fullmulti(n,q,N)
```

## Literature

To our knowledge, the particular algorithms implemented in this code do not yet appear in the literature, though they are a straightforward consequence of the existence of a spacelike or timelike configuration hyperplane discussed in Coll et al., Class.Quant.Grav. 27 (2010) 065013 (see also Fig. 3 of Coll et al., Phys. Rev. D 86, 084036 (2012)). The algorithm implemented here is chosen for its conceptual simplicity, and for its relative efficiency. It should be mentioned that the formula described in Coll et al., Class.Quant.Grav. 27 (2010) 065013 has been implemented in a code described in Puchades et al., Astrophys.Space Sci. 341 (2012) 631-643 (we also include in this repository an implementation in the file srl.jl, which is structured similarly to cereal.jl). Another algorithm which accomplishes the same is described in Kostić et al., Class. Quantum Grav. 32 (2015) 215004 and Čadež et al., Advanced Concepts Team, ARIADNA final report (09/1301), European Space Agency, 2010. Further details will be provided in an upcoming paper.
