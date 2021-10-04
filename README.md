# cereal.jl
This is a simple, Special Relativistic Location (SRL/cereal) tool for use in relativistic positioning systems. It computes the intersection of future pointing light cones from four distinct emission points which lie on a spacelike hyperplane.

The algorithm is conceptually simple:

  1. Find and perform a Lorentz transformation such that the four points *X<sub>i</sub>* (labeled with *i* ∈ {1,2,3,4}) have the same time coordinate in the transformed frame. In the transformed frame, the four points form a tetrahedron in space at some instant in time.

  2. Assuming the tetrahedron has finite volume, find the circumcenter and circumradius. The circumcenter provides the spatial coordinates of the intersection point in the new frame, and the circumradius provides the time coordinate (via time of flight).

  3. Transform back to obtain the coordinates of the intersection point in the original frame.

The computed point *P* should satisfy the four constraints *dX<sub>i</sub>*<sup>2</sup>:=*η<sub>μν</sub>*(*P<sup>μ</sup>-X<sup>μ</sup><sub>i</sub>*)(*P<sup>ν</sup>-X<sup>ν</sup><sub>i</sub>*)=0, where *μ*, *ν* are spacetime indices and *η<sub>μν</sub>* is the Minkowski metric.

To try out the code, open the Julia REPL in the directory containing ```cereal.jl``` , and run:

```julia
include("cereal.jl")
```
In the REPL, you can change directory using the command ```cd("[directory]")```, with ```[directory]``` being the directory containing ```cereal.jl```.

The following generates an 4x4 array array ```X```, the column vectors of which are the coordinates for each emission point:

```julia
X = cerealtest.epgen(Float64)
```
For higher precision calculations, one can use [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl):

```julia
using DoubleFloats
X = cerealtest.epgen(Double64)
```

Alternatively, ```X``` may be typed explicitly (by default, the elements will be interpreted as double precision [Float64] values):

```julia
X = [   0.444567   0.382671   0.422239   0.448663   ;
       -0.594917   0.338193  -0.255980  -0.662920   ;
       -0.287310   0.289673   0.345037  -0.216705   ;
       -0.693082  -0.274554  -0.381769   0.595151   ]
```

To compute the intersection point, run

```julia
P = cereal.locator(X)
```

To check the result, run the function:

```julia
cerealtest.single(1e-13,P,X)
```

The test function ```cerealtest.single(q,P,X)``` returns ```true``` if the constraints *dX<sub>i</sub>*<sup>2</sup>=0 are satisfied up a threshold value ```q```, and returns the array ```X``` if they are not.

Specifically, the test function computes the separation vectors ```V[i]=X[:,i]-P``` for each of the emission points ```X[:,i]``` and the intersection point ```P```, and compares the average of the ratio of *dX<sub>i</sub>*<sup>2</sup> (which should be zero) to the squared Euclidean norm *δ<sub>μν</sub>*(*P<sup>μ</sup>-X<sup>μ</sup><sub>i</sub>*)(*P<sup>ν</sup>-X<sup>ν</sup><sub>i</sub>*) for the separation vectors ```V[i]```. The test function returns ```true``` if the absolute value of the averaged ratio is less than ```q```.

For a more complete test, run the test function ```cerealtest.full(n,q)```, which performs the test described above for ```n``` randomly generated test cases. The following is an example which generates a million test cases with a threshold of ```1e-13```:

```julia
cerealtest.full(1e6,1e-13)
```

The default precision is Float64 (double precision), but one may increase the precision by using [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl) instead:

```julia
using DoubleFloats
n = Int64(1e6)
q = Double64(1e-26)
cerealtest.full(n,q)
```

With the increased precision, the cereal code passes this test with a much lower threshold of ```1e-26```.

To my knowledge, the particular algorithm implemented in this code does not yet appear in the literature, though it is a straightforward consequence of the existence of a spacelike configuration hyperplane discussed in Coll et al., Class.Quant.Grav. 27 (2010) 065013 (see also Fig. 3 of Coll et al., Phys. Rev. D 86, 084036 (2012)). The algorithm implemented here* is chosen for its conceptual simplicity. It should be mentioned that the formula described in Coll et al., Class.Quant.Grav. 27 (2010) 065013 has been implemented in Puchades et al., Astrophys.Space Sci. 341 (2012) 631-643. Another algorithm which accomplishes the same is described in Kostić et al., Class. Quantum Grav. 32 (2015) 215004 and Čadež et al., Advanced Concepts Team, ARIADNA final report (09/1301), European Space Agency, 2010.

*It is perhaps worth mentioning that the construction of the frame here is distinct from the concept of the central observer in Pozo, J. M, Journées 2005. Systèmes de référence spatio-temporels., Pologne (2005); AIP Conference Proceedings 841, 641 (2006), and Coll et al. Phys. Rev. D 86, 084036. The frame here differs in that no restriction is placed on the properties of the tetrahedron formed by the emission points.
