# Locator Functions

## Locator function selector

```@docs
cereal.locatorselect
```

## Basic method functions

### The `RTC21` method

The `RTC21` method is the preferred method.

```@docs
cereal.locator5RTC21
```

### The `CFM10` method

The `CFM10` method is the original four point relativistic location
formula given by Coll, Ferrando, and Morales-Lladosa: 

`` X_c^μ = X_4^μ + y_*^μ - \\frac{η(y_*,y_*)}{η(y_*,χ)±√|Δ|} χ^μ ``,

where the following are defined:

`` Δ := η(y_*,χ)^2 - η(y_*,y_*) η(χ,χ) ``

`` y_*^μ := (η(ξ,χ))^{-1} (i_ξ H)^μ ``

`` H := Ω_1 h(e_2,e_3) + Ω_2 h(e_3,e_1) + Ω_3 h(e_1,e_2) ``

`` Ω_i := η(e_i,e_i) ``

`` χ^μ = η^{μσ} ϵ_{σαβγ} e_1^α e_2^β e_3^γ``

`` e_i^μ = X_i^μ - X_4^μ ``

where ``h(U,V)`` yields the rank-2 tensor ``-ϵ_{αβμν} U^μ U^ν`` (so that
``H`` is a rank-2 quantity), and ``ξ^μ`` is an arbitrary vector which we
choose to satisfy ``η(ξ,χ)=1``. See See Coll et al., Class.Quant.Grav.
27 (2010) 065013 and Coll et al., Phys. Rev. D 86, 084036 (2012) for the
details of the derivation

```@docs
cereal.locator4CFM10
```

### The `FHC21` method

The `FHC21` method is a method developed by the authors, which shows
slight improvement over the `CFM10` method.

The algorithm is conceptually simple for a spacelike hyperplane:

1. Find and perform a Lorentz transformation such that the four points
    ``X_i`` (labeled with ``i∈``{1,2,3,4}) have the same time coordinate
    in the transformed frame. In the transformed frame, the four points
    form a tetrahedron in space at some instant in time.

2. Assuming the tetrahedron has finite volume, find the circumcenter and
    circumradius. The circumcenter provides the spatial coordinates of
    the intersection point in the new frame, and the circumradius
    provides the time coordinate (via time of flight).

3. Transform back to obtain the coordinates of the intersection point in
    the original frame.

The computed point ``X_c`` should satisfy the four constraints
``dX_I^2:=η_{μν}(X_c^μ-X^μ_I)(X_c^ν-X^ν_I)=0``, where ``μ``, ``ν`` are
spacetime indices and ``η_{μν}`` is the Minkowski metric. 

The case of a timelike hyperplane is conceptually similar, but more
intricate:

1. Find and perform a Lorentz transformation such that the four points
    ``X_I`` (labeled with ``I∈`` {1,2,3,4}) lie on a plane ``P_z``
    specified by a ``z`` coordinate in the transformed frame. In the
    transformed frame, the four points lie on an elliptic hyperboloid.

2. Find the vertex of the hyperboloid, and the Minkowski distance ``R``
    from the vertex to a point on the hyperboloid. There are two
    intersection points. The ``z`` coordinate for the intersection
    points is a distance ``R`` in the direction normal to the ``P_z``
    plane.

3. Transform back to obtain the coordinates of the intersection point in
    the original frame.

There are two intersection points in the case of a timelike hyperplane;
this is known as the bifurcation problem (Coll et al., Phys. Rev. D 86,
084036 (2012)). The addition of a fifth emission point will in most
cases permit the selection of a single emission point.

```@docs
cereal.locator4FHC21
```

## Multi locator function

Given a large number of emission points, one can reduce errors and avoid
the bifurcation problem. A multi locator function is implemented here,
which computes the intersection point ``X_c``. Given an ``m×n`` matrix
`X` of ``n`` emission points, and a function `locator` designed to work
with `Nbase```<n`` emission points, this function is designed to
minimize errors and solve the bifurcation problem when needed.

For `locator` functions designed to work with ``N = 4`` emission
points, the bifurcation problem is solved by way of a rudimentary
clustering algorithm. The clustering algorithm sorts points according to
their norms and identifies closely clustered points according to the
differences between the neighbors of the sorted points, and the tolerance
parameter `q`.

In all cases, the errors are minimized by sorting points according to
their Minkowski norms and selecting the point with the smallest norm.

```@docs
cereal.mlocator
```
