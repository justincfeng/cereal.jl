# cereal.jl
This is a Special Relativistic Location (SRL/cereal) code. It computes the intersection of future pointing light cones from four distinct spacelike separated emission points.

The algorithm is simple to describe. First, find and perform a Lorentz transformation such that all points lie on the same plane in the transformed frame. The four points form a tetrahedron. Assuming the tetrahedron has finite volume, find the circumcenter and circumradius; the circumcenter provides the spatial coordinates of the intersection point, and the circumradius provides the time coordinate. Transform back to obtain the coordinates of the intersection point.
