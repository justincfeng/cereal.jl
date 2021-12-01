# Home

## Introduction

Relativistic positioning refers to the concept of establishing spacetime
positions from proper time broadcasts emitted by a system of satellites.
Central to relativistic positioning is the relativistic location
location problem, which is the problem of finding the intersection of
future pointing light cones from a collection of at least four emission
points. `cereal.jl` contains a collection of functions for the
relativistic location problem in flat spacetime. 

## Short tutorial

To run the cereal code, one begins by generating a set of emission
points with the following:

    ( X , Xtar ) = cereal.ceval.pgen(Float64,5)

The quantity `Xtar` is a four component vector representing the true
intersection point, and `X` is a ``4×5`` matrix consisting set of `4`
column vectors representing the coordinates of the emission points. The
emission points are constructed by finding points on the past light cone
of the target point `Xtar`.

Three different methods for finding the intersection point have been implemented, which are represented by the strings `CFM10`, `FHC21` and `RTC21`. The method `RTC21` (see reference below) is recommended, but requires at least five emission points. The location function is defined by the following:

    locator = cereal.locatorfunc(5,"RTC21")

## References

`CFM10`: Coll, B. and Ferrando, J. J. and Morales-Lladosa, J. A., *Positioning Systems in Minkowski Space-Time: from Emission to Inertial Coordinates*, Class. Quant. Grav. **27**, 065013 (2010)  
doi:[10.1088/0264-9381/27/6/065013](https://doi.org/10.1088/0264-9381/27/6/065013) [\[arXiv:0910.2568\]](https://arxiv.org/abs/0910.2568)

`RTC21`: Ruggiero, M. L., Tartaglia, A., Casalino, L., *Geometric approach to the definition of emission coordinates*, (2021)  
[\[arXiv:2111.13423\]](https://arxiv.org/abs/2111.13423)

`FHC21`: Feng, J. C., Hejda, F., Carloni, S., *Relativistic location algorithm in curved spacetime*, (2021)  
\[In preparation\]

## Acknowledgements


