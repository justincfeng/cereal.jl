var documenterSearchIndex = {"docs":
[{"location":"locatorfuncs.html#Methods-and-Locator-Functions","page":"Locator Functions","title":"Methods and Locator Functions","text":"","category":"section"},{"location":"locatorfuncs.html#Locator-function-selector","page":"Locator Functions","title":"Locator function selector","text":"","category":"section"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"cereal.locatorselect","category":"page"},{"location":"locatorfuncs.html#cereal.locatorselect","page":"Locator Functions","title":"cereal.locatorselect","text":"Locator selection function.\n\nlocatorselect( N::Int , Method::String )\n\nThe function locatorselect selects the locator function to use. The first argument N determines the minimum number of emission points that the locator function will use. N=4 or N≥5 is recommended. The second argument Method selects the formula or algorithm used to compute the intersection point X_c. The available methods are \"FHC21\", \"CFM10\", and for N≥5, \"RTC21\".\n\nBy default, the locator function assumes N=5 and Method=\"RTC21\":\n\njulia> locatorselect() == locatorselect(5,\"RTC21\")\ntrue\n\n\n\n\n\n","category":"function"},{"location":"locatorfuncs.html#Basic-method-functions","page":"Locator Functions","title":"Basic method functions","text":"","category":"section"},{"location":"locatorfuncs.html#The-RTC21-method","page":"Locator Functions","title":"The RTC21 method","text":"","category":"section"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The RTC21 formula of Ruggiero, Tartaglia, and Casalino (2021) [arxiv:2111.13423] is the preferred method, since it is the fastest and most accurate, despite requring an additional emission point. The calculation is rather staightforward; first define the matrix mathcalM and the vector B:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"mathcalM = 2\n  left\n  beginarraycccc\n    t_1-t_2    x_2-x_1    y_2-y_1    z_2-z_1 \n    t_2-t_3    x_3-x_2    y_3-y_2    z_3-z_2 \n    t_3-t_4    x_4-x_3    y_4-y_3    z_4-z_3 \n    t_4-t_5    x_5-x_4    y_5-y_4    z_5-z_4 \n  endarray\n  right ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"B = \n  left\n  beginarraycccc\n    x_2^2-x_1^2 + y_2^2-y_1^2 + z_2^2-z_1^2 + t_1^2-t_2^2 \n    x_3^2-x_2^2 + y_3^2-y_2^2 + z_3^2-z_2^2 + t_2^2-t_3^2 \n    x_4^2-x_3^2 + y_4^2-y_3^2 + z_4^2-z_3^2 + t_3^2-t_4^2 \n    x_5^2-x_4^2 + y_5^2-y_4^2 + z_5^2-z_4^2 + t_4^2-t_5^2 \n  endarray\n  right ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The intersection point X_c is then given by:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"X_c = mathcalM^-1 B ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"This method is implemented in the following function:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"cereal.locator5RTC21","category":"page"},{"location":"locatorfuncs.html#cereal.locator5RTC21","page":"Locator Functions","title":"cereal.locator5RTC21","text":"Five point locator function (RTC21)\n\nlocator5RTC21( X::RealMtx )\n\nThe function presented here implements the five point relativistic  location formula given by Ruggiero, Tartaglia, and Casalino in  Ruggiero et al., arxiv:2111.13423 (2021).\n\n\n\n\n\n","category":"function"},{"location":"locatorfuncs.html#The-FHC21-method","page":"Locator Functions","title":"The FHC21 method","text":"","category":"section"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The FHC21 method is a method developed by Feng, Hejda, and Carloni,  which shows slight improvement over the CFM10 method.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The algorithm is conceptually simple for a spacelike hyperplane:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"Find and perform a Lorentz transformation such that the four points  X_i (labeled with i{1,2,3,4}) have the same time coordinate  in the transformed frame. In the transformed frame, the four points  form a tetrahedron in space at some instant in time.\nAssuming the tetrahedron has finite volume, find the circumcenter and  circumradius. The circumcenter provides the spatial coordinates of  the intersection point in the new frame, and the circumradius  provides the time coordinate (via time of flight).\nTransform back to obtain the coordinates of the intersection point in  the original frame.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The computed point X_c should satisfy the four constraints","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"dX_I^2=η_μν(X_c^μ-X^μ_I)(X_c^ν-X^ν_I)=0","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"where μ, ν are spacetime indices and η_μν is the Minkowski  metric. ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The case of a timelike hyperplane is conceptually similar, but more intricate:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"Find and perform a Lorentz transformation such that the four points  X_I (labeled with I {1,2,3,4}) lie on a plane P_z  specified by a z coordinate in the transformed frame. In the  transformed frame, the four points lie on an elliptic hyperboloid.\nFind the vertex of the hyperboloid, and the Minkowski distance R  from the vertex to a point on the hyperboloid. There are two  intersection points. The z coordinate for the intersection  points is a distance R in the direction normal to the P_z  plane.\nTransform back to obtain the coordinates of the intersection point in  the original frame.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"There are two intersection points in the case of a timelike hyperplane; this is known as the bifurcation problem (Coll et al., Phys. Rev. D 86, 084036 (2012)). The addition of a fifth emission point will in most cases permit the selection of a single emission point.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The method described above is implemented in the following function:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"cereal.locator4FHC21","category":"page"},{"location":"locatorfuncs.html#cereal.locator4FHC21","page":"Locator Functions","title":"cereal.locator4FHC21","text":"Four point locator function (FHC21)\n\nlocator4FHC21( X::RealMtx )\n\nThis function implements the four point relativistic location algorithm of the authors. It outputs a pair of location points.\n\n\n\n\n\n","category":"function"},{"location":"locatorfuncs.html#The-CFM10-method","page":"Locator Functions","title":"The CFM10 method","text":"","category":"section"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The CFM10 method is the original four point relativistic location formula given by Coll, Ferrando, and Morales-Lladosa: ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"X_c^μ = X_4^μ + y_*^μ - fracη(y_*y_*)η(y_*χ)sqrtΔ χ^μ ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"where the following are defined:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"Δ = η(y_*χ)^2 - η(y_*y_*) η(χχ) ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"y_*^μ = frac1η(ξχ) (i_ξ H)^μ ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"H = Ω_1 h(e_2e_3) + Ω_2 h(e_3e_1) + Ω_3 h(e_1e_2)","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"Ω_i = η(e_ie_i)","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"χ^μ = η^μσ ϵ_σαβγ e_1^α e_2^β e_3^γ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"e_i^μ = X_i^μ - X_4^μ ","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"where h(UV) yields the rank-2 tensor -ϵ_αβμν U^μ U^ν (so that H is a rank-2 quantity), and ξ^μ is an arbitrary vector which we choose to satisfy η(ξχ)=1. See See Coll et al., Class.Quant.Grav. 27 (2010) 065013 and Coll et al., Phys. Rev. D 86, 084036 (2012) for the details of the derivation.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"The method describe above is implemented in the following function:","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"cereal.locator4CFM10","category":"page"},{"location":"locatorfuncs.html#cereal.locator4CFM10","page":"Locator Functions","title":"cereal.locator4CFM10","text":"Four point locator function (CFM10)\n\nlocator4CFM10( X::RealMtx )\n\nThe function presented here implements the four point relativistic location formula given by Coll, Ferrando, and Morales-Lladosa. It outputs a pair of location points.\n\n\n\n\n\n","category":"function"},{"location":"locatorfuncs.html#Multi-locator-function","page":"Locator Functions","title":"Multi locator function","text":"","category":"section"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"Given a large number of emission points, one can reduce errors and avoid the bifurcation problem. A multi locator function is implemented here, which computes the intersection point X_c. Given an mn matrix X of n emission points, and a function locator designed to work with Nbasen emission points, this function is designed to minimize errors and solve the bifurcation problem when needed.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"For locator functions designed to work with N = 4 emission points, the bifurcation problem is solved by way of a rudimentary clustering algorithm. The clustering algorithm sorts points according to their norms and identifies closely clustered points according to the differences between the neighbors of the sorted points, and the tolerance parameter q.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"In all cases, the errors are minimized by sorting points according to their Minkowski norms and selecting the point with the smallest norm.","category":"page"},{"location":"locatorfuncs.html","page":"Locator Functions","title":"Locator Functions","text":"cereal.mlocator","category":"page"},{"location":"locatorfuncs.html#cereal.mlocator","page":"Locator Functions","title":"cereal.mlocator","text":"Multiple emission set locator\n\nmlocator( X::RealMtx , locator::Function , Nbase::Int , dual::Bool, q::Real )\n\nThe function mlocator computes a single intersection point X_c for a large number of emission points X, given a locator function locator.\n\n\n\n\n\n","category":"function"},{"location":"index.html#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"index.html#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Relativistic positioning refers to the concept of establishing spacetime positions from proper time broadcasts emitted by a system of satellites. Central to relativistic positioning is the relativistic location location problem, which is the problem of finding the intersection of future pointing light cones from a collection of at least four emission points. cereal.jl contains a collection of functions for the relativistic location problem in flat spacetime. ","category":"page"},{"location":"index.html#Short-tutorial","page":"Home","title":"Short tutorial","text":"","category":"section"},{"location":"index.html#Setup","page":"Home","title":"Setup","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The cereal.jl code was written for and tested in Julia 1.6; we recommend Julia 1.6 or newer. To add the code, run the following command in the package manager for the Julia REPL:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"pkg> add https://github.com/justincfeng/cereal.jl/","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Once added, one may access the cereal module with the following command:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> using cereal","category":"page"},{"location":"index.html#Relativistic-locator","page":"Home","title":"Relativistic locator","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"To run the cereal code, one begins by generating a set of emission points with the following:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> ( X , Xtar ) = cereal.ceval.pgen(Float64,5)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The quantity Xtar is a four component vector representing the true intersection point, and X is a 45 matrix consisting set of 4 column vectors representing the coordinates of the emission points. The emission points are constructed by finding points on the past light cone of the target point Xtar.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Three different methods for finding the intersection point have been implemented, which are represented by the strings CFM10, FHC21 and RTC21. The method RTC21 (see reference below) is recommended, but requires at least five emission points. ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"To select the locator function associated with the method RTC21, use the cereal.locatorselect function, which outputs the appropriate locator function:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> locator = cereal.locatorselect(5,\"RTC21\")","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The first argument is the number of emission points; for the method RTC21, the value should be at least 5 (larger values yield functions which take additional points into consideration). Once the locator function is selected, one may feed the emission point matrix X into the locator function to obtain the intersection point Xc:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> Xc = locator(X)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The intersection point may then be compared with Xtar:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> Xc - Xtar","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"In most cases, the differences in the components should be on the order of the machine precision (10^-15 for the default floating point type Float64).","category":"page"},{"location":"index.html#Other-methods","page":"Home","title":"Other methods","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The methods CFM10, FHC21 are four point methods, the first argument of cereal.locatorselect can have a value of at least 4. ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"To try out the method CFM10, one may use the following command:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> locator4a = cereal.locatorselect(4,\"CFM10\")","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"To try out the method FHC21, use:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> locator4b = cereal.locatorselect(4,\"FHC21\")","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Since four point methods generally suffer from the bifurcation problem (see Coll et al., Phys. Rev. D 86, 084036 (2012)), these locator functions return a tuple of points. It should be mentioned that if one feeds 45 matrix X, the functions locator4a and locator4b only use the first four emission points.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> Xca = locator4a(X)\n\njulia> Xcb = locator4b(X)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"In the tuple Xca, either Xca[1] or Xca[2] should be close to the point Xtar.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"If one increases the number of emission points, then the resulting functions take additional emission points into consideration for the purpose of minimizing errors:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> locator5a = cereal.locatorselect(5,\"CFM10\")\n\njulia> locator5b = cereal.locatorselect(5,\"FHC21\")","category":"page"},{"location":"index.html#Evaluation","page":"Home","title":"Evaluation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Routines have been written to evaluate the methods more comprehensively. The function cereal.ceval.main(locator,N,q,ne) takes a locator function locator generates N sets of emission points X on the past light cone of target points Xtar, feeds each set into locator, and checks that locator yields results Xc that differ from Xtar by a factor less than a threshold value q. One may run the following:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> cereal.ceval.main(cereal.locatorselect(4,\"CFM10\"),100000,1e-6,4)\n\njulia> cereal.ceval.main(cereal.locatorselect(4,\"FHC21\"),100000,1e-6,4)\n\njulia> cereal.ceval.main(cereal.locatorselect(5,\"RTC21\"),100000,1e-9,5)\n\njulia> cereal.ceval.main(cereal.locatorselect(6,\"RTC21\"),100000,5e-13,6)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"At most, one should encounter less than 10 failures in each case.","category":"page"},{"location":"index.html#References","page":"Home","title":"References","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"CFM10: Coll, B. and Ferrando, J. J. and Morales-Lladosa, J. A., Positioning Systems in Minkowski Space-Time: from Emission to Inertial Coordinates, Class. Quant. Grav. 27, 065013 (2010)   doi:10.1088/0264-9381/27/6/065013 [arXiv:0910.2568]","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"RTC21: Ruggiero, M. L., Tartaglia, A., Casalino, L., Geometric approach to the definition of emission coordinates, (2021)   [arXiv:2111.13423]","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"FHC21: Feng, J. C., Hejda, F., Carloni, S., Relativistic location algorithm in curved spacetime, (2021)   [In preparation]","category":"page"},{"location":"evaluation.html#Evaluation","page":"Evaluation","title":"Evaluation","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The details of evaluation functions are provided here in case one is interested in performing more detailed tests of the location algorithms.","category":"page"},{"location":"evaluation.html#Main-evaluation-function","page":"Evaluation","title":"Main evaluation function","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"cereal.ceval.main","category":"page"},{"location":"evaluation.html#cereal.ceval.main","page":"Evaluation","title":"cereal.ceval.main","text":"Evaluation function.\n\nmain( locator::Function , N::Number , q::Real , k::Number ,     \n      counter::Bool , usexgen::Bool )\n\nThe function main tests the user specified locator function for N stochastically generated test cases. The results produced by the locator are compared (by way of the comp function) to the generated intersection points up to a threshold value of q. The variable k is the number of emission points to generate for each test case, and the variable usexgen replaces the function pgen with xgen for test case generation.\n\nBasic example:\n\ncereal.ceval.main(cereal.locatorselect(5,\"RTC21\"),1e5,5e-13,5)\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html#Test-case-generator-functions","page":"Evaluation","title":"Test case generator functions","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"cereal.ceval.vgenerator","category":"page"},{"location":"evaluation.html#cereal.ceval.vgenerator","page":"Evaluation","title":"cereal.ceval.vgenerator","text":"Vector generator.\n\nvgenerator( tpfl::DataType=Float64 )\n\nThe function vgenerator generates a random 3-vector of unit length. Returns a three component vector.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"cereal.ceval.nullgen","category":"page"},{"location":"evaluation.html#cereal.ceval.nullgen","page":"Evaluation","title":"cereal.ceval.nullgen","text":"Null vector generator.\n\nnullgen( tpfl::DataType=Float64 )\n\nThe function nullgen generates a random past directed null vector. Returns a four component vector.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"cereal.ceval.pgen","category":"page"},{"location":"evaluation.html#cereal.ceval.pgen","page":"Evaluation","title":"cereal.ceval.pgen","text":"Intersection and emission point generator.\n\npgen( tpfl::DataType=Float64 , N::Int=4 )\n\nThe function pgen generates a point Xc and N random emission points (in a 4N matrix X) on the past null cone of Xc. Returns a tuple (X,Xc).\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"cereal.ceval.xgen","category":"page"},{"location":"evaluation.html#cereal.ceval.xgen","page":"Evaluation","title":"cereal.ceval.xgen","text":"Restricted intersection and emission point generator.\n\nxgen( xc::Real , r1::Real , r2::Real=r1 , N::Int=4 )\n\nThe function xgen generates a point Xc=[0;xc;0;0] and N random emission points (in a 4N matrix X) on the past null cone of Xc at a radius r such that r1<r<r2. Returns a tuple (X,Xc).\n\n\n\n\n\n","category":"function"}]
}
