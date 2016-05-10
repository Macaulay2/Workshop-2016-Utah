restart
debug needsPackage "NumericalAlgebraicGeometry"
R = CC[x,y,a]
PS = polySystem {x^2+y^2-1, x+a*y}
makeGateMatrix(PS,Parameters=>drop(gens R,2))  
PH := parametricSegmentHomotopy PS
a0 = 0; a1 = 1;
H = specialize (PH, transpose matrix{{a0,a1}})
s'sols = { {{0,1}},{{0,-1}} }/point
time sols = trackHomotopy(H,s'sols)
assert areEqual(sols,{{ { -.707107, .707107}, SolutionStatus => Regular }, { {.707107, -.707107}, SolutionStatus => Regular }} / point)
