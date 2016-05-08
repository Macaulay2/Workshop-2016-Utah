needsPackage "SpectralSequences"

tensor(RingMap,ChainComplex) := ChainComplex => 
opts -> (f,C) -> (
    differentials := apply(select(keys C.dd, 
	    i -> instance(i,ZZ)), j -> f ** C.dd_j);
    k := min C; 
    D := chainComplex(differentials);
    if even(k) then D[-k] else (-1)*D[-k]
    )
end

restart
load "FirstExample.m2"
-- Example 1
-- Resolves k over a polynomial ring and a quotient of it
k = QQ;
R = k[x];
S = R/ideal"x2";
N = S^1/ideal"x";
M = R^1/R_0;
C = res M;
C' = C ** S;
D = res(N,LengthLimit => 10);
E0 = C' ** (filteredComplex D);
E = prune spectralSequence E0
E_infinity
-- Example 2
restart
load "FirstExample.m2"
k = QQ;
R = k[a,b,c]/ideal"b2-ac";
S = k[s,t]
phi = map(S,R,{s^2,s*t,t^2})
M = R^1/ideal"a,b,c"
N = S^1/ideal(s,t)
C = res (M,LengthLimit => 10)
C' = tensor(phi,C)
D = res(N,LengthLimit => 10);
E0 = C' ** (filteredComplex D);
E = prune spectralSequence E0
E^infinity