needsPackage "SpectralSequences"
end

restart
load "SecondExample.m2"
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
E = prune spectralSequence E0;
E_infinity
E_2 .dd_{1,0}

-- Example 2
restart
load "SecondExample.m2"
k = QQ;
R = k[a,b,c]/ideal"b2-ac";
S = k[s,t]
phi = map(S,R,{s^2,s*t,t^2})
M = R^1/ideal"a,b,c"
N = S^1/ideal(s,t)
C = res (M,LengthLimit => 10);
C' = tensor(phi,C)
D = res(N,LengthLimit => 10);
E0 = C' ** (filteredComplex D);
ahE = prune spectralSequence E0
E^infinity
-- Example 3: more interesting version of Example 1
restart
load "SecondExample.m2"
k = QQ;
S = k[x];
M = S^1/ideal"x3";
R = S/ideal"x7";
N = R^1/ideal"x2";
C = res M;
C' = C ** S;
D = res(N,LengthLimit => 20);
E0 = C' ** (filteredComplex D);
E = prune spectralSequence E0
E_infinity
-- Example 4: two variable version of Example 1
restart
load "SecondExample.m2"
k = QQ;
S = k[x,y];
M = S^1/ideal"x,y";
R = S/ideal"x2,y2";
N = R^1/ideal"x,y";
C = res M;
C' = C ** S;
D = res(N,LengthLimit => 20);
E0 = C' ** (filteredComplex D);
E = prune spectralSequence E0
E_infinity
-- Example 5: three variable version of Example 1
restart
load "SecondExample.m2"
k = QQ;
S = k[x,y,z];
M = S^1/ideal"x,y,z";
R = S/ideal"x2,y2,z2";
N = R^1/ideal"x,y,z";
C = res M;
C' = C ** R;
D = res(N,LengthLimit => 7);
E0 = C' ** (filteredComplex D);
E = prune spectralSequence E0;
E_infinity
