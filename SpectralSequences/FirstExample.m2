unineedsPackage "SpectralSequences"

tensor(RingMap,ChainComplex) := ChainComplex => 
opts -> (f,C) -> (
    k := min C; 
    D := chainComplex(
	if even(k) then apply(
	    drop(select(keys C.dd, 
	    	i -> instance(i,ZZ)),1), 
	    j -> f ** C.dd_j)
	else apply(
	    drop(select(keys C.dd, 
	    	i -> instance(i,ZZ)),1), 
	    j -> (-1) * (f ** C.dd_j)));
    D[-k]
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
E = prune spectralSequence E0;
E_infinity
E_2 .dd_{1,0}

viewHelp netPage
netPage(E_3,{3,0},{5,1})
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
-- Example 3: more interesting version of Example 1
restart
load "FirstExample.m2"
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
load "FirstExample.m2"
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
load "FirstExample.m2"
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
-- Example with Ext-Ext (4):
restart
load "FirstExample.m2"
k = QQ;
R = k[a,b,c, Degrees => {2,2,2}]/ideal"b2-ac";
S = k[s,t]
phi = map(S,R,{s^2,s*t,t^2})
-- N is k over the Veronese
N' = R^1/ideal"a,b,c"
-- M is k over the polynomial ring
N = S^1/ideal(s,t)
-- C is Hom_R(S,N)
C = res (pushForward(phi,S^1), LengthLimit => 7)
Hom(C,N')
C' = tensor(phi,Hom(C,N'))
D = res N
K = Hom( filteredComplex D, C')
E = prune spectralSequence K
pushForward(phi,(E_infinity)_{-1,0})
Ext^1(N',N')
E_infinity
(E_infinity)_{0,0}
---
for i to 10 list (if E_i_{1,0} != E_infinity_{1,0} then E_i_{1,0} else continue
H = Hom(C,N)
H.dd
tensor(phi,H.dd_1)
methods symbol *
C.dd
prune pushForward(phi,S^1)
prune presentation pushForward(phi,S^1)



