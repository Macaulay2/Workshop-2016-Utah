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

--Other scratch Example
restart
needsPackage"SpectralSequences"
viewHelp"PushForward"
k=QQ
R=k[a,b,c]
S=k[s,t]
f = map(S,R,{s^2,s*t,t^2})
N = (pushFwd f)_0
cN = chainComplex(R)
cN_0 = N
M = coker vars R
C = res M

F = (filteredComplex C)** cN
E = prune spectralSequence F
E^0
E^1
E^2

FF = C ** (filteredComplex cN)
EE = prune spectralSequence FF
EE^0
EE^1
EE^2


E^infinity
EE^infinity

--
-- some other experimental things
--

restart
needsPackage"SpectralSequences"

-- possible rewrite  of pushFwd Chain Complex code 
-- What's written below might update the existing code
-- it relies on the pushFwd method from the package "PushForward.m2"

pushFwd(RingMap,ChainComplex):=o->(f,C) ->
(    pushFwdC := chainComplex(source f);
     maps := apply(spots C, i-> (i,pushFwd(f,C.dd_i)));
     for i from min C to max C do (
	 pushFwdC.dd_(maps#i_0) = maps#i_1 
	 );
    pushFwdC
    )


--possible rewrite of Change of Rings Tor code
-- What's written below might update the existing code

changeOfRingsTor = method()
changeOfRingsTor(Module,Module,RingMap) := (M,N,f) -> (
    -- f : R --> S a finite ring map, N an S module, M an R module
    F := complete res N;
    pushFwdF := pushFwd(f,F);
    G := complete res M;
    E := spectralSequence(filteredComplex(G)** pushFwdF);
    EE := spectralSequence((G) ** (filteredComplex pushFwdF));
    (E,EE) 
)


-- Try an example
k=QQ
R=k[a,b,c]
S=k[s,t]
f = map(S,R,{s^2,s*t,t^2})
kappa = coker vars S
kkappa = coker vars R
(E,EE) = changeOfRingsTor(kkappa,kappa,f)
e = prune E
ee = prune EE

e^0
e^1
e^2
e^infinity

ee^0
ee^1
ee^2
(ee^2).dd
ee^3
ee^infinity 
