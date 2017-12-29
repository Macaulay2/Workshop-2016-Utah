restart

needsPackage "ToricVolumeProject"
needsPackage "ToricCompleteIntersections"

eg11 = parseKS getKreuzerSkarke 11;
eg11 = parseKS get "~/my-projects/toric-volume/cohomology/m2-example/h11-11-KS.txt";
time P1 = convexHull matrixFromString last first eg11
time P2 = polar P1
time LP1 = select(latticePoints P1, v -> v != 0)
time LP2 = select(latticePoints P2, v -> v != 0)
time for i from 1 to dim P2 list (faces(i,P2))/vertices
faces(1,P2) / vertices
latticePoints P2
faces(2,P2)


  

-- the following are very slow
elapsedTime findInteriors P -- 11 seconds
P = convexHull matrixFromString last first eg11
elapsedTime findInteriors P2 -- 10.5 seconds
P = convexHull matrixFromString last first eg11
-- each findInteriors below takes 2 - 3.6 seconds
for i from 0 to dim P - 1 do elapsedTime findInteriors(i,P)
faces(2,P)
time for f in faces(2,P) list vertices f; -- faceIndices(P,f);

loadPackage "PolymakeInterface"



restart

needsPackage "ToricCompleteIntersections"

str = ///    1    0    0    0   -1    1    0    0    0    1   -1    1   -2
         0    1    0    0    1   -1    0    1    1   -1    1   -2    0
         0    0    1    0    1   -1    0    1    0    0   -1   -2    2
         0    0    0    1   -1    1   -1   -1   -1    1    1    0   -1
         ///
M = matrixFromString str
P = convexHull M
h11OfCY P
h21OfCY P
h11OfCY polar P
h21OfCY polar P
hodgeOfCYToricDivisor(P,{0,0,0,1})

P = P1
lp = select(latticePoints P,p->p!= 0)
LP = for p in lp list flatten entries lift(p,ZZ)
for p in LP list p => hodgeOfCYToricDivisor(P,p)
hashTable oo
#LP

time intP = interiorLattice P;


M =  matrixFromString last eg11_2
time P1 = convexHull M
elapsedTime intP = interiorLattice P1;
elapsedTime interiorLattice polar P1;
time h11OfCY polar P1
time h21OfCY polar P1
time h11OfCY P1
time h21OfCY P1
time hodgeOfCYToricDivisors(P1, intP)
time h11OfCY(P1, intP)
time h21OfCY(P1, intP)
isFavorable(P, intP)

elapsedTime V = reflexiveToSimplicialToricVariety M
# rays V
ring V
dual monomialIdeal V
j = 0
for i from 0 to 4 list HH^i(V, OO V_j)

for i from 0 to #rays X - 1 list 
for i from 1 to 100 list (
    M := matrixFromString last eg11_i;
    P1 := convexHull M;
    elapsedTime intP = interiorLattice P1;
    << "example " << i << ": " << M << endl;
    << "  header: " << eg11_i_0 << endl;
    << " (h11, h12, favorable) = ";
    h11 = h11OfCY(P1, intP);
    h21 = h21OfCY(P1, intP);
    isfav = isFavorable(P1,intP);
    << (h11, h21, isfav) << endl;
    )


