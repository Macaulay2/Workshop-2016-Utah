-- in polymake
--   rows are points
--   points are given in homogenized form, 1st coordinate is the homogenizing variable.
--   

loadPackage "ToricVolumeProject"
--loadPackage "PolyhedralObjects"
--loadPackage "PolymakeInterface"

findPolarDual = (f, P, P2) -> (
    -- given a face 'f' of P2
    -- this returns the list of indices (with respect to 'vertices P') which give
    -- the face dual to f.
    positions(entries (transpose (vertices P) * (vertices P2)_f), x -> all(x, x1 -> x1 == -1))
    )
faces Polyhedron := (P) -> (
    v := entries transpose vertices P;
    -- careful! v#i is currently a vector over QQ.
    -- indexing using a list of ZZ elements won't give the same element.
    V := hashTable for i from 0 to #v-1 list v#i => i;
    hashTable for i from 1 to dim P list 
      (dim P)-i => for f in faces(i,P) list (
        for p in entries transpose vertices f list V#p
        )
    )

findInteriors = (facedim,P) -> (
   -- consider faces of P.  Return a hash table of
   -- (vertices,vertices of dual face) => {interior lattice points, dual interior lattice points}
   -- over all faces of P of dimension facedim.
   -- the indices mathc the order in 'vertices P' and 'vertices P2'
   P2 := polar P;
   v := entries transpose vertices P;
   V := hashTable for i from 0 to #v-1 list v#i => i;
   hashTable for f in faces(dim P - facedim,P) list (
       L := interiorLatticePoints f;
       verticesof := for p in entries transpose vertices f list V#p;
       dualface := findPolarDual(verticesof,P2,P);
       (verticesof, dualface) => {L, interiorLatticePoints(convexHull (vertices P2)_dualface)}
       )
   )

-- Batyrev formula   
h11 = method()
h11 Polyhedron := (P) -> (
    interior3 := findInteriors(0,P);
    interior2 := findInteriors(1,P);
    contrib3 := (values interior3)/(x -> #x#0 * #x#1)//sum;
    contrib2 := (values interior2)/(x -> #x#0 * #x#1)//sum;
    # (latticePoints polar P) - 1 - 4 + contrib2 - contrib3
    )

batyrev = method()
batyrev Polyhedron := (P) -> (
    interior3 := findInteriors(0,P);
    interior2 := findInteriors(1,P);
    contrib3 := (values interior3)/(x -> #x#0 * #x#1)//sum;
    contrib2 := (values interior2)/(x -> #x#0 * #x#1)//sum;
    (# (latticePoints polar P) - 1 - 4, contrib2, - contrib3)
    )

end

restart    
installPackage "PolyhedralObjects"
installPackage "PolymakeInterface"
viewHelp PolymakeInterface

restart
load "mes-polymake.m2"
str = ///    1    0    0    0   -1    1    0    0    0    1   -1    1   -2
    0    1    0    0    1   -1    0    1    1   -1    1   -2    0
    0    0    1    0    1   -1    0    1    0    0   -1   -2    2
    0    0    0    1   -1    1   -1   -1   -1    1    1    0   -1///
M = matrixFromString str
M1 = matrix{for i from 0 to numColumns M - 1 list 1} || M
P = new Polyhedron from {"Points" => transpose M1}
runPolymake(P, "AmbientDim")
runPolymake(P, "Vertices")
runPolymake(P, "InteriorLatticePoints") -- bug if none
runPolymake(P, "BoundaryLatticePoints")

P1 = convexHull M
vertices P1
latticePoints P1
interiorLatticePoints P1
P2 = polar P1
vertices P2
latticePoints P2
interiorLatticePoints P2

findInteriors(0,P2)
h11 P2
batyrev P2
h11 P1
