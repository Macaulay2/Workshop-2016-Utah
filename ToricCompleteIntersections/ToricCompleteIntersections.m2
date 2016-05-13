newPackage(
        "ToricCompleteIntersections",
        Version => "0.1", 
        Date => "",
        Authors => {{Name => "", 
                  Email => "", 
                  HomePage => ""}},
        Headline => "Hodge Diamond of Complete Intersections in Toric Varieties",
        DebuggingMode => true,
	Reload => true,
	PackageExports => {
	    --"FourTiTwo",
            --"Normaliz",
            --"gfanInterface",
            --"SimplicialComplexes",
            "NormalToricVarieties",
            --"Schubert2",
	    "Polyhedra"
	    --"PolyhedralObjects",
	    --"PolymakeInterface"
        })

export {"vertexIndices",
    "dualFaces",
    "findInteriors",
    "findPolarDual",
    "faceIndices",
    "polarDualFace",
    "latticePointFaces",
    -- The following are the ones we care about
    "matrixFromString",
    "hodgeOfCYToricDivisor",
    "hodgeOfCYToricDivisors",
    "h11OfCY",
    "h21OfCY",
    "interiorLattice", -- this should be internal?  Right now, we need to see it.
    "isFavorable"
    }

-- Code here
matrixFromString = method()
matrixFromString String := (str) -> (
    -- expect input in one of the three forms:
    {*
      (a)
              [ 1 -1 -1  1 -1 -1 -1  1  1  0]
              [ 0  1  0 -1  1  0  0 -1 -1  0]
              [-1  0  0  0  1  1  0  0  0  0]
              [ 2  0  0  1 -1 -1 -1 -1  0  0]
      (b)
           [[1, 0, -1, 2], [-1, 1, 0, 0], [-1, 0, 0, 0], [1, -1, 0, 1], [-1, 1, 1, -1], [-1, 0, 1, -1], [-1, 0, 0, -1], [1, -1, 0, -1]]    
      (c) 
        "         1   0   0   0   1   1   0  -1  -1  -2  -4
                  0   1   1   0  -2   2   3  -1  -4   1  -1
                  0   0   2   0  -2   4   4  -1  -4  -2  -4
                  0   0   0   1   0  -2  -2   2   2   0   2"
    *}
    s1 := replace(", ", " ", str);
    s2 := replace("] ", "]\n", s1);
    s3 := replace("\\[", "", s2);
    s4 := replace("]", "", s3);
    rows := for s in lines s4 list (
            t := separateRegexp(" +", s); 
            if first t == "" then t = drop(t,1);
            if last t == "" then t = drop(t,-1);
            if #t == 0 then continue;
            t/value
            );
    matrix rows
    )

isFavorable = method();

findPolarDual = (f, P, P2) -> (
    -- given a face 'f' of P2
    -- this returns the list of indices (with respect to 'vertices P') which give
    -- the face dual to f.
    positions(entries (transpose (vertices P) * (vertices P2)_f), x -> all(x, x1 -> x1 == -1))
    )

{*faces Polyhedron := (P) -> (
    v := entries transpose vertices P;
    -- careful! v#i is currently a vector over QQ.
    -- indexing using a list of ZZ elements won't give the same element.
    V := hashTable for i from 0 to #v-1 list v#i => i;
    hashTable for i from 1 to dim P list 
      (dim P)-i => for f in faces(i,P) list (
        for p in entries transpose vertices f list V#p
        )
    )*}

vertexIndices = (P) -> (
    if not P.cache.?vertexIndices then P.cache.vertexIndices = (
    	v := entries transpose lift(vertices P,ZZ);
    	hashTable for i from 0 to #v-1 list v#i => i
	);
    P.cache.vertexIndices
    );

{*findInteriors = (facedim,P) -> (
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
   )*}
findInteriors = method();
findInteriors(ZZ,Polyhedron) := (facedim,P) -> (
    hashTable for f in faces(dim P - facedim,P) list (
	fi := faceIndices(P,f);
	gi := polarDualFace(P,fi);
	(fi,gi) => (interiorLatticePoints(P,fi),interiorLatticePoints(polar P,gi))
	)
    );

findInteriors(Polyhedron) := (P) -> (
    hashTable for i from 0 to dim P - 1 list (
	i => findInteriors(i,P)
	)
    );

latticePointFaces = method();
latticePointFaces(Polyhedron) := (P) -> (
    H := findInteriors(P);
    result := new MutableHashTable;
    resultpolar := new MutableHashTable;
    for i in keys(H) do (
	L := H#i;
	for j in keys(L) do (
	    for lp in L#j#0 do result#lp = (j#0,i);
	    for lp in L#j#1 do resultpolar#lp = (j#1,dim P - i - 1);
	    )
	);
    (result,resultpolar)
    );

{*pointInterior = (p,n) -> (
    v := entries transpose vertices P
    if member(n,v) then return 0;
    f1 := *}
    
interiorLatticePoints(Polyhedron,List) := (P,f) -> (
    if not P.cache.?interiorLatticePoints then P.cache.interiorLatticePoints = new MutableHashTable;
    if not P.cache.interiorLatticePoints#?f then (
	P.cache.interiorLatticePoints#f = interiorLatticePoints convexHull (vertices P)_f;
	);
    return P.cache.interiorLatticePoints#f
    );

faceIndices = method();
faceIndices(Polyhedron,Polyhedron) := (P,f) -> (
    v := entries transpose lift(vertices f,ZZ);
    V := vertexIndices(P);
    for v1 in v list V#v1
    );

polarDualFace = method();
polarDualFace(Polyhedron,List) := (P,f) -> (
    vp := vertices P;
    vpolar := vertices polar P;
    positions(entries (transpose (vpolar) * (vp)_f), x -> all(x, x1 -> x1 == -1))
    );

hodgeOfCYToricDivisor = method();
hodgeOfCYToricDivisor(Polyhedron,List) := (P,l) -> (
    P2 := polar P;
    l1 := transpose matrix{l};
    lp := latticePointFaces(P2);
    n := #interiorLatticePoints(P,polarDualFace(P2,lp#0#l1#0));
    if lp#0#l1#1 == 0 then {1,0,n} 
      else if lp#0#l1#1 == 1 then {1,n,0}
      else if lp#0#l1#1 == 2 then {1+n,0,0}
    )

hodgeOfCYToricDivisors = method()
hodgeOfCYToricDivisors Polyhedron := (P) -> (
    if dim P != 4 then error "expected a 4d reflexive polytope";
    intP := interiorLattice P;
    flatten for i from 0 to 3 list (
        flatten for f in keys intP#(3-i) list (
            pts := intP#(3-i)#f;
            for g in pts#1 list g => (
                n := #pts#0;
                if i == 0 then {1,0,n} 
                else if i == 1 then {1,n,0}
                else if i == 2 then {1+n,0,0}
                )
            )
        )
    )

h11OfCY = method()
h21OfCY = method()

-- original code, changed by Mike to avoid use of faces(...,P)
h11OfCY(Polyhedron,Thing) := (P,notused) -> (
    np := #(latticePoints polar P);
    l := for f in faces(1,polar P) list #(interiorLatticePoints f);
    t := sum l;
    l1 := for f in faces(2, polar P) list (#(interiorLatticePoints f))*(#(interiorLatticePoints(P,polarDualFace(polar P,faceIndices(polar P,f)))));
    t1 := sum l1;
    np - 5 - t + t1
    );

-- original code, changed by Mike to avoid use of faces(...,P)
h21OfCY(Polyhedron,Thing) := (P,notused) -> (
    np := #(latticePoints P);
    l := for f in faces(1,P) list #(interiorLatticePoints f);
    t := sum l;
    l1 := for f in faces(2,P) list (#(interiorLatticePoints f))*(#(interiorLatticePoints(polar P,polarDualFace(P,faceIndices(P,f)))));
    t1 := sum l1;
    np - 5 - t + t1
    );

h11OfCY Polyhedron := (P) -> (
    interiors := interiorLattice P;
    np1 := 1 + (values interiors)/values//flatten/last/length//sum;
    np := #(latticePoints polar P);
    if np != np1 then error "oops";
    t := (values interiors#0)/last/length//sum;
    t1 := (values interiors#1)/(v -> #v#0 * #v#1)//sum;
    np - 5 - t + t1
    );

h21OfCY Polyhedron := (P) -> (
    interiors := interiorLattice P;
    np1 := 1 + (values interiors)/values//flatten/first/length//sum;
    np := #(latticePoints P);
    if np != np1 then error "oops";
    t := (values interiors#3)/first/length//sum;
    t1 := (values interiors#2)/(v -> #v#0 * #v#1)//sum;
    np - 5 - t + t1
    );

isFavorable Polyhedron := (P) -> (
    -- This is from the Batyrev formula for h^11, 
    -- The term ell^*(theta) * ell^*(theta^*) gives new divisors.
    --  here theta is an edge of P, theta^* is a (dim P)-2 face of (polar P)
    --  and ell^* is the number of interior lattice points.
    interiors := interiorLattice P;
    t1 := (values interiors#1)/(v -> #v#0 * #v#1)//sum;
    t1 == 0
    )

-- TODO: this function has an awful name
-- This function takes a polytope, and returns a hash table with keys the dimension i of faces
-- whose key is also a hashtable with keys being the faces of P of dimension i (given by 
--  an ascending list of integer indices of the vertices of the face), 
-- and the value at that index is a pair:
--   (list of interior points to that face of P, list of interior points to the dual face in polar P)
-- If both of these are empty lists, then that face is not placed in as a key.
-- It is stashed inside the Polyhedron object.
-- This function uses the following potentially low functions:
--  latticePoints
--  faces.  We could potentially avoid this by recomputing the dimension of
interiorLattice = method()
interiorLattice Polyhedron := (P1) -> (
    if not P1.cache.?interiorLattice then P1.cache.interiorLattice = (
        P2 := polar P1;
        LP1 := select(latticePoints P1, v -> v != 0);
        LP2 := select(latticePoints P2, v -> v != 0);
        facedims := hashTable flatten for i from 0 to dim P1-1 list for f in faces(dim P1-i,P1) list faceIndices(P1,f) => i;
        lp1 := for p in LP1 list (
            pf := positions(flatten entries ((transpose vertices P2) * p), x -> x == -1);
            f := polarDualFace(P2,pf);
            (f,p)
            );
        interiors1 := applyPairs(partition(f->f#0, lp1), (k,v) -> (k,v/last));
        lp2 := for p in LP2 list (
            pf := positions(flatten entries ((transpose vertices P1) * p), x -> x == -1);
            --f := polarDualFace(P1,pf);
            (pf,p)
            );
        interiors2 := applyPairs(partition(f->f#0, lp2), (k,v) -> (k,v/last));
        actualfaces := sort toList (set keys interiors1 + set keys interiors2);
        result := for k in actualfaces list k => (if interiors1#?k then interiors1#k else {}, if interiors2#?k then interiors2#k else {});
        result = partition(f->facedims#(first f), result);
        applyPairs(result, (k,v) -> (k,hashTable v))
        );
    P1.cache.interiorLattice
    )
    
beginDocumentation()

doc ///
Key
  ToricCompleteIntersections
Headline
  Hodge Diamond of Complete Intersections in Toric Varieties
Description
  Text
  Example
Caveat
SeeAlso
///

TEST ///
  -- We work on one example in 4 dimensions, where we know the answers (or have computed them elsewhere).
  -- Second polytope (index 1) on h11=3 Kreuzer-Skarke list of 4d reflexive polytopes for h11=3.
  -- 4 5  M:48 5 N:8 5 H:3,45 [-84]
  str = "  1   0   2   4  -8
         0   1   5   3  -9
         0   0   6   0  -6
         0   0   0   6  -6
  "
  M = matrixFromString str
  assert(M == matrix {{1, 0, 2, 4, -8}, {0, 1, 5, 3, -9}, {0, 0, 6, 0, -6}, {0, 0, 0, 6, -6}})

  P = convexHull M  
  assert(# latticePoints P == 48)
  assert(# latticePoints polar P == 8)
  assert(h11OfCY P == 3)
  assert(h21OfCY P == 45)
  assert(h11OfCY polar P == 45)
  assert(h21OfCY polar P == 3)
  
  -- Now compute all of the cohomologies of the (irreducible) toric divisors 
  LP = select(latticePoints polar P, p -> p != 0)
  LP = for lp in LP list flatten entries lift(lp,ZZ)
  assert(#LP == 7)
  cohoms = for v in LP list hodgeOfCYToricDivisor(P, v)
  cohomH = hashTable hodgeOfCYToricDivisors P
  cohoms1 = for v in LP list cohomH#(transpose matrix {v})
  assert(cohoms == cohoms1)
  assert isFavorable P

  -- Now compute all of the cohomologies of the (irreducible) toric divisors for the polar dual
  LP = select(latticePoints P, p -> p != 0)
  LP = for lp in LP list flatten entries lift(lp,ZZ)
  assert(#LP == 47)
  cohoms = for v in LP list hodgeOfCYToricDivisor(polar P, v)
  cohomH = hashTable hodgeOfCYToricDivisors polar P
  cohoms1 = for v in LP list cohomH#(transpose matrix {v})
  assert(cohoms == cohoms1)
  assert not isFavorable polar P
///

end--

doc ///
Key
Headline
Usage
Inputs
Outputs
Consequences
Description
  Text
  Example
  Code
  Pre
Caveat
SeeAlso
///

TEST ///
-- test code and assertions here
-- may have as many TEST sections as needed
{*vertexIndices = (P) -> (
    v := entries transpose vertices P;
    hashTable for i from 0 to #v-1 list v#i => i
    );*}
restart
loadPackage "ToricCompleteIntersections"
str = ///    1    0    0    0  -11   -3   -1   -3
    0    1    0    0   -6   -2   -2   -6
    0    0    1    0   -2   -2   -2   -2
    0    0    0    1   -2    2    4    6 ///
M = matrixFromString str
P = convexHull M
v = vertexIndices(P)
l = latticePoints P
l = select(l,p->p!= 0)
lp = for p in l list flatten entries lift(p,ZZ)
hodgeOfCYToricDivisor(P,{0,0,0,1})
for p in lp list p => hodgeOfCYToricDivisor(P,p)
hashTable oo
latticePointFaces(P)
new HashTable from first oo
findInteriors(P)
h11OfCY(P)
h21OfCY(P)
