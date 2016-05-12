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
    "matrixFromString",
    "findInteriors",
    "findPolarDual",
    "faceIndices",
    "polarDualFace",
    "latticePointFaces",
    "hodgeOfCYToricDivisor",
    "hodgeOfCYToricDivisors",
    "h11OfCY",
    "h21OfCY",
    "interiorLattice",
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
    l1 := transpose matrix{l};
    lp := latticePointFaces(P);
    n := #interiorLatticePoints(polar P,polarDualFace(P,lp#0#l1#0));
    if lp#0#l1#1 == 0 then {1,0,n} 
      else if lp#0#l1#1 == 1 then {1,n,0}
      else if lp#0#l1#1 == 2 then {1+n,0,0}
    )
hodgeOfCYToricDivisors = method()
hodgeOfCYToricDivisors(Polyhedron, HashTable) := (P, intP) -> (
    flatten for i from 0 to 3 list (
        flatten for f in keys intP#i list (
            pts := intP#i#f;
            for g in pts#0 list g => (
                n := #pts#1;
                if i == 0 then {1,0,n} 
                else if i == 1 then {1,n,0}
                else if i == 2 then {1+n,0,0}
                )
            )
        )
    )

h11OfCY = method()
h11OfCY(Polyhedron) := (P) -> (
    np := #(latticePoints polar P);
    l := for f in faces(1,polar P) list #(interiorLatticePoints f);
    t := sum l;
    l1 := for f in faces(2, polar P) list (#(interiorLatticePoints f))*(#(interiorLatticePoints(P,polarDualFace(polar P,faceIndices(polar P,f)))));
    t1 := sum l1;
    np - 5 - t + t1
    );

h11OfCY(Polyhedron, HashTable) := (P, interiors) -> (
    -- interiors is the result of 'interiorLattice P'
    np1 := 1 + (values interiors)/values//flatten/last/length//sum;
    np := #(latticePoints polar P);
    if np != np1 then error "oops";
    t := (values interiors#0)/last/length//sum;
    t1 := (values interiors#1)/(v -> #v#0 * #v#1)//sum;
    np - 5 - t + t1
    );

h21OfCY = method()
h21OfCY(Polyhedron) := (P) -> (
    np := #(latticePoints P);
    l := for f in faces(1,P) list #(interiorLatticePoints f);
    t := sum l;
    l1 := for f in faces(2,P) list (#(interiorLatticePoints f))*(#(interiorLatticePoints(polar P,polarDualFace(P,faceIndices(P,f)))));
    t1 := sum l1;
    np - 5 - t + t1
    );
h21OfCY(Polyhedron, HashTable) := (P, interiors) -> (
    -- interiors is the result of 'interiorLattice P'    
    np1 := 1 + (values interiors)/values//flatten/first/length//sum;
    np := #(latticePoints P);
    if np != np1 then error "oops";
    t := (values interiors#3)/first/length//sum;
    t1 := (values interiors#2)/(v -> #v#0 * #v#1)//sum;
    np - 5 - t + t1
    );
isFavorable(Polyhedron, HashTable) := (P, interiors) -> (
    t1 := (values interiors#1)/(v -> #v#0 * #v#1)//sum;
    t1 == 0
    )

interiorLattice = method()
interiorLattice Polyhedron := (P1) -> (
    time P2 := polar P1;
    time LP1 := select(latticePoints P1, v -> v != 0);
    time LP2 := select(latticePoints P2, v -> v != 0);
    facedims := time hashTable flatten for i from 0 to dim P1-1 list for f in faces(dim P1-i,P1) list faceIndices(P1,f) => i;
    time lp1 := for p in LP1 list (
        pf := positions(flatten entries ((transpose vertices P2) * p), x -> x == -1);
        f := polarDualFace(P2,pf);
        (f,p)
        );
    interiors1 := applyPairs(partition(f->f#0, lp1), (k,v) -> (k,v/last));
    time lp2 := for p in LP2 list (
        pf := positions(flatten entries ((transpose vertices P1) * p), x -> x == -1);
        --f := polarDualFace(P1,pf);
        (pf,p)
        );
    interiors2 := applyPairs(partition(f->f#0, lp2), (k,v) -> (k,v/last));
    actualfaces := sort toList (set keys interiors1 + set keys interiors2);
    result := for k in actualfaces list k => (if interiors1#?k then interiors1#k else {}, if interiors2#?k then interiors2#k else {});
    result = partition(f->facedims#(first f), result);
    applyPairs(result, (k,v) -> (k,hashTable v))
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
