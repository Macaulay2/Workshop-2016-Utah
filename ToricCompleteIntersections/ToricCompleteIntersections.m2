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

export {"vertexIndices","dualFaces","matrixFromString","findInteriors","findPolarDual","faceIndices","polarDualFace","latticePointFaces","hodgeOfCYToricDivisor"}

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

{*hodgeOfCYToricDivisor = method();
hodgeOfCYToricDivisor(Polyhedron,List) := (p,n) -> (
    n := transpose matrix {n}
    f := fa
    vlist := for p in f list vertices p;
    v := vertexIndices(P);
    n1 := v#n;
    pts := findPolarDual({n1},p,polar p) - 
    if member(n,vlist) then pts := interiorLatticePoints polar findPolarDual({n1},p,polar p);
    	hodgenums := {1,0,#pts};
    f1 := faces(3,p);
    f1list := for p1 in f1 list interiorLatticePoints p1;
    for i in f1list do if member(n,i) then pts := interiorLatticePoints polar findPolarDual({n1},p,polar p);
    	hodgenums := {1,#pts,0};
    i := position(f1list, 
    f2 := faces(2,p);
    f2list := for p2 in f2 list interiorLatticePoints p2;
    for i in p2 if member(n,i) then pts := interiorLatticePoints polar findPolarDual({n1},p,polar p);
    	hodgenums := {1 + #pts, 0,0};
    hodgenums
    );*}

hodgeOfCYToricDivisor = method();
hodgeOfCYToricDivisor(Polyhedron,List) := (P,l) -> (
    l1 := transpose matrix{l};
    lp := latticePointFaces(P);
    n := #interiorLatticePoints(polar P,polarDualFace(P,lp#0#l1#0));
    if lp#0#l1#1 == 0 then {1,0,n} 
      else if lp#0#l1#1 == 1 then {1,n,0}
      else if lp#0#l1#1 == 2 then {1+n,0,0}
    )


{*dualFaces(ZZ,Polyhedron) := (n,P) -> (
    Pv := polar P;
    
    );*}
    
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
str = ///    1    0    0    0   -1    1    0    0    0    1   -1    1   -2
    0    1    0    0    1   -1    0    1    1   -1    1   -2    0
    0    0    1    0    1   -1    0    1    0    0   -1   -2    2
    0    0    0    1   -1    1   -1   -1   -1    1    1    0   -1///
M = matrixFromString str
P = convexHull M
v = vertexIndices(P)
l = latticePoints P
l = select(l,p->p!= 0)
lp = for p in l list flatten entries lift(p,ZZ)
hodgeOfCYToricDivisor(P,{0,0,0,1})
for p in lp list p => hodgeOfCYToricDivisor(P,p)
hashTable oo
