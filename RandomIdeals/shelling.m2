--goal: construct shellable complexes at random

newPackage ( "shelling",
    Version => "1.0",
    Date => "07 May 2016",
    Authors => {
        {Name => "Katy"},
	{Name => "David Eisenbud",
         Email => "de@msri.org",
         HomePage => "http://www.msri.org/~de"},
	{Name => "Robert"},
	{Name => "Jay"}
	},
    Headline => "Package for constructing random simplicial complex",
    Reload => true,
    DebuggingMode => true
    )


export {
	"randomAddition", 
	"randomChain",
	"randomLink",
	"testNewSimplex",
        "idealFromSC",
	"idealChainFromSC",
        "isShelling",
	"isLicci",
	"minimalRegularSequence",
	"linkageBound",
	"UseNormalModule",
	"randomRegularSequence"
        };

testNewSimplex = method()
testNewSimplex(List, List) := (P, D) ->(
--given a pure, d-dimensional simplicial complex (sc) as a list of ordered lists of d+1 vertices in [n], and
--a simplex D as such a list, tests whether the intersection of D with P is a union of facets of D.
     d := #D-1; --dimension
     ints := unique apply(P, D' -> intersectLists(D',D));
     addSimplex := (L,S) -> (
         if any(S,F -> subsetList(L,F))
         then S
         else select(S,F -> not (subsetList(F,L))) | {L}
         );
     ints = fold(ints,{}, addSimplex);
     all(ints, L -> #L == d)
)

subsetList = (A,B) -> (
    --checks if A\subset B. requires that both lists be sorted
    lenA := #A;
    lenB := #B;
    a := 0;
    b := 0;
    while (lenA-a)<=(lenB-b) and a<lenA and b<lenB do (
        if A_a<B_b then return false;
        if A_a==B_b then a=a+1;
        b = b+1;
    );
    a==lenA
)

intersectLists = (D',D) -> D - set(D-set D')

randomSubset = (n,m) -> (
    L := new MutableList from toList (0..m-1);
    for i from m to n-1 do (
	j := random(i+1);
    	if j < m then L#j = i;
	);
    sort toList L
    )

randomAddition = method()
randomAddition(ZZ,ZZ,List) := (n,m,P) ->(
    if #P == 0 then return {randomSubset(n,m+1)};
    Plarge := select(P, D-> #D >= m); -- the facets big enough to be glued to
    if #Plarge == 0 then error "m is too large";
    t := false;
    D' := {null};
    D := Plarge#(random(#Plarge)); -- a random facet from Plarge
    compD := toList(0..n-1) - set D;
    count := 0;
    while not t and count < 20 do (
    	i := random (#compD);
    	J := randomSubset(#D,#D-m);
    	D' = sort(D - set apply(J, j->D_j) | {compD_i});
    	t = (testNewSimplex(P,D'));
	count = count+1);
    if count == 20 then return P;
    unique (P|{D'})
    )

randomAddition(Ring,ZZ,List) := (R,m,L) -> (
    P := monomialsToLists(L,R);
    listsToMonomials(randomAddition(numgens R,m,P),R)
    )

idealFromSC = method()
idealFromSC (List,Ring) := (P,S) -> (
    Delta := toList (0..numgens S - 1);
    V := vars S;
    intersect apply(P, D -> ideal(V_(Delta - set D)))
    )

idealFromSC List := P -> (
    n := (max flatten P)+1;
    x := symbol x;
    S := QQ[x_0..x_(n-1)];
    idealFromSC(P,S)
    )

idealChainFromSC = method()
idealChainFromSC List := P -> toList apply(#P,i->idealFromSC(take(P,i+1)))

isShelling = method()
isShelling(List) := P -> all(#P, i-> i==0 or testNewSimplex(take(P,i),P#i))

randomChain = method()
-- random chain of shellable complexes on n vertices, with pure dim m, up to the complete m skeleton

randomChain(ZZ,ZZ) := (n,m) -> randomChain(n,m,binomial(n,m+1))
-- random chain of shellable complexes on n vertices, with pure dim m, and k facets

--Should we change the following to start with {{0..m}, {0..m-1,m} to diminish autos?
randomChain(ZZ,ZZ,ZZ) := (n,m,k) -> (
    if k > binomial(n,m+1) then error "k is too large";
    P := {};
    while #P < k do P = randomAddition(n,m,P);
    P
    )

randomChain(Ring,ZZ,ZZ) := (R,m,k) -> listsToMonomials(randomChain(numgens R,m,k),R)
randomChain(Ring,ZZ)    := (R,m)   -> listsToMonomials(randomChain(numgens R,m),R)

--this is NOT the Reisner association
listsToMonomials = (P,R) -> apply(P, D->product apply(D,d->R_d))
monomialsToLists = (L,R) -> apply(L, m->select(numgens ring m, i->((listForm m)#0#0#i > 0)))


------------------------------------------------------------
-- DOCUMENTATION randomChain
------------------------------------------------------------

doc ///
     Key
          randomChain
	  (randomChain,ZZ,ZZ)
	  (randomChain,ZZ,ZZ,ZZ)
	  (randomChain,Ring,ZZ)
	  (randomChain,Ring,ZZ,ZZ)
     Headline
          produces a random chain of shellable complexes
     Usage
          P=randomChain(n,m)
	  P=randomChain(n,m,k)
	  P=randomChain(R,m)
	  P=randomChain(R,m,k)
     Inputs
          n:ZZ
	       the number of vertices
	  R:Ring
	       a polynomial ring with a variable for each vertex
	  m:ZZ
	       the dimension of the facets
	  k:ZZ
	       the number of facets (if ommited, the number will be {\tt n} choose {\tt m+1})
	      
     Outputs
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.  If called with a Ring {\tt R} instead of an integer {\tt n}, each facet is represented by a square-free monomial instead of a list.
     Description
          Text
               
          Example
               P = randomChain(6,3,10)
     Caveat
	  No claim is made on the distribution of the random chain.
///





------------------------------------------------------------
-- DOCUMENTATION isShelling
------------------------------------------------------------
doc ///
     Key
          isShelling
	  (isShelling,List)
     Headline
          determines whether a list represents a shelling of a simplicial complex.
     Usage
          b = isShelling(P)
     Inputs
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a possible shelling.
     Outputs
          B:Boolean
	       true if and only if P is a shelling.
     Description
          Text
              Determines if a list of faces is a shelling order of the simplicial complex. 
          Example
	      P = {{1, 2, 3}, {1, 2, 5}};
	      isShelling(P)
	      Q = {{1,2,3},{3,4,5},{2,3,4}};
	      isShelling(Q)
	     
      
///


------------------------------------------------------------
-- DOCUMENTATION randomAddition
------------------------------------------------------------
doc ///
     Key
          randomAddition
	  (randomAddition,ZZ,ZZ,List)
     Headline
          Adds a random facet to a shellable complex
     Usage
          p=randomAddition(n,m,P)
     Inputs
     	  n:ZZ
	       the number of vertices
	  m:ZZ
	       the dimension of the new facet
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.
     Outputs
          p:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.
     Description
          Text
            This function randomly chooses a facet of size m+1 and checks whether the facet can be shellably added to the shelling. If it can be shellably added to the shelling, it is added to the shelling and the new shelling is returned. Otherwise, the process repeats up to 20 times.  
          Example
            P={{1,2,3}}
	    L=randomAddition(6,3,P)
     Caveat
	  If the input is not a shellable simplicial complex, the new complex will not be shellable.
///

------------------------------------------------------------
-- DOCUMENTATION testNewSimplex
------------------------------------------------------------

doc ///
     Key
          testNewSimplex
	  (testNewSimplex,List,List)
     Headline
          Tests whether a facet can be shellably added to a shelling.
     Usage
          b=testNewSimplex(P,S)
     Inputs
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.
          S:List
               A list of integers. This list is the new facet to add.
     Outputs
          b:Boolean
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.
     Description
          Text
               
          Example
            P={{1,2,3}}
	    b=testNewSimplex(P,{2,3,4});
     Caveat
          We do not test if P is a shelling in the first place.
///
         

TEST///
assert(#randomChain(5,2,6)==6)
assert(#randomChain(5,2)==binomial(5,3))
R=QQ[x1,x2,x3,x4,x5];
assert(#randomChain(R,2,6)==6)
///


TEST///
assert(isShelling({}))
assert(isShelling({{1,2,3}}))
assert(isShelling({{1,2,3},{2,3,4}}))
assert(isShelling(randomChain(5,3,5)))
--non pure shellings
assert(isShelling({{1,2,3},{2,4}}))
assert(isShelling({{1},{2}}))
assert(not isShelling({{1,3},{2,4}}))
assert(isShelling({{1,2},{3}}))
assert(not isShelling({{3},{1,2}}))
///


TEST///
setRandomSeed(0);
assert(#randomAddition(6,3,{{1,2,3}})==2)
assert(#randomAddition(6,3,{{1,2,3,4}})==2)
///

TEST///
needsPackage "SimplicialComplexes"
needsPackage "SimplicialDecomposability"
R=QQ[x1,x2,x3,x4,x5];
assert(isShellable simplicialComplex randomChain(R,2,6))
///

end--
restart
uninstallPackage "RandomIdeals"
installPackage "RandomIdeal"

restart
installPackage "shelling"

restart
loadPackage("shelling", Reload=>true)
check "shelling"

R = ZZ/32003[x_0..x_4]
P = randomChain(5,1)
L = apply(#P, i->idealFromSC take(P,i+1));
netList L
apply(L, I->linkageBound(I,UseNormalModule=>true))
