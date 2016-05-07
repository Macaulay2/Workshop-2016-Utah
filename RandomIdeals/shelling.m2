--goal: construct shellable complexes at random

newPackage ( "shelling",
    Version => "1.0",
    Date => "07 May 2016",
    Authors => {
	{Name => "David Eisenbud",
         Email => "de@msri.org",
         HomePage => "http://www.msri.org/~de"},
     	 {Name => "Robert,Katy,Robert, Jay"}
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
        "isShelling",
	"isLicci",
	"minimalRegularSequence"
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

isShelling = method()
isShelling(List) := P -> all(#P, i-> i==0 or testNewSimplex(take(P,i),P#i))

randomChain = method()
-- random chain of shellable complexes on n vertices, with pure dim m, up to the complete m skeleton
randomChain(ZZ,ZZ) := (n,m) -> randomChain(n,m,binomial(n,m+1))
-- random chain of shellable complexes on n vertices, with pure dim m, and k facets
randomChain(ZZ,ZZ,ZZ) := (n,m,k) -> (
    if k > binomial(n,m+1) then error "k is too large";
    P := {};
    while #P < k do P = randomAddition(n,m,P);
    P
    )



///
Q = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}, {0, 3, 4}}
P = {1,3,4}
(t,smalls, facets) = testNewSimplex(Q,P)
smalls
facets
apply(smalls, e ->apply(facets, E -> #(e-set E)))
///


-----Toward the test for licci


randomLink = method()
randomLink (ZZ,Ideal) := (c,I) ->(
{*
c:ZZ
 codim of I
I:Ideal
 homogeneous
*}
if numgens I <= c then return ideal(1_(ring I));
--sgens := sort gens I;
--n :=numcols sgens;
--rsgens  := sgens * random(source sgens, source sgens);
--regseq := ideal rsgens_{n-c..n-1};
regseq := minimalRegularSequence(c,I);
trim(regseq : I)
)

linkageBound = method()
linkageBound Ideal := I ->(
    --2(mu N(I) - (codim I +3))
N := prune Hom(I, (ring I)^1/I);
n := numgens N;
2*(n-codim I)
)

minimalRegularSequence = method()
minimalRegularSequence(ZZ,Ideal) := (c,I) ->(
{*
c:ZZ
 codim of I
I:Ideal
 homogeneous
Description
 Text
  finds a maximal regular sequence in I of minimal degree.
*}
if numgens I == c then return I;
    --takes care of I = 0 and I principal;
sgens := sort gens I;
rsgens := sgens * random(source sgens, source sgens);
n :=numcols sgens;
J := ideal sgens_{0};
K := J;
count := 1; -- next element to add
c' := 1; -- current codim J
while c'<c do(
    if codim (K = J + ideal sgens_{count}) > c' then (J = K; c' = c'+1)
    else if codim (K = J + ideal rsgens_{count}) > c' then (J = K; c' = c'+1);
    count = count+1;
    );
J
)

isLicci = method()
isLicci(ZZ, ZZ, Ideal) := (b,c,I) -> (
    --I homogeneous ideal
    --b = linkageBound I
    --c = codim I
    --output is list of up to b integers, the numbers of generators of the
    --successive random links
    J := I;
    p := numgens J;
    <<p<<endl;flush;
    apply(b, i -> (
	    J = randomLink(c,J);
	    <<numgens J<<endl;flush;
	    numgens J))
    )
   
///
restart
loadPackage("shelling", Reload =>true)
S = ZZ/101[x_0..x_3]
I = ideal(x_0*x_1,x_1^2, x_2^3, x_3^5)
isLicci(3, codim I, I)

I = minors(3, random(S^3, S^{-2,-3,-4,-4}));
isLicci(3, codim I, I)

I = minors(2, random(S^2, S^{4:-1}))
isLicci(3, codim I, I)

--b = linkageBound I
b = 2*c
///

------------------------------------------------------------
-- DOCUMENTATION randomChain
------------------------------------------------------------

doc ///
     Key
          randomChain
	  (randomChain,ZZ,ZZ)
	  (randomChain,ZZ,ZZ,ZZ)
     Headline
          produces a random chain of shellable complexes
     Usage
          b=isShelling(P)
     Inputs
          P:List
	       the number of vertices
     Outputs
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a shelling.
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
          b = isShelling P
     Inputs
          P:List
	       A list of lists of integers.  Each list of integers is a facet of the complex and the order is a possible shelling.
     Outputs
          b:Boolean
	       Returns true if the given list gives a shelling of the complex and false otherwise. 
     Description
          Text
               
          Example
  
///

TEST///
assert(#randomChain(5,2,6)==6)
assert(#randomChain(5,2)==binomial(5,3))
///


TEST///
assert(isShelling(randomChain(5,3,5)))
///


end--


restart
installPackage "shelling"

restart
loadPackage("shelling", Reload=>true)
check "shelling"