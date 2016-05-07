--goal: construct shellable complexes at random

testNewSimplex = method()
testNewSimplex(List, List) := (P, D) ->(
--given a pure, d-dimensional simplicial complex (sc) as a list of ordered lists of d+1 vertices in [n], and
--a simplex D as such a list, tests whether the intersection of D with P is a union of facets of D.
     d := #D-1; --dimension
     ints := apply(P, D' -> intersectLists(D',D));
     facets := unique select(ints, E -> #E==d);
     if facets == {} then return false;
     smalls := unique select(ints, E -> #E<d);
     if sum apply(smalls, e ->product apply(facets, E ->  #(e-set E)))===0 then t=true else t=false;
--error();
(t,smalls,facets)
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

randomAddition = (n,m,P) ->(
    if #P == 0 then return {randomSubset(n,m)};
    Plarge := select(P, D-> #D >= m-1); -- the facets big enough to be glued to
    if #Plarge == 0 then error "m is too large";
    t := false;
    D' := {null};
    D := Plarge#(random(#Plarge)); -- a random facet from Plarge
    compD := toList(0..n-1) - set D;
    count := 0;
    while not t and count < 20 do (
    	i := random (#compD);
    	J := randomSubset(#D,#D-m+1);
    	D' = sort(D - set apply(J,j->D_j) | {compD_i});
    	t = (testNewSimplex(P,D'))_0;
	count = count+1);
    if count == 20 then return P;
    unique (P|{D'})
    )

idealFromSC = (P) ->(
    numverts := #unique flatten P;
    S := ZZ/101[x_0..x_(numverts-1)];
    Delta := toList (0..numgens S -1);
    V = vars S;
    intersect apply(P, D -> ideal(V_(Delta - set D)))
	    )

        

randomChain = method()
-- random chain of shellable complexes on n vertices, with pure dim m, up to the complete m skeleton
randomChain(ZZ,ZZ) := (n,m) -> randomChain(n,m,binomial(n,m))
-- random chain of shellable complexes on n vertices, with pure dim m, and k facets
randomChain(ZZ,ZZ,ZZ) := (n,m,k) -> (
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
load "shelling.m2"
S = ZZ/101[x_0..x_3]
I = ideal(x_0*x_1,x_1^2, x_2^3, x_3^5)
isLicci(3, codim I, I)

I = minors(3, random(S^3, S^{-2,-3,-4,-4}));
isLicci(3, codim I, I)

I = minors(2, random(S^2, S^{4:-1}))
isLicci(3, codim I, I)

c = codim I
b = linkageBound I
b = 2*c
count = 0

while(numgens I > c and count < b) list (
    I = randomLink(c, I);
    count = count+1;
    <<numgens I<<endl;
    flush;)

///

end--
viewHelp

restart
uninstallPackage "RandomIdeal"
installPackage "RandomIdeal"

restart
load "shelling.m2"
P = {{1,2,3}}
netList (L =  manyRandomAdditions(500,6,P))
netList (L1 =  manyRandomAdditions(500,6,P))
#L
apply(19, i -> (
    <<betti res idealFromSC(L_i)<< endl;
    <<betti res idealFromSC(L1_i)<< endl;
<<endl;<<endl;))
#unique(L|L1)
netList apply(L, P->betti res idealFromSC(P))

degree I
betti res I
betti res (I^3)
Q = subsets(5,3)
betti res idealFromSC(Q)



--an apparently maximal shellable with non-linear resolution:
P = {{1,2,3},
 {1,2,5},
 {2,3,7},
 {1,5,6},
 {0,2,3},
 {0,2,7},
 {1,4,5},
 {4,5,6}
 }
I = idealFromSC(P)
betti res I

netList (L =  manyRandomAdditions(500,8,P))
betti res idealFromSC(last L)

M = apply(1000, i->(
	L = {};
	P = {{1,2,3}};
	manyRandomAdditions(500,6,P)));
#unique flatten M

P = {{1,2,3}}
# (L=manyRandomAdditions(500,6,P))
L


















