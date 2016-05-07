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

randomAddition = (n,P) ->(
    --P must be non-empty
    d := #P_0 - 1;
    t := false;
    D' := {null};
    r := random (#P);
    D := P_r;
    compD := toList(0..n-1) - set D;
    count := 0;
    while not t and count < 20 do (
    	i := random (#compD);
    	j := random (#D);
    	D' = sort(D - set {D_j} | {compD_i});
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

manyRandomAdditions = (N,n,P) -> (
    unique apply(N, i-> P = randomAddition(n,P)))

///
Q = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}, {0, 3, 4}}
P = {1,3,4}
(t,smalls, facets) = testNewSimplex(Q,P)
smalls
facets
apply(smalls, e ->apply(facets, E -> #(e-set E)))
///

end
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


















