newPackage ( "ResidualIntersections",
    Version => "1.0",
    Date => "07 May 2016",
    Authors => {
	{Name => "David Eisenbud",
         Email => "de@msri.org",
         HomePage => "http://www.msri.org/~de"},
     	 {Name => "Robert,Katy,Robert, Jay"}
	},
    Headline => "Package for studying conditions associated to Residual Intersection theory",
    Reload => true,
    DebuggingMode => true
    )

export {
	"isLicci",
	"minimalRegularSequence",
	"linkageBound",
	"UseNormalModule",
	"randomRegularSequence",
	"genericResidual",
	"genericArtinNagata",
	"numgensByCodim",
	"maxGd",
	"residualCodims",
        "koszulDepth",
        "isStronglyCM"
        };

--Generic Artin-Nagata Code
genericArtinNagata = method()
genericArtinNagata(ZZ,Ideal) := (s,I) -> (
    needsPackage "MCMApproximations";
    S := ring I;
    K := genericResidual(s,I);
    s' := codim K;
    if s' === s then 
      codepth := numgens (ring K) - depth ((ring K)^1/K)
    else codepth = -1;
    {s',codepth,K}
    )
    --tests whether the generic link is CM

genericResidual = method()
genericResidual(ZZ,Ideal):= (s,I) ->(
    if s>= numgens I then return ideal(1_(ring I));
    sgens := sort gens I;
    rgens := (sgens)*random(source sgens, source sgens);
    n := numcols rgens;
    (ideal (rgens_{n-s..n-1})): I
    )

///
restart
loadPackage("RandomIdeal", Reload => true)
loadPackage("ResidualIntersections", Reload =>true)
S = ZZ/32003[x_0..x_5]
--6 vars
I = randomShellableIdeal(S,2,4)
S = ZZ/32003[x_0..x_3]
I = minors(3, random(S^3, S^{-2,-3,-4,-5}));
--codim 3
codim I
s = 2;
codim genericResidual(3,I)
L = genericArtinNagata(s,I);
L_{0,1}
///
---Licci code
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
randomLink Ideal := I->randomLink(codim I, I)

linkageBound = method(Options => {UseNormalModule =>false})
linkageBound Ideal := opts -> I -> (
if opts.UseNormalModule == false then
    2*(codim I)*(degree I -1) -6
 else (
	N := prune Hom(I, (ring I)^1/I);
	2*(numgens N - codim I))
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
minimalRegularSequence Ideal := I -> minimalRegularSequence(codim I, I)

isLicci = method(Options => {UseNormalModule =>false})
isLicci(ZZ, ZZ, Ideal) := opts -> (b,c,I) -> (
    --I homogeneous ideal
    --b = linkageBound I
    --c = codim I
    --output is list of up to b integers, the numbers of generators of the
    --successive random links
    J := I;
    p := numgens J;
--    <<p<<endl;flush;
    apply(b, i -> (
	    J = randomLink(c,J);
	    if numgens J == c then break;
	    <<numgens J<<endl;flush;
	    numgens J))
    )
isLicci(ZZ,Ideal) := opts -> (b,I) -> isLicci(b,codim I, I)
isLicci Ideal := opts -> I -> (
isLicci(linkageBound(I, UseNormalModule => opts.UseNormalModule), I
    ))

depth := profondeur;


profondeur = method()
profondeur(Ideal, Module) := (I,M) ->(
    --requires R to be an affine ring (eg NOT ZZ[x])
    R := ring M;
    d := max(1,dim M); -- d=0 causes a crash
    if not isCommutative R then error"profondeur undefined for noncommutative rings";
    F := M**dual res (R^1/I, LengthLimit => d);
    i := 0;
    while HH_i F == 0 do i=i-1;
    -i)

profondeur Module := M -> (
    --profondeur of a module with respect to the max ideal, via finite proj dim
    --gives error if the ultimate coeficient ring of R = ring M is not a field.
    R := ring M;
    if not isCommutative R then error"profondeur undefined for noncommutative rings";
    (S,F) := flattenRing R;
    if not isField coefficientRing S then error"input must be a module over an affine ring";
    S0 := ring presentation S;
    r := F*map(S,S0);
    MM := pushForward(r,M);
    numgens S0 - pdim MM)

profondeur Ring := R -> profondeur R^1

koszulDepth = method()
koszulDepth(Ideal) := I -> (
    C := koszul mingens I;
    for i in 0..(numColumns(mingens I)-codim I) list depth HH_i(C)
    )

isStronglyCM = method()
isStronglyCM(Ideal) := I -> (
    d := dim I;
    all(koszulDepth I,i -> i==d)
    )

-------------------------------------
-- G_d Code
-------------------------------------

numgensByCodim = method()	
numgensByCodim (Ideal,ZZ) := (J,k) -> (
    R := ring J;
    n := numgens R;
    max for A in subsets(n,k) list (
	M := new MutableList from (n:1_R);
	for a in A do M#a = R_a;
	M = map(R,R,matrix{toList M});
	numgens trim M J
	)
    )

numgensByCodim Ideal := J -> (
    n := numgens ring J;
    toList apply(n, i->numgensByCodim(J,i+1))
    )

maxGd = method()
maxGd Ideal := J -> (
    for i from 1 to numgens ring J do if numgensByCodim(J,i) > i then return i;
    infinity
    )

residualCodims = method()
residualCodims Ideal := J -> (
    toList select((codim J + 1..numgens ring J + 1), i->numgensByCodim(J,i-1) <= i)
    )

------------------------------------------------------------
-- DOCUMENTATION isLicci
------------------------------------------------------------
doc ///
   Key
    isLicci
    (isLicci,ZZ,ZZ,Ideal)
    (isLicci,ZZ,Ideal)
    (isLicci,Ideal)
    [isLicci,UseNormalModule]
   Headline
    Tests whether an ideal is licci
   Usage
    L = isLicci(b,c,I)
    L = isLicci(b,I)
    L = isLicci I
   Inputs
    b:ZZ
     upper bound on how many successive random links to try
    c:ZZ
     codim of I
    I:Ideal
   Outputs
    L:List
     of integers, the numbers of generators of the successive links
   Description
    Text
     Computes up to b successive random links,
     using a regular sequence among the generators of I, and outputs the 
     numbers of generators. If I is licci, such a sequence must terminate
     in an ideal with c = codim I generators in at most
     linkageBound I steps.
    Example
   Caveat
    linkageBound I can be very large; linkageBound(I, UseNormalModule => true) can be very slow.
   SeeAlso
    linkageBound
    randomRegularSequence
    randomLink
///
------------------------------------------------------------
-- DOCUMENTATION UseNormalModule
------------------------------------------------------------
doc ///
   Key
    UseNormalModule
   Headline
    option for linkageBound and isLicci
   Description
    Text
     Default value is false. When true, it
     enables a more refined computation of the bound on the number of general links of
     an ideal I
     that must be taken to definitively test the licci propery. 
     When UseNormalModule = true the computation of the 
     normal module Hom(I, (ring I)/I) is required and this can be slow;
     if UseNormalModule = false the computation is fast, but the bound is large.
   SeeAlso
    isLicci
    linkageBound
///

------------------------------------------------------------
-- DOCUMENTATION linkageBound
------------------------------------------------------------

doc ///
   Key
    linkageBound
    (linkageBound,Ideal)    
    [linkageBound, UseNormalModule]
   Headline
    computes a bound on the number of general links of an ideal to test the licci property
   Usage
    b = linkageBound I
   Inputs
    I:Ideal
   Outputs
    b:ZZ
   Description
    Text
    Example
   Caveat
   SeeAlso
    isLicci
///

------------------------------------------------------------
-- DOCUMENTATION minimalRegularSequence
------------------------------------------------------------
doc ///
   Key
    minimalRegularSequence
    (minimalRegularSequence,ZZ,Ideal)    
   Headline
    finds a maximal regular sequence of minimal degree in an ideal
   Usage
    J=minimalRegularSequence(n,I)
   Inputs
    n:ZZ
    I:Ideal
   Outputs
    J:Ideal
   Description
    Text
    Example
   Caveat
   SeeAlso
///


end--
   
///
restart
loadPackage("ResidualIntersections", Reload =>true)
installPackage"ResidualIntersections"

viewHelp isLicci
S = ZZ/101[x_0..x_3]
installPackage "MCMApproximations"
I = ideal(x_0*x_1,x_1^2, x_2^3, x_3^5)
isLicci(3, codim I, I)
linkageBound (I, UseNormalModule => false)
linkageBound (I, UseNormalModule => true)

I = minors(2, random(S^2, S^{4:-1}))
isLicci(3, codim I, I)
linkageBound (I, UseNormalModule => false)
linkageBound (I, UseNormalModule => true)

I = minors(3, random(S^3, S^{-2,-3,-4,-4}));
isLicci(3, codim I, I)

linkageBound (I, UseNormalModule => false)
time linkageBound (I, UseNormalModule => true)


--b = linkageBound I

restart
loadPackage "ResidualIntersections"
loadPackage "RandomIdeal"
J = idealChainFromSC randomChain(10,5,20);
J/maxGd
J/residualCodims

///
