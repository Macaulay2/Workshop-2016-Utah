newPackage ( "ResidualIntersections",
    Version => "1.0",
    Date => "07 May 2016",
    Authors => {
	{Name => "Katie Ansaldi",
	    Email => "kansaldi@gmail.com"},
	{Name => "David Eisenbud",
	    Email => "de@msri.org",
	    HomePage => "http://www.msri.org/~de"},
	{Name => "Robert Krone",
	    Email => "rckrone@gmail.com",
	    HomePage => "http://rckr.one"},
	{Name => "Jay Yang",
	    Email => "jkelleyy@gmail.com"}
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
	"genericResidual",
	"genericArtinNagata",
	"numgensByCodim",
	"maxGs",
	"residualCodims",
        "koszulDepth",
        "hasSlidingDepth",
        "isStronglyCM",
	"depthsOfPowers"
        };

--
depthsOfPowers = method()
depthsOfPowers(ZZ,ZZ,Ideal) := (s,c,I) ->(
    --c should be codim I
    S := ring I;
    apply(s-c+1, j->profondeur(S^1/I^(j+1)))
    )
depthsOfPowers(ZZ,Ideal) := (s,I) -> depthsOfPowers(s,codim I, I)


--generic Artin-Nagata Code
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
codim I
depthsOfPowers(numgens ring I,I)

S = ZZ/32003[x_0..x_3]
I = minors(3, random(S^3, S^{-2,-3,-3,-3}));
--codim 3
codim I
s = 2;
depthsOfPowers(3,I)
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
    max(0, 2*(codim I)*(degree I -1) -6)
 else (
	N := prune Hom(I, (ring I)^1/I);
	max(0, 2*(numgens N - codim I)-6))
)

{*
minimalRegularSequence = method()
minimalRegularSequence(ZZ,Ideal) := (c,I) ->(
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
*}

minimalRegularSequence = method()
minimalRegularSequence(ZZ,Ideal) := (c,I) ->(
if numgens I == c then return I;
    --takes care of I = 0 and I principal;
sgens := sort gens I;
n :=numcols sgens;
J := ideal sgens_{0};
K := J;
c' := 1;
c'' := c;
for i from 1 to n-1 do(
    c'' = codim(K = J + ideal(sgens_{i}));
    if c''>c' then (
        J = K;
	c' = c''));
if c' == c then return J;
rgens := sgens * random(source sgens, source sgens);
for i from 0 to n-1 do(
    c'' = codim(K = J + ideal(rgens_{i}));
    if c''>c' then(
        J = K;
	c' = c'';
	if c' ==c then break));
J)
minimalRegularSequence Ideal := I -> minimalRegularSequence(codim I, I)

///
restart
loadPackage "ResidualIntersections"
S = ZZ/101[a,b,c]
I = ideal"cb,b2,ab,a2"
codim I 
minimalRegularSequence(codim I, I)
     I = ideal"cb,b2,a2"
     minimalRegularSequence I     
     I = ideal"ab,ac,bc"
     minimalRegularSequence( I)
///

isLicci = method(Options => {UseNormalModule =>false, Verbose =>false})
isLicci(ZZ, ZZ, Ideal) := opts -> (b,c,I) -> (
    --I homogeneous ideal
    --b = linkageBound I
    --c = codim I
    --output is list of up to b integers, the numbers of generators of the
    --successive random links
    J := I;
    p := numgens J;
    count := 0;
--    <<p<<endl;flush;
    scan(b, i -> (count = i+1;
	    J = randomLink(c,J);
	    if numgens J == c then break;
	    if Verbose === true then <<numgens J<<endl<<flush;
	    ));
   if Verbose === true then <<" done in "<< count << " steps"<<endl;
   c == numgens J)

isLicci(ZZ,Ideal) := opts -> (b,I) -> isLicci(b,codim I, I)
isLicci Ideal := opts -> I -> (
isLicci(linkageBound(I, UseNormalModule => opts.UseNormalModule), I
    ))
///
restart
loadPackage("ResidualIntersections", Reload=>true)
installPackage("RandomIdeal")
--viewHelp RandomIdeal
     setRandomSeed 0     
     S = ZZ/32003[x_0..x_6]
     L = idealChainFromShelling(S,randomShelling(7,3,6))
     
     apply(L, I-> {linkageBound I, linkageBound(I, UseNormalModule =>true)})
     scan(L, I ->print isLicci(I, UseNormalModule => true))
     numgens prune Hom(module I, S^1/I)
     
     (codim I)*(degree I)
///
--depth but faster
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
    if I==0 then return {};
    C := koszul mingens I;
    for i in 0..(numColumns(mingens I)-codim I) list profondeur HH_i(C)
    )

koszulDepth(ZZ,Ideal) := (k,I) -> (
    if I==0 then return -1;
    C := koszul mingens I;
    profondeur HH_k(C)
    )


isStronglyCM = method()
isStronglyCM(Ideal) := I -> (
    d := dim I;
    all(koszulDepth I,i -> i==d)
    )

hasSlidingDepth = method()

hasSlidingDepth(ZZ,Ideal) := (k,I) -> (
    d := dim I;
    s := numColumns(mingens I)-codim I;
    all(k+1, i -> (koszulDepth(s-i,I))>=d-i)
    )

hasSlidingDepth(Ideal) := I -> (
    s := numColumns(mingens I)-codim I;
    hasSlidingDepth(s,I)
    )

-------------------------------------
-- G_s Code
-------------------------------------

numgensByCodim = method()	
numgensByCodim (MonomialIdeal,ZZ) := (J,k) -> (
    R := ring J;
    n := numgens R;
    max for A in subsets(n,k) list (
	M := new MutableList from (n:1_R);
	for a in A do M#a = R_a;
	M = map(R,R,matrix{toList M});
	numgens trim M J
	)
    )

numgensByCodim MonomialIdeal := J -> (
    n := numgens ring J;
    toList apply(n, i->numgensByCodim(J,i+1))
    )

maxGs = method()
maxGs MonomialIdeal := J -> (
    for i from 1 to numgens ring J do if numgensByCodim(J,i) > i then return i;
    infinity
    )

residualCodims = method()
residualCodims MonomialIdeal := J -> (
    toList select((codim J + 1..numgens ring J + 1), i->numgensByCodim(J,i-1) <= i)
    )

------------------------------------------------------------
-- DOCUMENTATION ResidualIntersections
------------------------------------------------------------
doc ///
   Key
    ResidualIntersections
   Headline
    Tests for the conditions used in the theory of residual intersections
   Description
    Definition: If I \subset S is an ideal in a polynomial ring (or Gorenstein ring) and
    a_1..a_s are elements of I, then K = (a_1..a_s):I is called an
    s-residual intersection of I if the codimension of K is at least s.
    
    In the simplest case, s == codim I, the ideal K is said to be linked to I
    if also I = (a_1..a_s):K; this is automatic when S/I is Cohen-Macaulay,
    and in this case S/K is also Cohen-Macaulay; see Peskine-Szpiro,
    Liaison des variétés algébriques. I. Invent. Math. 26 (1974), 271–302).

    The theory for s>c, which has been used in algebraic geometry since the 19th century,
    was initiated in a commutative algebra setting by Artin and Nagata in the paper
    Residual intersections in Cohen-Macaulay rings. 
    J. Math. Kyoto Univ. 12 (1972), 307–323.
    
    Craig Huneke (Strongly Cohen-Macaulay schemes and residual intersections,
    Trans. Amer. Math. Soc. 277 (1983), no. 2, 739–763)
    proved that an s-residual intersection K is Cohen-Macaulay
    if I satisfies the G_d condition and
    is strongly Cohen-Macaulay, and successive authors have weakened the latter
    condition to sliding depth, and, most recently, Bernd Ulrich
    (Artin-Nagata properties and reductions of ideals. 
    Commutative algebra: syzygies, multiplicities, and 
    birational algebra,
    Contemp. Math., 159, 1994) showed that
    the weaker condition
    depth( S/(I^t) ) >= dim(S/I) - (t-1) for t = 1..s-codim I +1
    suffices. All these properties are true if I is licci.
    
    This package implements tests for most of these properties.
   SeeAlso
///

------------------------------------------------------------
-- DOCUMENTATION isLicci
------------------------------------------------------------
doc ///
   Key
    isLicci
    (isLicci,ZZ,ZZ,Ideal)
    (isLicci,ZZ,Ideal)
    (isLicci,Ideal)
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
     
     Every perfect codimension 2 ideal (nxn minors of an (nx(n+1) matrix) is licci,
     but other ideals of minors are generally not, as illustrated below.
    Example
     setRandomSeed 0     
     needsPackage "RandomIdeal"
     needsPackage "ResidualIntersections"
     S = ZZ/32003[x_0..x_6]
     L = idealChainFromShelling(S,randomShelling(7,3,8))
     apply(L, I-> {linkageBound I, linkageBound(I, UseNormalModule =>true)})
     scan(L, I ->print isLicci(I, UseNormalModule => true))
   Caveat
    linkageBound I can be very large; linkageBound(I, UseNormalModule => true) can be slow.
   SeeAlso
    linkageBound
///
------------------------------------------------------------
-- DOCUMENTATION UseNormalModule
------------------------------------------------------------
doc ///
   Key
    UseNormalModule
    [isLicci,UseNormalModule]
    [linkageBound, UseNormalModule]
   Headline
    option for linkageBound and isLicci
   Description
    Text
     Default value is false. When true, it
     enables a more refined computation of the bound on the number of general links of
     an ideal I
     that must be taken to definitively test the licci propery. 
     When UseNormalModule == true the computation of the 
     normal module Hom(I, (ring I)/I) is required and this can be slow;
     if UseNormalModule == false the computation is fast, but the bound is large.
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
     An ideal I in a polynomial ring S is licci if it Cohen-Macaulay and is linked in finitely many steps
     I --> (F):I, where F is a maximal regular sequence in I,
     to a complete intersection. Bernd Ulrich showed that if I is licci and each
     step of the linkage
     is done via a regular sequence F that is a subset of a minimal set of generators,
     then the linkage process will terminate after at most b steps, where
     
     b = 2(codim I)*(degree I -1) -6.
     
     (Theorem 2.4 of "On Licci Ideals", Contemp. Math 88 (1989).
     This is computed by linkageBound I.
     He did this via a more refined formula; the (generally sharper) 
     intermediate result gives the bound 
     
     b = 2(numgens(Hom(I, S/I) - codim I).
	 
     The call linkageBound(I, UseNormalModule =>true) computes this refined bound.
     See isLicci for examples.
   Caveat
    The crude bound can be quite large; computing the refined bound (which is often large
    as well) can be quite slow.
   SeeAlso
    isLicci
    UseNormalModule
///


------------------------------------------------------------
-- DOCUMENTATION minimalRegularSequence
------------------------------------------------------------
doc ///
   Key
    minimalRegularSequence
    (minimalRegularSequence,ZZ,Ideal)
    (minimalRegularSequence,Ideal)
   Headline
    finds a maximal regular sequence of minimal degree in an ideal
   Usage
    J=minimalRegularSequence(n,I)
    J=minimalRegularSequence(I)
   Inputs
    n:ZZ
    I:Ideal
///

------------------------------------------------------------
-- DOCUMENTATION maxGs
------------------------------------------------------------

doc ///
   Key
      maxGs
      (maxGs,MonomialIdeal)    
   Headline
      maximum G_s of a monomial ideal
   Usage
      d = maxGs I
   Inputs
      I:MonomialIdeal
   Outputs
      d:ZZ
         the maximum value of {\tt s} such that {\tt I} has property G_s (possibly infinity).
   Description
      Text
      Example
      	  R = QQ[x_1,x_2,x_3];
	  I = monomialIdeal(x_1^2,x_1*x_2,x_1*x_3,x_2^2,x_2*x_3);
	  maxGs(I)
   Caveat
   SeeAlso
      numgensByCodim
      residualCodims
///

------------------------------------------------------------
-- DOCUMENTATION genericArtinNagata
------------------------------------------------------------

doc ///
   Key
    genericArtinNagata
    (genericArtinNagata,ZZ,Ideal)    
   Headline
    Generic Artin nagata
   Usage
    L = genericArtinNagata(n,I)
   Inputs
    n:ZZ
    I:Ideal
   Outputs
    L:List
   Description
    Text
     Produces an "efficient" regular sequence that is among a set of minimal generators of I
     by the following algorithm
     Sorts the generators of I by degree-rev-lex to get sgens, and
     takes as many elements from this list as possible. Then 
     makes a general triangular change of generators to get rgens.
     and take the rest of the regular sequence from this list.
    Example
     setRandomSeed 0
     S = ZZ/101[a,b,c]
     I = ideal"ab,b2,c2"
     minimalRegularSequence I
     I = ideal"cb,b2,a2"
     minimalRegularSequence I     
   Caveat
   SeeAlso
///

------------------------------------------------------------
-- DOCUMENTATION numgensByCodim
------------------------------------------------------------
doc ///
   Key
      numgensByCodim
      (numgensByCodim,MonomialIdeal)
      (numgensByCodim,MonomialIdeal,ZZ)
   Headline
      maximum number of generators of localizations of a monomial ideal
   Usage
      d = numgensByCodim(I,k)
      L = numgensByCodim(I)
   Inputs
      I:MonomialIdeal
      k:ZZ
         an integer between 1 and the dimension of the ring
   Outputs
      d:ZZ
         the maximum number of generators of {\tt I} localized at a prime {\tt P} of codimension {\tt k}.
      L:List
         a list of the numbers of generators for each codimension from 1 to the dimension of the ring
   Description
      Text
         Because {\tt I} is monomial, we can check the number of generators of {\tt I} localized at a prime {\tt P} over only monomial primes {\tt P}.
      Example
         R = QQ[x_0..x_4];
	 I = monomialIdeal{x_0^2,x_1*x_2,x_3*x_4^2}
	 numgensByCodim(I,2)
	 numgensByCodim I
   SeeAlso
      residualCodims
      maxGs
///

------------------------------------------------------------
-- DOCUMENTATION residualCodims
------------------------------------------------------------

doc ///
   Key
      residualCodims
      (residualCodims,MonomialIdeal)
   Headline
      a list of possible residual intersection codimensions
   Usage
      L = residualCodims I
   Inputs
      I:MonomialIdeal
   Outputs
      L:List
         a list of integers {\tt s} such that {\tt I} localized at any prime of codimension {\tt s-1} has at most s generators.
	 The range of values is {\tt codim I} + 1 and the dimension of the ring + 1.
   Description
      Text
         For each {\tt s} computes the maximum over all monomial primes {\tt P} with codimension {\tt s-1} 
	 of the minimal size of a generating set of {\tt I} localized at {\tt P}.  If this number is 
	 less than {\tt s}, then {\tt s} is included in the list.
      Text
         The values {\tt s} returned are candidates for {\tt I} possibly being an {\tt s}-rsidual intersection.
      Example
         R = QQ[a,b,c];
	 I = monomialIdeal{a*b,b*c^2}
	 residualCodims I
   SeeAlso
      numgensByCodim
      maxGs
///

------------------------------------------------------------
-- DOCUMENTATION isStronglyCM
------------------------------------------------------------

doc ///
   Key
      isStronglyCM
      (isStronglyCM,Ideal)
   Headline
      Checks if the given ideal is Strongly Cohen Macaulay
   Usage
      b = isStronglyCM I
   Inputs
      I:Ideal
   Outputs
      b:Boolean
         true if {\tt I} is Strongly Cohen Macaulay
   Description
      Text
         Checks whether {\tt I} is Strongly Cohen Macaulay. We compute the depths of the Koszul homology by using {\tt koszulDepth} and compares it to {\tt codim I}.
      Example
         R = QQ[x_1..x_5];
	 I = ideal{x_1*x_3,x_2*x_4,x_3*x_4,x_1*x_5,x_3*x_5};
         isStronglyCM I
   SeeAlso
       koszulDepth
       hasSlidingDepth
///

------------------------------------------------------------
-- DOCUMENTATION koszulDepth
------------------------------------------------------------

doc ///
   Key
      koszulDepth
      (koszulDepth,Ideal)
      (koszulDepth,ZZ,Ideal)
   Headline
      Computes the depths of the Koszul homology
   Usage
      L = koszulDepth I
      d = koszulDepth(k,I)
   Inputs
      I:Ideal
      k:ZZ
         the homological index to compute
   Outputs
      L:List
         a list of the depths of Kozul homology
      d:ZZ
         the depth of the k-th Koszul homology
   Description
      Text
         The one parameter version computes the depths of the non-vanishing Koszul homology of {\tt I}.
         The two parameter version computes only the depth of the {\tt k}-th Koszul homology.
      Example
         R = QQ[x_1..x_6];
	 I = ideal{x_1*x_2,x_1*x_3,x_2*x_4*x_5,x_1*x_6,x_4*x_6,x_5*x_6};
         koszulDepth I
         koszulDepth(2,I)
   SeeAlso
       isStronglyCM
       hasSlidingDepth
///

------------------------------------------------------------
-- DOCUMENTATION hasSlidingDepth
------------------------------------------------------------

doc ///
   Key
      hasSlidingDepth
      (hasSlidingDepth,Ideal)
      (hasSlidingDepth,ZZ,Ideal)
   Headline
      Checks if an ideal has the sliding depth property
   Usage
      b = hasSlidingDepth I
      b = hasSlidingDepth(k,I)
   Inputs
      I:Ideal
      k:ZZ
   Outputs
      b:Boolean
         true if {\tt I} has sliding depth
   Description
      Text
         This computes whether the ideal {\tt I} has sliding depth.
      Text
         For an ideal $I$ with minimal generating set ${\bf f}=(f_1,\ldots,f_n)$, we say $I$
         has k-sliding depth if for all $i\leq k$ we have $depth(H_{n-codim(I)-i}({\bf f}))\geq dim I - i$.
         Note that since $H_{n-codim(I)}({\bf f})$ is the canonical module which always has
         depth equal to $dim I$, every ideal has 0-sliding depth. We say that a module has
         sliding depth without a parameter if it has $(n-codim(I))$-sliding depth
      Example
         R = QQ[x_1..x_6];
	 I = ideal{x_1*x_2,x_1*x_3,x_2*x_4*x_5,x_1*x_6,x_4*x_6,x_5*x_6};
         hasSlidingDepth I
         hasSlidingDepth(1,I)
         hasSlidingDepth(2,I)
   SeeAlso
       isStronglyCM
       koszulDepth
///

------------------------------------------------------------
-- DOCUMENTATION depthsOfPowers
------------------------------------------------------------

doc ///
   Key
      depthsOfPowers
      (depthsOfPowers,ZZ,ZZ,Ideal)
      (depthsOfPowers,ZZ,Ideal)
   Headline
      Computes depth of powers of an ideal
   Usage
      L = depthsOfPowers(s,c,I)
      L = depthsOfPowers(s,I)
   Inputs
      s:ZZ
      	  number of powers to compute
      c:ZZ
      	  (If omitted, it will use c = codim I)
      I:Ideal
   Outputs
      L:List
      	 The depths of the powers of I from 1 to s-c+1. 
   Description
      Text
      	  Computes the depth of $S/I^k$ for $k$ from 1 to $s-c+1$.
      Example
      	  R = QQ[a,b,c,d,e,f];
	  I = ideal (b*c, b*d, b*e, d*e, a*d*f, e*f);
	  depthsOfPowers(6,3,I)
   Caveat
   SeeAlso
///


------------------------------------------------------------
-- DOCUMENTATION genericResidual
------------------------------------------------------------

doc ///
   Key
      genericResidual
      (genericResidual,ZZ,Ideal)
   Headline
      Computes generic residual intersections of an ideal
   Usage
      J = genericResidual(s,I)
   Inputs
      s:ZZ
      I:Ideal
   Outputs
      J:Ideal
   Description
      Text
      Example
   Caveat
   SeeAlso
///

TEST ///
R = QQ[x_1..x_6];
I = ideal{x_1*x_2,x_1*x_3,x_2*x_4*x_5,x_1*x_6,x_4*x_6,x_5*x_6};
assert(koszulDepth I == {3,0,2,3})
assert(koszulDepth(2,I) == 2)
assert(not hasSlidingDepth I)
assert(hasSlidingDepth(1,I))
assert(not hasSlidingDepth(2,I))
assert(not isStronglyCM I)
///

TEST ///
R = QQ[x_1..x_5];
I = ideal{x_1*x_3,x_2*x_4,x_3*x_4,x_1*x_5,x_3*x_5};
assert(isStronglyCM I)
///
end--

linkageBound (I, UseNormalModule => false)
time linkageBound (I, UseNormalModule => true)


--b = linkageBound I

restart
loadPackage "ResidualIntersections"
loadPackage "RandomIdeal"
uninstallPackage "ResidualIntersections"
installPackage "ResidualIntersections"

J = idealChainFromSC randomChain(10,5,20);
J/maxGs
J/residualCodims
