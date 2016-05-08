newPackage ( "ResidualIntersections",
    Version => "1.0",
    Date => "07 May 2016",
    Authors => {
	{Name => "David Eisenbud",
         Email => "de@msri.org",
         HomePage => "http://www.msri.org/~de"},
     	 {Name => "Robert,Katy,Robert, Jay"}
	},
    Headline => "Package for studying conditions associated to Residual Intersection theory"
    Reload => true,
    DebuggingMode => true
    )


export {
	"isLicci",
	"minimalRegularSequence",
	"linkageBound",
	"UseNormalModule",
	"randomRegularSequence"
        };

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

///
