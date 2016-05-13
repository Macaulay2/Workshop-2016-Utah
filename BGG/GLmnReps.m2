newPackage(
	"GLmnReps",
    	Version => "1.0", 
    	Date => "May 8, 2016",
    	Authors => {
	     {Name => "Michael Perlman", Email => "mperlman@nd.edu", HomePage => "https://www3.nd.edu/~mperlman/"},
	     {Name => "Claudiu Raicu", Email => "craicu@nd.edu", HomePage => "https://www3.nd.edu/~craicu/"}
	     },
    	Headline => "Representations of gl(mn) and syzygies",
	PackageExports => {"BGG","PieriMaps","SchurRings","SimpleDoc"},
    	DebuggingMode => false
    	)
    
export{"characterLlambda", "characterKlambda", 
    "bettiTableIlambda", "conjecturalBettiTableIlambda", 
    "hilbertFunctionLlambda",
    "makeLinearComplexLGenericlambda"}

--the character of the exterior algebra
charE = method()
charE(ZZ,ZZ) := (m,n) -> (
    parsmn := flatten for i from 0 to m*n list select(partitions i,x->#x <= n and (try(x#0 <= m) else true));
    apply(parsmn,par -> {toList conjugate par,toList par})   
    )

--the character of the Kac module associated with lambda
characterKlambda = method()
characterKlambda(List,List,SchurRing,SchurRing) := (mn,ll,S,T) -> (
    cE := charE(mn#0,mn#1);
    chrE := sum apply(cE,x -> S_(x#0) * T_(x#1));
    S_ll * T_ll * chrE
    )
undocumented (characterKlambda,List,List,SchurRing,SchurRing)

characterKlambda(ZZ,ZZ,List) := (m,n,ll) -> (
    t := local t;
    s := local s;
    T := schurRing(t,n);
    S := schurRing(T,s,m);
    chrK := characterKlambda({m,n},ll,S,T);
    flatten apply(listForm(chrK),ter -> apply(listForm(ter#1), x -> {ter#0,x#0,x#1}))
    )

--half sum of the positive weights
delta := k -> for i from 0 to k-1 list k-1-i
delta = memoize delta

--the dot action of the symmetric group on partitions
dotact := (per,mu,k) -> (dk := delta k; (for i from 0 to k-1 list mu#(per#i) + dk#(per#i))-dk)
dotact = memoize dotact

--the length function for permutations
leng := (per,k) -> (lng := 0; for i from 0 to k-2 do for j from i+1 to k-1 do if per#i>per#j then lng = lng + 1;lng)
leng = memoize leng

generateMus := (lam,k) ->
(
    prs := {};
    for i from 0 to k-2 do
    (
    	noexit := true;
    	for j from i+1 to k-1 do if (lam#i - lam#j < j - i and noexit) then prs = prs | {(i,j)} else noexit = false;
    	);
    perms := select(permutations k,per -> (include := true; for c in prs do if per#(c#0)>per#(c#1) then include = false; include));
    perms = apply(perms, per -> apply(sort(apply(per,i->(per#i,i))),j->j#1));
    musleng := for per in perms list (pl := dotact(per,lam,k); lst := last pl; (leng(per,k),reverse(for i from 0 to k-1 list (lst = max(lst,pl#(k-1-i));lst))));
    musleng / (x->x#1)
    )

--the character of the irreducible glmn-module associated to the partition lambda
characterLlambda = method()
characterLlambda(ZZ,ZZ,List) := (m,n,lam) -> (
    lam = lam | splice{n-#lam:0};
    mus := generateMus(lam,n);
    MAX := lam + splice{n:m};
    rez := 0;
    s := local s;
    t := local t;
    S := schurRing(s,m);
    T := schurRing(S,t,n);
    for mu in mus do 
        for ll in mu..MAX do if ll == reverse sort(ll) then 
	    rez = rez + (-1)^(sum(ll-lam)) * characterKlambda({m,n},ll,S,T);
    selsmall := select(listForm rez,x->all(splice{0..#(x#0)-1},i-> (x#0#i <= MAX#i)));
    flatten apply(selsmall,x -> apply(listForm(x#1),ter -> {ter#0,x#0,ter#1}))
    )

--the betti table of the ideal I_lambda
--p is the prime characteristic of the residue field
bettiTableIlambda = method()
bettiTableIlambda(ZZ,ZZ,ZZ,List) := (n,m,p,lam) -> (
     r := local r;
     s := local s;
     x := local x;
     R := schurRing(r,n);
     T := schurRing(s,m);
     conjlam := toList conjugate( new Partition from lam);
     d := dim r_lam;
     e := dim s_lam;
     kk := ZZ/p;
     S := kk[x_(1,1)..x_(n,m)];
     M := genericMatrix(S,m,n);
     lis := for i from 0 to d*e-1 list
     (
    A := random(kk^m,kk^m);
    B := random(kk^n,kk^n);
    N := A * M * B;
    product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1})
     );
    J := ideal lis;
    I := ideal mingens J;
    if (numgens I != e*d) then error"wrong number of generators"
    else (F:= res(module I, FastNonminimal => true);betti(F, Minimize => true))
    )

bettiTableIlambda(ZZ,ZZ,List) := (n,m,lam) -> (
    bettiTableIlambda(n,m,32003,lam)
    )

--the dimension of the degree d component of the irreducible glmn module Llambda
dimLlamd = method()
dimLlamd(ZZ,ZZ,List,ZZ) := (n,m,lam,d) -> (
    F := select(characterLlambda(n,m,lam), x-> sum(x#0) == d);
    t := local t;
    s := local s;
    T1 := schurRing(t,n);
    T2 := schurRing(s,m);
    sum apply(F, x-> dim(T1_(x#0)) * dim(T2_(x#1)) * x#2)
    )

--the list of dimensions of graded components of the simple glmn module Llambda
hilbertFunctionLlambda = method()
hilbertFunctionLlambda(ZZ,ZZ,List) := (n,m,lam) -> (
    t := local t;
    s := local s;
    T1 := schurRing(t,n);
    T2 := schurRing(s,m);
    charL := characterLlambda(n,m,lam);
    G := apply(charL,x->sum(x#0));
    for i from sum lam to max G list (
	F := select(charL, x -> sum(x#0) == i);
	sum apply(F, x-> dim(T1_(x#0)) * dim(T2_(x#1)) * x#2)
	)
    )

--replaces x_i*y_j via e_(i,j) -- internal, not to be exported!
makeLinearMap := (pol,m,n,e,x,y) -> (
    sum flatten for i from 1 to m list 
                    for j from 1 to n list e_(i,j) * lift(contract(x_i*y_j,pol),coefficientRing ring pol)
    )
--make the linear complex associated via BGG to the E-module L_lambda
makeLinearComplexLGenericlambda = method()
makeLinearComplexLGenericlambda(ZZ,ZZ,List) := (m,n,lam) -> (
    kk := ZZ/32003;
    z := local z;
    S := kk[z_(1,1)..z_(m,n)];
    makeLinearComplexLGenericlambda(m,n,lam,S)
    )

makeLinearComplexLGenericlambda(ZZ,ZZ,List,PolynomialRing) := (m,n,lam,S) -> (
    kk := coefficientRing S;
    a := local a;
    b := local b;
    x := local x;
    y := local y;
    e := local e;
    A := kk[a_1..a_m];
    B := kk[b_1..b_n];
    C := kk[x_1..x_m,y_1..y_n];
    E := kk[e_(1,1)..e_(m,n),SkewCommutative => true];
    fAC := map(C,A,{x_1..x_m});
    fBC := map(C,B,{y_1..y_n});
    if #lam < n then lam = lam | {0};
    mu := {lam#0 + 1} | drop(lam,1);
    pA := pieri(mu,{1},A);
    pB := pieri(mu,{1},B);
    mat := fAC(pA) ** fBC(pB);
    for i from 1 to #lam-1 list if lam#(i-1) > lam#i then 
        (
	    ei := splice{i:0} | {1} | splice{#lam-i-1:0};
	    mu = lam + ei;
	    pA = pieri(mu,{i+1},A);
	    pB = pieri(mu,{i+1},B);
	    mat = mat | (fAC(pA) ** fBC(pB));
	    ) else continue;
    matE := matrix apply(entries mat,row -> apply(row,pol -> makeLinearMap(pol,m,n,e,x,y)));
    dual bggComplex(coker matE,S)
    )


-----Conjectural Betti table for I_lambda
-----

borderStrip = (n,lam) -> (
    b0 := (lam#0,0,lam#0);
    conjlam := conjugate new Partition from lam;
    listBoxes := {b0};
    for i from 1 to min(conjlam#0,n-1) do
    (
	li := for j from lam#i to lam#(i-1) list (j,i,i+j);
	listBoxes = listBoxes | (reverse li);
	);
    listBoxes
    )

admissiblePatterns := admissiblePatterns;

generateDyck = method()
generateDyck(ZZ,List,List,List) := (n,lamD,dyckPat,limits) ->
(
    starts := limits#0;
    ends := limits#1;
    admissiblePatterns = admissiblePatterns | {(lamD,set dyckPat)};
    bS := borderStrip(n,lamD); 
    for i from 0 to #bS-2 do if member((bS#i)_1,starts) then
    	for j from i+1 to #bS-1 do if (bS#j)_2 < (bS#i)_2 then break else
	    if member((bS#j)_1,ends) and (bS#i)_2 == (bS#j)_2 then
	    (
		newD := bS_(toList(i..j));
		ri := (bS#i)_1;
		rj := (bS#j)_1;
		newlamD := (for k from 0 to ri list if lamD#k>lamD#ri then lamD#k else lamD#k+1) |
		    	   (for k from ri+1 to rj list lamD#(k-1)+1) |
			   (for k from rj+1 to #lamD-1 list lamD#k);
		newstarts := starts - set(toList(ri+1..rj));
		newends := ends - set(toList(ri..rj-1));
		generateDyck(n,newlamD,dyckPat | {newD},{newstarts,newends});
		)
    )


conjecturalBettiTableIlambda = method()
conjecturalBettiTableIlambda(ZZ,ZZ,List) := (m,n,lam) ->
(
    l := set(toList(0..n-1));
    admissiblePatterns = {};
    if #lam < n then lam = lam | splice{n-#lam:0};
    generateDyck(n,lam,{},{l,l});
    bet := apply(toList(set admissiblePatterns),x-> (x#0,sum apply(toList(x#1),y->#y)));
    sum apply(bet,t -> new BettiTally from 
	(
	    hF := hilbertFunctionLlambda(m,n,t_0);
	    for i from 0 to #hF-1 list (t_1+i,{sum(t#0)+i},sum(t#0)+i) => hF#i)
	    )
    )


beginDocumentation()


doc ///
    Key
       GLmnReps
    Headline
       Representations of the general linear Lie superalgebra, and syzygies
    Description
      Text
        This package implements basic functionality for the representation theory of
	the general linear Lie superalgebra {\tt gl(m|n)}. We compute characters of irreducible 
	{\tt gl(m|n)}-modules, as well as those of Kac modules. As an application, we investigate
	the syzygies of some thickenings of determinantal ideals.
--    Caveat
--       It's far from completed
///
  
doc ///
    Key
       characterKlambda
       (characterKlambda,ZZ,ZZ,List)
    Headline
       Compute the character of the Kac module corresponding to a partition lambda
    Usage
       lis = characterKlambda(m,n,lambda)    
    Inputs
       m:ZZ
       n:ZZ
       lambda:List
    Outputs
       lis:List
           a list of triples, where the first two entries in each triple are partitions, and the third is a non-negative integer
    Description
      Text
        For a given a partition {\tt \lambda}, we let {\tt E=\Lambda(C^m\otimes C^n)} be the 
	exterior algebra and consider the Kac module 
	{\tt K_{\lambda} = S_{\lambda}C^m\otimes S_{\lambda}C^n\otimes E}. As an {\tt E}-module
	{\tt K_{\lambda}} is free, but it's structure as a {\tt gl(m|n)}-module is more
	interesting. A first approximation of this structure is its decomposition into irreducible
	representations of {\tt gl(m)\oplus gl(n)}. This is computed by the function
	{\tt characterKlambda} as a list of triples {\tt (\mu,\delta,r)}, where {\tt r} denotes
	the multiplicity of the irreducible {\tt gl(m)\oplus gl(n)}-representation
	{\tt S_{\mu}C^m\otimes S_{\delta}C^n} inside {\tt K_{\lambda}}.
	
	When {\tt \lambda=\{\}} is the empty partition, {\tt K_{\lambda} = E} and its character
	is described by Cauchy's formula:
      
      Example
    	c = characterKlambda(3,2,{})
	netList c --irreducibles are classified by pairs of conjugate partitions, and they appear with multiplicity one

      Text
      	In general, the partition {\tt \lambda} should have at most {\tt min(m,n)} parts.
	
      Example
        c = characterKlambda(2,2,{3,1})
	netList c
    
    SeeAlso
      characterLlambda
///  

doc ///
    Key
    	makeLinearComplexLGenericlambda
	(makeLinearComplexLGenericlambda,ZZ,ZZ,List)
	(makeLinearComplexLGenericlambda,ZZ,ZZ,List,PolynomialRing)
    Headline
    	The linear complex associated via the BGG correspondence to a (generic) simple gl(m|n)-module
    Usage
    	C = makeLinearComplexLGenericlambda(m,n,lambda,S)
    Inputs
    	m:ZZ
	n:ZZ
	lambda:List
	S:PolynomialRing
    Outputs
        C:ChainComplex
    Description
      Text
        Given a partition {\tt \lambda}, the irreducible {\tt gl(m|n)}-module {\tt L_{\lambda}} 
	is a module over the exterior algebra. Via the BGG correspondence, this module
	is transformed into a linear complex over the polynomial ring {\tt S} in {\tt m*n}
	variables. 
	
      Example
        m = 3;
	n = 2;
	lambda = {3,1};
        S = QQ[z_(1,1)..z_(m,n)];
	makeLinearComplexLGenericlambda(m,n,lambda,S)
	hilbertFunctionLlambda(m,n,lambda) --check to see that the ranks of the free modules in the complex are computed correctly

      Text
        When {\tt n=1} and {\tt \lambda=(d)} is a partition with only one part, {\tt L_{\lambda}}
	coincides with the minimal resolution of the ideal generated by degree d polynomials
	in {\tt m} independent variables.
      
      Example
        d = 5
	m = 3
        S = QQ[z_1..z_m]
	I = ideal(z_1..z_m)
	res module I^d
	makeLinearComplexLGenericlambda(m,1,{d},S)
	
    Caveat
       At the moment, this only produces the correct answer for generic enough partitions lambda.
///

doc ///
    Key
       characterLlambda
       (characterLlambda,ZZ,ZZ,List)
    Headline
       Compute the character of the irreducible gl(m|n)-module corresponding to a partition lambda
    Usage
       c = characterLlambda(m,n,lambda)    
    Inputs
       m:ZZ
       n:ZZ
       lambda:List
    Outputs
       c:List
         a list of triples, where the first two entries in each triple are partitions, 
	 and the third is a non-negative integer
    Description
      Text
        For a given a partition {\tt \lambda}, {\tt L_{\lambda}} is the irreducible 
	{\tt gl(m|n)}-module of lowest weight {\tt \lambda}. One concrete realization of {\tt L_{\lambda}}
	is as the unique simple quotient of the Kac module {\tt K_{\lambda}}. 
	The function {\tt characterLlambda} computes the decomposition of {\tt L_{\lambda}} into
	a direct sum of irreducible representations of {\tt gl(m)\oplus gl(n)}. The output is
	given as a list of triples {\tt (\mu,\delta,r)}, where {\tt r} denotes
	the multiplicity of the irreducible {\tt gl(m)\oplus gl(n)}-representation
	{\tt S_{\mu}C^m\otimes S_{\delta}C^n} inside {\tt L_{\lambda}}.
	
      Example
    	c = characterLlambda(3,1,{4})
	netList c
        d = characterLlambda(2,2,{3,1})
	netList d
    
    SeeAlso
      characterKlambda
/// 

doc ///
    Key
       bettiTableIlambda
       (bettiTableIlambda,ZZ,ZZ,List)	 
       (bettiTableIlambda,ZZ,ZZ,ZZ,List)
    Headline
        compute the Betti table of the GL_m X GL_n invariant ideal I_lambda
    Usage
        f = bettiTableIlambda(m,n,lam)
        f = bettiTableIlambda(m,n,p,lam)
    Inputs
        m:ZZ
        n:ZZ
        p:ZZ
	  a prime characteristic 
        lam:List 
    Outputs
        f:BettiTally
    Description
      Text
         The function takes as input three integers
         m,n,p (if the prime p is not provided, we take p=32003) and a partition {\tt lam}. 
	 It computes the Betti table of the {\tt GL_m\times GL_n}-invariant 
	 ideal {\tt I_{lam}} in characteristic p (assuming that p is generic enough). 
    
      Example
       m=3;
       n=3;
       p=101;
       lam={2,1};
       bettiTableIlambda(m,n,p,lam)
             
      Example
       m=3;
       n=2;
       lam={3};
       bettiTableIlambda(m,n,lam)
   
   Caveat
     This function only works in positive characteristic: it uses speed-ups for computing the
     Betti table which are currently not available in characteristic zero.
   SeeAlso
     conjecturalBettiTableIlambda
///

doc ///
    Key
       conjecturalBettiTableIlambda
       (conjecturalBettiTableIlambda,ZZ,ZZ,List)	 
    Headline
       compute the Betti table of the GL_m X GL_n invariant ideal I_lambda, based on a conjectural relationship with gl(m|n)-representations
    Usage
        b = conjecturalBettiTableIlambda(m,n,lam)
    Inputs
        m:ZZ
        n:ZZ
        lam:List 
    Outputs
        b:BettiTally
    Description
      Text
         The function takes as input the integers
         m,n and a partition {\tt lam}. 
	 It computes the betti table of the {\tt GL_m\times GL_n} invariant ideal {\tt I_{lam}}
    
      Example
       m = 3;
       n = 3;
       lam = {2,1};
       b = conjecturalBettiTableIlambda(m,n,lam)

      Example
       m=3;
       n=2;
       lam={3};
       conjecturalBettiTableIlambda(m,n,lam)

    SeeAlso
      bettiTableIlambda       
///

doc ///
    Key
       hilbertFunctionLlambda
       (hilbertFunctionLlambda,ZZ,ZZ,List)	 
    Headline
        The hilbert function of the simple gl(m|n)-module L_{lambda}
    Usage
        f = hilbertFunctionLlambda(m,n,lam)
    Inputs
        m:ZZ
        n:ZZ
        lam:List 
    Outputs
        f:List
    Description
      Text
         This function computes the hilbert function of the simple gl(m|n)-module L_{lambda} in
	 all degrees for which it is not zero. 
    
      Example
       m=2;
       n=2;
       lam={2,1};
       hilbertFunctionLlambda(m,n,lam)
      
      Text
       Another example:
       
      Example
       m=3;
       n=2;
       lam={2,2};
       hilbertFunctionLlambda(m,n,lam)
///
 
end

restart
uninstallPackage"GLmnReps"
installPackage"GLmnReps"
viewHelp GLmnReps

debug needsPackage"GLmnReps"

characterLlambda(3,2,{2,1})
characterKlambda(3,2,{1})
charE(4,3)
hilbertFunctionLlambda(3,3,{3,1})
dimLlamd(3,3,{3,1},5)
time conjecturalBettiTableIlambda(4,3,{2,1,0})
time bettiTableIlambda(4,3,{2,1,0})

--TODO
--list the composition factors of K_lambda
--construct the linear complex corresponding to L_lambda for arbitrary lambda




