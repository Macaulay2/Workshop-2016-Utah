--Approximates the F-pure Threshold
--Gives a list of nu_I(p^d)/p^d for d=1,...,e
FPTApproxList = method();

FPTApproxList (Ideal,ZZ) := (I,e) ->
(
     p := char ring I;
     nus := nuList(I,e);
     apply( nus, 1..e, (n,k) -> n/p^k )
)

FPTApproxList (RingElement,ZZ) := (f,e) -> FPTApproxList(ideal(f),e)

--Approximates the F-Threshold with respect to an ideal J
--More specifically, this gives a list of nu_I^J(p^d)/p^d for d=1,...,e

FTApproxList = method();

FTApproxList (Ideal,Ideal,ZZ) := (I,J,e) ->
(
    if not isSubset( I, radical(J) ) then error "F-threshold undefined.";
     p := char ring I;
     nus := nuList(I,J,e);
     apply( nus, 1..e, (n,k) -> n/p^k )
)

FTApproxList (RingElement,Ideal,ZZ) := (f1,J1,e1) -> FTApproxList(ideal(f1),J1,e1)

FTHatApproxList = method();

FTHatApproxList (Ideal,Ideal,ZZ) := (I,J,e) ->
(
    if not isSubset( I, radical(J) ) then error "F-threshold undefined.";
     p := char ring I;
     nus := nuHatList(I,J,e);
     apply( nus, 1..e, (n,k) -> n/p^k )
)

FTHatApproxList (RingElement,Ideal,ZZ) := (f1,J1,e1) -> FTHatApproxList(ideal(f1),J1,e1)


--Guesses the FPT of ff.  It returns a list of all numbers in 
--the range suggested by nu(ff,e1) with maxDenom as the maximum denominator
guessFPT ={OutputRange=>false}>>o -> (ff, e1, maxDenom) ->(
     nn := nu(ff, e1);
     pp := char ring ff;
     if (o.OutputRange == false) then 
          findNumberBetween({nn/(pp^e1-1), (nn+1)/(pp^e1)}, maxDenom)
     else
          {{ nn/(pp^e1-1), (nn+1)/(pp^e1)}, findNumberBetween({nn/(pp^e1-1), (nn+1)/(pp^e1)}, maxDenom)}
)

--F-pure threshold estimation, at the origin
--e is the max depth to search in
--FinalCheck is whether the last isFRegularPoly is run (it is possibly very slow) 
--If MultiThread is set to true, it will compute the two F-signatures simultaneously
estFPT={FinalCheck=> true, Verbose=> false, MultiThread=>false, DiagonalCheck=>true, BinomialCheck=>true, NuCheck=>true} >> o -> (ff,ee)->(
     print "starting estFPT";
     
     maxIdeal := ideal( first entries vars( ring ff) );   --the maximal ideal we are computing the fpt at  

     foundAnswer := false; --this will be set to true as soon as we found the answer.  Setting it to true will stop further tests from being run
     answer := null; --this stores the answer until it can be returned.
     
     --first check if it is diagonal:
     if ( (o.DiagonalCheck==true) and (foundAnswer == false) ) then (
	   if (isDiagonal(ff)==true) then ( 
		if (o.Verbose==true) then print "Polynomial is diagonal."; 
		answer = diagonalFPT(ff); 
		foundAnswer = true
	   )
     );

     --now check if it is binomial:
     if ( (o.BinomialCheck==true) and (foundAnswer == false) ) then (
	  if (isBinomial(ff)==true) then ( 
	       if  (o.Verbose==true) then print "Polynomial is binomial.";
	       answer = binomialFPT(ff);
	       foundAnswer = true
	  )
     );
     
     --compute nu's
     if (foundAnswer == false) then (
     	  pp:=char ring ff;
     	  nn:=nu(ff,ee);
	  if  (o.Verbose==true) then print "nu's have been computed";

     	  --if our nu's aren't fine enough, we just spit back some information
       	  if nn==0 then (
	       answer = {0,1/pp};
	       foundAnswer = true
	   )
      );
 
      --check to see if nu/(p^e-1) is the fpt
      if ((o.NuCheck==true) and (foundAnswer == false)) then (
	   if (isFRegularPoly(ff,(nn/(pp^ee-1)),maxIdeal)==false) then ( 
		if  (o.Verbose==true) then print "Found answer via nu/(p^e-1)."; 
		answer = nn/(pp^ee-1);
		foundAnswer = true
	   ) 
      	   else (
	   	if  (o.Verbose==true) then print "nu/(p^e - 1) is not the fpt.";
	   )
      );
	 
	--check to see if (nu+1)/p^e is the FPT
	if ((o.NuCheck==true) and (foundAnswer == false)) then(
		if (isFPTPoly(ff, (nn+1)/pp^ee,Origin=>true) == true) then (
			answer = (nn+1)/pp^ee;
			foundAnswer = true
		)
	);

     --do the F-signature computation
     if (foundAnswer == false) then (
	   ak := 0;
	   if (o.MultiThread==false ) then (ak=threshInt(ff,ee,(nn-1)/pp^ee,fSig(ff,nn-1,ee),nn) ) else(
		if (o.Verbose==true) then print "Beginning multithreaded F-signature";
		allowableThreads = 4;
		numVars := rank source vars (ring ff);
		YY := local YY;
		myMon := monoid[  toList(YY_1..YY_numVars), MonomialOrder=>RevLex,Global=>false];
		--print myMon;
     		R1:=(coefficientRing ring ff) myMon;
		rMap := map(R1, ring ff, vars myMon);
		gg := rMap(ff);
		
		
		H := (fff,aaa,eee) -> () -> fSig(fff,aaa,eee);
		newSig1 := H(gg,nn-1,ee);
		t1 := schedule newSig1;
	     	s2 := fSig(ff,nn,ee);	
		if (o.Verbose==true) then print "One signature down";
		while ((not isReady t1)) do sleep 1;
		s1 := taskResult t1;
     	      --  print s1; print s2;
		ak = {s2,xInt( (nn-1)/pp^ee, s1, (nn)/pp^ee,s2)};
		--print nn;		
	   );
	   if  (o.Verbose==true) then print "Computed F-signatures.";
	   --now check to see if we cross at (nu+1)/p^e, if that's the case, then that's the fpt.
	   if ( (nn+1)/pp^ee == (ak#1) ) then (
		if  (o.Verbose==true) then print "F-signature line crosses at (nu+1)/p^e."; 
		answer = ak#1;
		foundAnswer = true
	   )
      );	  
      	
      --if we run the final check, do the following
      if ( (foundAnswer == false) and (o.FinalCheck == true)) then ( 
	  if  (o.Verbose==true) then print "Starting FinalCheck."; 
          	if ((isFRegularPoly(ff,(ak#1),maxIdeal)) ==false ) then (	
	      		if  (o.Verbose==true) then print "FinalCheck successful"; 
	      		answer = (ak#1);
	      		foundAnswer = true 
      	  	)
	  		else ( 
	      		if  (o.Verbose==true) then print "FinalCheck didn't find the fpt."; 
	      		answer = {(ak#1),(nn+1)/pp^ee};
	      		foundAnswer = true
	  		)
       );
       
       --if we don't run the final check, do the following
       if ((foundAnswer == false) and (o.FinalCheck == false) ) then (
	  if  (o.Verbose==true) then print "FinalCheck not run.";
	  answer = {(ak#1),(nn+1)/pp^ee};
      	  foundAnswer = true
       );
     
     --return the answer
     answer
)

--isFPTPoly, determines if a given rational number is the FPT of a pair in a polynomial ring. 
--if Origin is specified, it only checks at the origin. 

isFPTPoly ={Verbose=> false,Origin=>false}>> o -> (f1, t1) -> (
	pp := char ring f1;
	if (o.Origin == true) then org := ideal(vars (ring f1));
	funList := divideFraction(t1, pp);
	--this writes t1 = a/(p^b(p^c-1))
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	mySigma := ideal(f1);
	myTau := ideal(sub(1, ring f1));
	myA := aa;
	myA2 := 0;
	
	if (cc != 0) then (
		myA = floor(aa / (pp^cc - 1));
		myTau = tauPoly( f1, (aa%(pp^cc-1))/(pp^cc-1) )
	);
	
	if (o.Verbose==true) then print "higher tau Computed";

	--first we check whether this is even a jumping number.
	if (cc == 0) then (
		myA2 = aa-1;
		mySigma = sigmaAOverPEMinus1Poly(f1, (pp-1), 1)
	)
	else (
		myA2 = floor((aa-1)/(pp^cc-1));
		mySigma = (sigmaAOverPEMinus1Poly(f1, ((aa-1)%(pp^cc-1))+1, cc))
	);
	if (o.Verbose==true) then print "higher sigma Computed";

	returnValue := false;
	
	if ( isSubset(ideal(sub(1, ring f1)), ethRootSafe(f1, mySigma, myA2, bb) )) then (
		if (o.Verbose==true) then print "we know t1 <= FPT";
		if (not isSubset(ideal(sub(1, ring f1)), ethRootSafe(f1, myTau, myA, bb) ))  then returnValue = true 
	);
		
	returnValue
)

--isFJumpingNumberPoly determines if a given rational number is an F-jumping number
--***************************************************************************
--This needs to be speeded up, like the above function
--***************************************************************************
isFJumpingNumberPoly ={Verbose=> false}>> o -> (f1, t1) -> (
	pp := char ring f1;
	funList := divideFraction(t1, pp);
	--this writes t1 = a/(p^b(p^c-1))
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	mySigma := ideal(f1);
	myTau := ethRoot(tauPoly(f1, t1*pp^bb), bb);
	if (o.Verbose==true) then print "higher tau Computed";

	--first we check whether this is even a jumping number.
	if (cc == 0) then
		mySigma = ethRoot((ideal(f1^(aa-1)))*((sigmaAOverPEMinus1Poly(f1, (pp-1), 1))), bb)
	else 
		mySigma = ethRoot((sigmaAOverPEMinus1Poly(f1, aa, cc)),bb);
	if (o.Verbose==true) then print "sigma Computed";

	not (isSubset(mySigma, myTau))
)


--Computes the F-pure threshold of a diagonal hypersurface 
--x_1^(a_1) + ... +x_n^(a_n) using Daniel Hernandez' algorithm
diagonalFPT = f ->
(
     p := char ring f;
     w := apply(terms f, g->first degree(g));
     y := 0; if firstCarry(reciprocal(w),p)==-1 then for i from 0 to #w-1 do y = y + 1/w#i else
     (
	  x := 0; for c from 0 to #w-1 do x = x + truncation(firstCarry(reciprocal(w),p)-1, 1/w#c, p); 
	  y = x+1/p^(firstCarry(reciprocal(w),p)-1);
     );
     y
)

--Given a polynomial f, outputs a list of multi-degrees (under the usual grading)
--of f as lists
multiDegree = f ->
(
     variables := first entries vars ring f;
     apply(terms f, g -> apply(#variables, i ->degree(variables#i,g)))
)

--Determines whether a polynomial f is diagonal; i.e., of the form 
--x_1^(a_1)+...+x_n^(a_n) (up to renumbering variables)
isDiagonal = f ->
(
     d := multiDegree(f);
     alert1 := true;
     alert2 := true;
     for i from 0 to #d-1 do
     (
	  for j from 0 to #(d#0)-1 do
	  (
	       if (d#i)#j!=0 and alert1==false then alert2=false;
	       if (d#i)#j!=0 and alert1==true then alert1=false;
	  );
     alert1=true;
     );
     for j from 0 to #(d#0)-1 do
     (
	  for i from 0 to #d-1 do 
	  (
     	       if alert1==false and (d#i)#j!=0 then alert2=false;
     	       if alert1==true and (d#i)#j!=0 then alert1=false;
	  );
     alert1=true;
     );
     alert2
)

--Given input vectors v={a_1,...,a_n} and w={b_1,...,b_n}, gives the
--corresponding vectors that omit all a_i and b_i such that a_i=b_i
factorOutMonomial = (v,w) ->
(
     v1 := new MutableList;
     w1 := new MutableList;
     c := 0; i := 0; for i from 0 to #v-1 do (if v#i != w#i then (v1#c = v#i; w1#c = w#i; c = c+1; ); );
     (v1,w1)
)

--Given input vectors v={a_1,...,a_n} and w={b_1,...,b_n}, gives the
--vector of the a_i for which a_i=b_i
monomialFactor = (v,w) ->
(
     a := new MutableList;
     c := 0; i := 0; for i from 0 to #v-1 do (if v#i == w#i then (a#c = v#i; c = c+1; ); );
     a
)

--Given two vectors v={v0,v1} and w={w0,w1} in the real plane, finds 
--the intersection of the associated lines v0*x+w0*y=1 and v1*x+w1*y=1
twoIntersection = (v,w) ->
(
     if v#0*w#1-v#1*w#0 != 0 then 
     (
	  x := (w#1-w#0)/(v#0*w#1-v#1*w#0);
	  y := (v#0 - v#1)/(v#0*w#1-v#1*w#0);
	  P := {x,y};
     ) else P = {0,0};
P
)

--Given two vectors v={v0,...vn} and w={w0,...,wn}, list the 
--intersections of all lines vi*x+wi*y=1 and vj*x+wj*y=1
allIntersections = (v,w) ->
(
     L := new MutableList;
     c := 0;
     for i from 0 to #v-1 do
     (
	  for j from i+1 to #v-1 do 
     	  (
     	       if twoIntersection({v#i,v#j}, {w#i,w#j}) != {0,0} then 
     	       (
	  	    L#c = twoIntersection({v#i,v#j}, {w#i,w#j});
	  	    c = c+1;
     	       );
	  );
     );
     for i from 0 to #v-1 do
     (
	  if v#i !=0 then  
	  (
	       L#c = {1/(v#i), 0};
	       c = c + 1;
	  );
     );
     for i from 0 to #v-1 do
     (
	  if w#i !=0 then  
	  (
	       L#c = {0, 1/(w#i)};
	       c = c + 1;
	  );
     ); 
     K := new MutableList;
     c = 0; for i from 0 to #L-1 do
     (
	  if (L#i)#0 >= 0 and (L#i)#1 >=0 then (K#c = {(L#i)#0, (L#i)#1}; c = c+1);
     );
     K
)

--Given a point a=(x,y) in the real plane and two vectors v={v0,...,vn} and w={w0,...,wn}, checks whether a is in the polytope defined by the equations vi*x+wi*y<=1
isInPolytope = (a,v,w) ->
(
     alert := true;
     for i from 0 to #v-1 do
     (
	  if v#i*a#0 + w#i*a#1 > 1 then alert = false;
     );
     alert
)


--Given a point a=(x,y) in the real plane and two vectors
--v={v0,...,vn} and w={w0,...,wn}, checks whether a is in
--the polytope defined by the equations vi*x+wi*y<=1
isInInteriorPolytope = (a,v,w) ->
(
     alert := true;
     for i from 0 to #v-1 do
     (
	  if v#i*a#0 + w#i*a#1 >= 1 then alert = false;
     );
     alert
)

--Given two vectors v and w of the same length, outputs 
--a list of the defining vectors of the polytope as in isInPolytope
polytopeDefiningPoints = (v,w) ->
(
     L := allIntersections(v,w);
     K := new MutableList;
     c := 0;
     for j from 0 to #L-1 do
     (
	  if isInPolytope(L#j,v,w) == true then (K#c = {(L#j)#0, (L#j)#1}; c = c+1;)
     );
     K
)

--Given a list of coordinates in the real plane, 
--outputs the one with the largest coordinate sum
maxCoordinateSum = L ->
(
     K := new MutableList from {0,0};
     for i from 0 to #L-1 do if (L#i)#0 + (L#i)#1 > K#0 + K#1 then K = {(L#i)#0, (L#i)#1};
     K
)

--Finds the "delta" in Daniel Hernandez's algorithm
--for F-pure thresholds of binomials
dCalculation = (w,N,p) ->
(
     d := 0; for j from 0 to #w-1 do  d = d + digit(N+1,w#j,p);
     i := N; while d > p-2 do 
     (
	  d = 0; for j from 0 to #w-1 do  d = d + digit(i,w#j,p);
	  i = i - 1;
     );
     i + 1
)

--Given the "truncation" point (P1,P2) and two vectors 
--defining the binomial v and w, outputs the "epsilon" in 
--Daniel Hernandez's algorithm for F-thresholds of binomials
calculateEpsilon = (P1,P2,v,w) ->
(
     X := new MutableList;
     Y := new MutableList;
     c:=0; d := 0; for i from 0 to #v-1 do 
     (
	  if w#i != 0 then 
     	  (
	       X#c = (1 - (v#i)*(P2#0) - (w#i)*(P2#1))/(w#i);
	       c = c+1;
	  );
          if v#i != 0 then 
	  (
	       Y#d = (1 - (v#i)*(P1#0) - (w#i)*(P1#1))/(v#i);
	       d = d+1;
	  );
     );
     i:=0;
     epsilon:=0;
     if isInInteriorPolytope(P1,v,w)==false and isInInteriorPolytope(P2,v,w)==false then epsilon = -1 else
     (
	  if isInInteriorPolytope(P1,v,w)==false then for i from 0 to #v-1 do X#1 = 0;
	  if isInInteriorPolytope(P2,v,w)==false then for i from 0 to #v-1 do Y#1 = 0;
	  M := X#0; 
	  N := Y#0;
	  for i from 1 to #X-1 do M = min(M, X#i);
	  for j from 1 to #Y-1 do N = min(N, Y#j);
	  epsilon = max(M,N); 
     );
     epsilon
)

--Computes the FPT of a binomial, based on the work of Daniel Hernandez 
--(implemented by Emily Witt)
binomialFPT = g ->
(
     p := char ring g;
     v := (multiDegree(g))#0;
     w := (multiDegree(g))#1;
     FPT := 0;
     f := monomialFactor(v,w);
     x := factorOutMonomial(v,w);
     v = x#0;
     w = x#1;
     Q := maxCoordinateSum(polytopeDefiningPoints(v,w));
     if Q#0+Q#1 > 1 then FPT = 1 else
     (
	  L :=  firstCarry(Q,p);
	  if L == -1 then FPT = Q#0+Q#1 else
     	  (
     	       d := dCalculation(Q,L-1,p);
     	       P := (truncation(d,Q#0,p),  truncation(d,Q#1,p));
     	       P1 := {P#0, P#1+1/p^d};
     	       P2 := {P#0+1/p^d,P#1};
     	       FPT = truncation(L-1,Q#0+Q#1,p);
     	       if calculateEpsilon(P1,P2, v, w) != -1 then FPT = FPT +  calculateEpsilon(P1, P2, v, w);
     	  );
     );
     monFPT := infinity;
     for i from 0 to #f-1 do (if f#i!=0 then monFPT = min(monFPT, 1/(f#i)););
     if monFPT != 0 then FPT = min(FPT, monFPT);
     FPT
)

--Returns true if the polynomial is binomial.
isBinomial = f ->
(
     alert := true;
     if #(terms f)>2 then alert = false;
     alert
)

---------------------------------------------------------------------
--*****************************************************************--
--Functions for computing F-thresholds of forms in two variables   --
--over finite fields. Based on the work of Hernandez and Teixeira. -- 	                                           --                      
--*****************************************************************--
---------------------------------------------------------------------

{*
    Remark: At this point, only commands for computations of F-pure thresholds are
    implemented. Eventually computations of F-thresholds with respect to more general
    ideals will be implemented, and perhaps of more general polynomials. Some structures 
    and functions below are already designed to handle such greater generality. 
*}
    
{*
    Types and auxiliary commands
*}

--FTData is a HashTable that stores the data necessary in F-threshold calculations
--(for conveniently passing those data from one function to another).
--It contains the following keys:
--    "ring": the ring of the polynomial in question;
--    "char": the characteristic of ring;
--    "ideal": the ideal with respect to which we want to compute the F-threshold;
--    "gens": the generators of the ideal;
--    "polylist": a list of the (non-associated) factors of the polynomial in question;
--    "numpolys": the number of factors.
FTData = new Type of HashTable

--setFTData takes a list of generators of the ideal or the ideal itself and a list
--    of polynomials, and builds an FTData from them.
setFTData = method()

setFTData (List,List) := (gen,polylist) -> 
(
    	A:=ring gen_0;
    	p:= char A;	
	new FTData from {"char"=>p,"ring"=>A, "ideal"=>ideal gen, "gens" => gen,
	    "numpolys"=>#polylist,"polylist"=>polylist}
)

setFTData (Ideal,List) := (I,polylist) -> setFTData(I_*,polylist)

{*
    Tests and auxiliary functions
*}

--isInUpperRegion(a,q,S)/isInUpperRegion(u,S) test if the point u=a/q is in the
--upper region attached to S. Suppose I is the ideal of the FTData S under consideration 
--and L={L_1,...,L_n} is the "polylist". Then a point a/q (where a=(a_1,...,a_n) is a 
--nonnegative integer vector and q a power of "char") is in the "upper region" if 
--L_1^(a_1)...L_n^(a_n) is in I^[q]; otherwise it is in the lower region.
isInUpperRegion = method()

isInUpperRegion (List,ZZ,FTData) := (a,q,S) -> 
(
    frob:=ideal apply(S#"gens",f->f^q);
    F:=product(S#"polylist",a,(f,i)->fastExp(f,i));
    (F % frob) == 0
)

isInUpperRegion (List,FTData) := (u,S) ->
    isInUpperRegion append(getNumAndDenom(u),S)

--isInLoweRegion(a,q,S)/isInLoweRegion(u,S) test if the point u=a/q is in the
--lower region attached to S.
isInLowerRegion = method()

isInLowerRegion (List,ZZ,FTData) := (a,q,S) -> not isInUpperRegion(a,q,S)

isInLowerRegion (List,FTData) := (u,S) -> not isInUpperRegion(u,S)

--neighborInUpperRegion(a,q,S)/neighborInUpperRegion(u,S): auxiliary commands that, 
--given a point u=a/q in the upper region, try to find a "neighbor" of the form 
--(a-e_i)/q that also lies in the upper region. If the search is successful, they return
--the first such neighbor found; otherwise they return nothing.
neighborInUpperRegion = method()

neighborInUpperRegion (List,ZZ,FTData) := (a,q,S) ->
(
    if isInLowerRegion(a,q,S) then (error "Expected point in the upper region.");
    n := S#"numpolys";
    posEntries := positions(a,k->(k>0));
    found := false;
    i:=0;
    local candidate;
    local neighbor;
    while ((not found) and (i<#posEntries)) do 
    (
	candidate=a-canVector(posEntries_i,n);
	if isInUpperRegion(candidate,q,S) then (found=true; neighbor=candidate);
	i=i+1;
    );
    if (not found) then null else (neighbor,q)
)

neighborInUpperRegion (List,FTData) := (u,S) -> 
(
    nbr:=neighborInUpperRegion append(getNumAndDenom(u),S);
    if nbr===null then nbr else (nbr_0)/(nbr_1)
)

--isCP(a,q,S)/isCP(u,S) test if u=a/q is a critical point, that is, if u is in the
--upper region but each neighbor (a-e_i)/q (where a_i>0) is not.
isCP = method()

isCP (List,ZZ,FTData) := (a,q,S) -> 
(
    if isInLowerRegion(a,q,S) then return false;
    neighborInUpperRegion(a,q,S)===null
)

isCP (List,FTData) := (u,S) -> isCP append(getNumAndDenom(u),S)

--findCPBelow(u,S) takes a point u in the upper region attached to S and finds a 
--critical point <= u with the same denominator.
findCPBelow = method()

findCPBelow (List,FTData) := (pt,S) ->
(
    if isInLowerRegion(pt,S) then (error "The point must be in the upper region.");
    nbr:=neighborInUpperRegion(pt,S);
    if nbr===null then return pt else findCPBelow(nbr,S)
)

{*
    Computation of FPTs
*}

--FPT2VarHomogInternal({a1,...an},S): if S#"polylist={L1,...,Ln} is a list of linear
--forms, FPT2VarHomogInternal({a1,...an},S) finds the FPT of the polynomial
--F=L1^(a1)...Ln^(an)
FPT2VarHomogInternal = method(Options => {MaxExp => infinity, PrintCP => false, Nontrivial => false})

FPT2VarHomogInternal (List,FTData) := opt -> (a,S) ->
(
    deg:=taxicabNorm(a);
    pos:=positions(a,k->(k>=deg/2));
    if (pos!={}) then return(1/a_(pos_0)); 
       -- if some multiplicity a_i is "too big", return 1/a_i
    p:=S#"char";
    den:=denom(2/deg);
    local mult;
    if (opt.Nontrivial) then mult = infinity
    else
    ( 
    	if gcd(S#"char",den)==1 then mult = multOrder(p,den)
	else
	(
	    F:=product(S#"polylist",a,(f,i)->f^i);
	    if isFPTPoly(F,2/deg) then (return (2/deg))
	    else mult = infinity
	)
    );    
    rng:=S#"ring";
    polys:=S#"polylist";
    I:=S#"ideal";
    ideals:={I};
    e:=0;
    dgt:=0;
    u:=2*a/deg;
    while (I != ideal(1_rng) and e < (opt.MaxExp) and e < mult) do 
    (
	e=e+1;
	dgt=digit(e,u,p);
	I=frobeniusPower(I,1):product(polys,dgt,(f,k)->f^k);
	ideals=append(ideals,I)
    );
    if I!=ideal(1_rng) then 
    (
	if e == mult then (return (2/deg)) 
	else error "Reached MaxExp."
    );    
    e0:=e-1;
    S1:=setFTData(ideals_e0,polys);
    cp:=findCPBelow(dgt/p,S1); 
    	--if some coordinate of cp is 0, its magnification may not be a CP
    while product(cp)==0 do 
    (
	e0=e0-1;
        -- zoom out one step and look for CP again
    	S1=setFTData(ideals_e0,polys);
	cp=findCPBelow(cp/p+digit(e0+1,u,p)/p,S1) 
    );
    cp=cp/p^e0+truncation(e0,u,p); -- "zoom out"
    if opt.PrintCP then print(toString cp);
    max apply(cp,a,(c,k)->c/k)
)

-----------------------
FPT2VarHomog = method(Options => {MaxExp => infinity, PrintCP => false})

--FPT2VarHomog(RingElement)
--FPT(F) computes the F-pure threshold of a form F in two variables. 
--KNOWN ISSUE: if the splitting field of F is too big, factor will not work.
FPT2VarHomog (RingElement) :=  opt ->  F ->
(    
   if not isNonConstantBinaryForm(F) then (
	error "FPT2VarHomog expects a nonconstant homogeneous polynomial in 2 variables."
    );
    -- because factoring is the weakness of this algorithm, we try to avoid it
    -- by first checking if fpt=lct
    deg:=(degree F)_0;
    if isFPTPoly(F,2/deg) then return 2/deg;
    R:=ring F;
    vv:=R_*;
    kk:=splittingField(F);
    a:= symbol a;
    b:= symbol b;
    S:=kk[a,b];
    G:=sub(F,{(vv#0)=>a,(vv#1)=>b});
    (L,m):=toSequence transpose factorList(G);
    FPT2VarHomogInternal(m,setFTData(S_*,L),MaxExp=>(opt.MaxExp),PrintCP=>(opt.PrintCP),Nontrivial=>true)
)

--FPT2VarHomog(List,List)
--Given a list L={L_1,...,L_n} of linear forms in 2 variables and a list m={m_1,...,m_n}
--of multiplicities, FPT2VarHomog(L,m) returns the F-pure threshold of the polynomial 
--L_1^(m_1)*...*L_n^(m_n). 
FPT2VarHomog (List,List) :=  opt -> (L,m) -> 
    FPT2VarHomogInternal(m,setFTData(gens ring L_0,L),MaxExp=>(opt.MaxExp),PrintCP=>(opt.PrintCP))


{*
    Miscellaneous.
*}

-- Some commands for dealing with vectors --

--canVector(i,n) returns the i-th canonical basis vector in dimension n
--Warning: for convenience, this uses Macaulay2's convention of indexing lists starting 
--with 0; so, for example, {1,0,0,0} is canVector(0,4), not canVector(1,4).
canVector = method()

canVector (ZZ,ZZ) := (i,n) -> 
(
    if ((i<0) or (i>=n)) 
        then (error "canVector(i,n) expects integers i and n with 0<=i<n.");   
    apply(n,j->if i==j then 1 else 0)
)
 
-- getNumAndDenom(u) takes a rational vector u and returns a pair (a,q), where a 
--is an integer vector and q an integer such that u=a/q.
getNumAndDenom = method()

getNumAndDenom (List) := u -> 
(
    den := lcm apply( u, denom );
    a := apply( u, n-> lift( n*den, ZZ ) );
    ( a, den )        
)

--Computes the taxicab norm of a vector.
taxicabNorm = method()

taxicabNorm (List) := u -> sum( u, abs )

-- Factorization of polynomials and splitting fields --

--factorList(F) factors the RingElement F and returns a list of pairs of the form
--{factor,multiplicity}.
factorList = method()

factorList (RingElement) := F -> apply( toList( factor(F) ), toList )

--splittingField returns the splittingField of a polynomial over a finite field
splittingField = method()

splittingField (RingElement) := F -> 
(
    if not isPolynomialOverFiniteField(F) 
        then (error "splittingField expects a polynomial over a finite field");
    p:=char ring F;
    ord:=(coefficientRing(ring F))#order;
    factors:=first transpose factorList(F);
    deg:=lcm select(flatten apply(factors,degree),i->i>0);
    GF(p,deg*floorlog(p,ord))
)

-- Some tests

--isBinaryForm(F) checks if F is a homogeneous polynomial in two variables.
--WARNING: what we are really testing is if the *ring* of F is a polynomial ring in two 
--variables, and not whether F explicitly involves two variables. (For example, if F=x+y 
--is an element of QQ[x,y,z], this test will return "false"; if G=x is an element of 
--QQ[x,y], this test will return "true".)
isBinaryForm = method()

isBinaryForm (RingElement) := F ->
(
    R:=ring F;
    isPolynomialRing(R) and numgens(R)==2 and isHomogeneous(F)
)

--isNonconstantBinaryForm(F) checks if F is a nonconstant homogeneous polynomial in two 
--variables. See warning under "isBinaryForm".
isNonConstantBinaryForm = method()

isNonConstantBinaryForm (RingElement) := F -> (isBinaryForm(F) and (degree(F))_0>0)

--isLinearBinaryForm(F) checks if F is a linearform in two variables. See warning 
--under "isBinaryForm".
isLinearBinaryForm = method()

isLinearBinaryForm (RingElement) := F -> (isBinaryForm(F) and (degree(F))_0==1)

--isPolynomialOverFiniteField(F) checks if F is a polynomial over a finite field.
isPolynomialOverFiniteField = method()

isPolynomialOverFiniteField (RingElement) := F ->
(
    R:=ring F;
    kk:=coefficientRing(R);
    try kk#order then (isPolynomialRing(R) and isField(kk))
    	else false   
)    
