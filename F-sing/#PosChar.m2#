newPackage( "PosChar",
Version => "0.2a", 
Date => "May 30th, 2015", 
Authors => {
     {Name => "Erin Bela",
     Email=> "ebela@nd.edu"
     },
     {Name => "David J. Bruce",
     Email => "djbruce@math.wisc.edu",
     HomePage => "http://www.math.wisc.edu/~djbruce/"
     },
     {Name => "Daniel Hernandez",
     Email => "dhernan@math.utah.edu",
     HomePage => "http://math.utah.edu/~dhernan/"
     },
     {Name => "Zhibek Kadyrsizova",
     Email => "zhikadyr@umich.edu"
     },
     {Name => "Mordechai Katzman",
     Email=> "m.katzman@sheffield.ac.uk",
     HomePage=> "http://www.katzman.staff.shef.ac.uk/"
     },
     {Name => "Sara Malec",
     Email=> "smalec@gsu.edu"
     },
     {Name => "Karl Schwede",
     Email => "schwede@math.psu.edu",
     HomePage => "http://math.utah.edu/~schwede/"
     },
     {Name => "Pedro Teixeira",
     Email => "pteixeir@knox.edu",
     HomePage => "http://www.knox.edu/academics/faculty/teixeira-pedro.html"
     },
     {Name=> "Emily Witt",
     Email=> "ewitt@umn.edu",
     HomePage => "http://math.umn.edu/~ewitt/"
     }
},
Headline => "A package for calculations of singularities in positive characteristic", 
DebuggingMode => true, 
Reload => true 
)
export{
    "aPower",
    "ascendIdeal", -- EthRoots.m2
    "ascendIdealSafe", -- EthRoots.m2
    "ascendIdealSafeList", -- EthRoots.m2
    "AscentCount", -- EthRoots.m2
    "basePExp",
    "basePExpMaxE",
    "BinomialCheck", -- FThresholds.m2
    "binomialFPT", -- FThresholds.m2
    "canonicalIdeal",
    "canVector", -- FThresholds.m2
    "carryTest",
    "digit", 	 
    "denom",
    "DiagonalCheck", -- FThresholds.m2
    "diagonalFPT", -- FThresholds.m2
    "divideFraction",
    "estFPT", -- FThresholds.m2
    "ethRoot", -- EthRoots.m2
    "ethRootSafe", -- EthRoots.m2
    "ethRootSafeList", -- EthRoots.m2   
    "factorList", -- FThresholds.m2
    "fancyEthRoot", -- EthRoots.m2 
    "fastExp",  --frobeniousPowers.m2
    "findCPBelow", -- FThresholds.m2
    "findGeneratingMorphisms",     --MK
    "findHSLloci",                 --MK
    "findTestElementAmbient",
    "FinalCheck", -- FThresholds
    "findAllCompatibleIdeals", 	--- MK
    "findQGorGen",
    "finduOfIdeal",
    "firstCarry", 
    "FPTApproxList", -- FThresholds.m2    
    "FPT2VarHomog", -- FThresholds.m2    
    "FPT2VarHomogInternal", -- FThresholds.m2
    "fracPart",
    "frobenius",
    "frobeniusPower",  --frobeniousPowers.m2 
    "fSig",
    "FTApproxList",  -- FThresholds.m2
    "FTHatApproxList",  -- FThresholds.m2
    "FullMap",--specifies whether the full data should be returned
    "getNumAndDenom",
    "genFrobeniusPower",   --frobeniousPowers.m2 
    "guessFPT",  -- FThresholds.m2
    "HSL",
    "imageOfRelativeCanonical",
    "imageOfTrace", --doesn't work!
    "isBinomial",  -- FThresholds.m2
    "isCP",  -- FThresholds.m2
    "isDiagonal",  -- FThresholds.m2
    "isFJumpingNumberPoly",  -- FThresholds.m2
    "isFPTPoly",  -- FThresholds.m2
    "isFPure",
    "isFRegularPoly",
    "isFRegularQGor",
    "isInLowerRegion",  -- FThresholds.m2
    "isInUpperRegion",  -- FThresholds.m2
--    "isJToAInIToPe",
    "isSharplyFPurePoly",
    "isMapSplit",
    "MaxExp",  -- FThresholds.m2
    "maxIdeal",
    "minimalCompatible", -- EthRoots.m2 
---    "Mstar",			--- MK
    "multOrder",
    "MultiThread", --FThresholds
    "nonFInjectiveLocus",   --MK
    "Nontrivial", --SpecialFThresholds
    "nu", -- FThresholds.m2
    "nuAlt", -- FThresholds.m2
    "NuCheck", -- FThresholds.m2
    "nuHat", -- FThresholds.m2
    "nuHatList", -- FThresholds.m2
    "nuList", -- FThresholds.m2
    "nuListAlt", -- FThresholds.m2
    "nuListAlt1", -- FThresholds.m2
    "num",
    "Origin", --FThresholds
    "OutputRange", -- FThresholds
    "paraTestModule",
    "paraTestModuleAmbient",
    "PrintCP",  -- FThresholds.m2
    "setFTData",  -- FThresholds.m2
    "sigmaAOverPEMinus1Poly", 
    "sigmaQGorAmb", --needs optimization
    "sigmaAOverPEMinus1QGor",      --needs optimization
    "splittingField",  -- FThresholds.m2
    "tauPoly",
    "tauNonPrincipalAOverPEPoly",
    "tauAOverPEMinus1Poly",
    "tauGor",--needs optimization
    "tauGorAmb",--needs optimization
    "tauQGor",--needs optimization
    "tauQGorAmb",--needs optimization
    "taxicabNorm",  -- FThresholds.m2
    "truncation",
    "truncationBaseP"
}
--This file has "finished" functions from the Macaulay2 workshop at Wake 
--Forest in August 2012.  Sara Malec, Karl Schwede and Emily Witt contributed
--to it.  Some functions, are based on code written by Eric Canton and Moty
-- Katzman
--
--UPDATE January 2014 at Macaulay2 workshop at MSRI:  Daniel Hernandez, Moty 
--Katzman, Karl Schwede, Pedro Teixeira, Emily Witt added more functionality.
--
--UPDATE May 2015 at Macaulay2 workshop atBoise State:  Erin Bela, DJ Bruce,
-- Daniel Hernandez, Zhibek Kadyrsizova, and Emily Witt improved/fixed/added 
--functionality.

----------------------------------------------------------------
--************************************************************--
-- Miscellaneous auxiliary functions 	     	     	      --
--************************************************************--
----------------------------------------------------------------
   

-- maxIdeal returns the ideal generated by the variables of a polynomial ring
maxIdeal = method()

maxIdeal ( PolynomialRing ) := R -> monomialIdeal R_*




---------------------------------------------------------------
--***********************************************************--
--
--                                  
--***********************************************************--
---------------------------------------------------------------
 
-- This function computes the element in the ambient ring S of R=S/I such that
-- I^{[p^e]}:I = (f) + I^{[p^e]}
-- If there is no such unique element, the function returns zero

findQGorGen=method();
findQGorGen (Ring,ZZ) := (Rk,ek) -> (
     Sk := ambient Rk; -- the ambient ring
     Ik := ideal Rk; -- the defining ideal
     pp := char Sk; --the characteristic
     Ikpp := frobeniusPower(Ik,ek);
     
     J1 := trim (Ikpp : Ik); --compute the colon
     Tk := Sk/Ikpp; --determine the ideal in 
     
     J2 := trim sub(J1, Tk);
     
     Lk := first entries gens J2;
     
     nk := #Lk;
     val := 0_Sk;
     
     if (nk != 1) then (
	  error "findGorGen: this ring does not appear to be (Q-)Gorenstein, or
	   you might need to work on a smaller chart.  Or the index may not divide p^e-1
	   for the e you have selected.";
     )
     else (
	  val = lift(Lk#0, Sk);
     );    
     val 
)
findQGorGen(Ring) := (R2) -> ( findQGorGen(R2, 1) )




--Finds a test element of a ring R = k[x, y, ...]/I (or at least an ideal 
--containing a nonzero test element).  It views it as an element of the ambient ring
--of R.  It returns an ideal with some of these elements in it.
--One could make this faster by not computing the entire Jacobian / singular locus
--instead, if we just find one element of the Jacobian not in I, then that would also work
--and perhaps be substantially faster
findTestElementAmbient = Rk -> (
     --Sk := ambient Rk;
     -- Ik := ideal Sk;
     
     Jk := ideal singularLocus(Rk);
     if (isSubset(Jk, ideal Rk) == true) then 
          error "findTestElementAmbient: No test elements found, is the ring non-reduced?";
	  
     
     Jk          
)


--Outputs the test ideal of a (Q-)Gorenstein ring (with no pair or exponent)
--ek is the number such that the index divides (p^ek - 1)
--It actually spits out the appropriate stable/fixed ideal inside the ambient ring
tauQGorAmb = (Rk, ek) -> (
     Jk := findTestElementAmbient(Rk);
     hk := findQGorGen(Rk, ek);

     sub(ascendIdeal(Jk,hk,ek),Rk)
)

--Computes the test ideal of an ambient Gorenstein ring
tauGorAmb = (Rk) -> (tauQGorAmb(Rk, 1))

--Computes the test ideal of (R, f^(a/(p^e - 1)))
--when R is a polynomial ring.  This is based upon ideas of Moty Katzman.
tauAOverPEMinus1Poly = (fm, a1, e1) -> (
     Rm := ring fm;
     pp := char Rm;
     a2 := a1 % (pp^e1 - 1);
     k2 := a1 // (pp^e1 - 1); --it seems faster to use the fact that tau(f^(1+k)) = f*tau(f^k) 
     --this should be placed inside a try, and then if it fails we should be smarter...
     --fpow := fastExp(fm,a2);
     --IN := eR(ideal(fpow*fm),e1);  --the idea contained inside the test ideal.
     IN := ethRootSafe(fm, ideal(fm), a2, e1);
     
     IN = ascendIdealSafe(IN, fm, a2, e1);
     -- this is going to be the new value.  The *fm is a test element
     --5/0;
     --return the final ideal
     IN*ideal(fm^k2)
)

--Computes the test ideal of (R, f^t) when R 
--is a polynomial ring over a perfect field.
tauPolyOld = (fm, t1) -> (
     Rm := ring fm; 
     pp := char Rm;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     local I1;
     --first we compute tau(fm^{a/(p^c-1)})
     if (L1#2 != 0) then 
          I1 = tauAOverPEMinus1Poly(fm,L1#0,L1#2) else I1 = ideal(fm^(L1#0));     
	  
       
     
     --now we compute the test ideal using the fact that 
     --tau(fm^t)^{[1/p^a]} = tau(fm^(t/p^a))
     if (L1#1 != 0) then 
          ethRoot(I1, L1#1) else I1
)

--a slightly faster tauPoly
tauPoly = (fm, t1) -> (
     Rm := ring fm; 
     pp := char Rm;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     local I1;
     --first we compute tau(fm^{a/(p^c-1)})
     if (L1#2 != 0) then (
     	I1 = tauAOverPEMinus1Poly(fm,L1#0,L1#2);
     	if (L1#1 != 0) then
     		I1 = ethRoot(I1, L1#1)
     	)
     else (
     	if (L1#1 != 0) then
     		I1 = ethRootSafe(fm, ideal( sub(1, Rm)), L1#0, L1#1 )
     	else
 	    		I1 = ideal(fm^(L1#0))
 	    	);
     I1
)

--This is an internal function
--It is used to compute the test ideals of pairs (R, fm^(a1/p^e1-1)) where
--R = Sk/Ik.
--Inputs are Jk, a nonzero ideal contained in the test ideal
--hk, the multiple used to form phi of the ambient ring.  ek is the power associated with hk
--a1 and e1 and fm are as above
tauAOverPEMinus1QGorAmbOld = (Sk, Jk, hk, ek, fm, a1, e1) -> (
     pp := char Sk;
     et := lcm(ek, e1);
     hk1 := (hk)^(numerator ((pp^et - 1)/(pp^ek - 1)));  
       --hk for higher powers are simply appropriate powers of hk for lower powers, 
       --may as well take advantage of that
     a3 := numerator (a1*(pp^et - 1)/(pp^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
     
     a2 := a3 % (pp^et - 1);
     k2 := a3 // (pp^et - 1); --it seems faster to use the fact 
                              --that tau(f^(1+k)) = f*tau(f^k) 
     fpow := fm^a2; 
     
     Iasc := ascendIdeal(Jk*ideal(fm), fpow*hk1, et);
    
     Iasc*ideal(fm^k2)
)

tauAOverPEMinus1QGorAmb = (Sk, Jk, hk, ek, fm, a1, e1) -> (
     pp := char Sk;
     et := lcm(ek, e1);
     
     ak1 := numerator ((pp^et - 1)/(pp^ek - 1)); --an exponent for hk
     a3 := numerator (a1*(pp^et - 1)/(pp^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
                                               
	a2 := a3 % (pp^et - 1);
     k2 := a3 // (pp^et - 1); --it seems faster to use the fact that we can do simple Skoda for tau
     
     Jl := ascendIdealSafe(Jk, hk, 1, ek);
                      
        --          assert false;                             
     Iasc := ascendIdealSafeList(Jk*ideal(fm)^(ceiling(a3/(pp^et - 1))), (fm, hk), (a2, numerator ((pp^et - 1)/(pp^ek - 1))), et);
     
--     assert false;
     
     Iasc*ideal(fm^k2)
)


--Computes the test ideal of (Rk, fk^t1).  Here we assume the index of the canonical
--divides (p^ek - 1)
tauQGor = (Rk, ek, fk, t1) -> (
     Sk := ambient Rk;
     pp := char Sk;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     hk := findQGorGen(Rk, ek); --the term in the test ideal
     Jk := findTestElementAmbient(Rk); --this finds some test elements (lifted on the ambient
                                       --ring).  Right now it is slow because of a call to 
				       --singularLocus (ie, computing a Jacobian).
     I1 := ideal(0_Sk); I2 := ideal(0_Sk);
     fm := lift(fk, Sk); --we lift our f to the ambient polynomial ring
     a1 := L1#0; e1 := L1#2; pPow := L1#1; --t1 = a1 / (pp^pPow * (pp^e1 - 1))
     
     --before continuing, we need to rewrite a1 / (pp^pPow * (pp^e1 - 1)) as 
     --                                      a3 / (pp^(n1*ek) * (pp^e1 - 1))
     --the point is that ek needs to divide pPow
     remain := pPow % ek;
     dualRemain := ek - remain;
     
     pPow = pPow + dualRemain; --note in the end here, ek divides pPow
     a3 := a1*pp^(dualRemain);
     
     if (e1 != 0) then assert (t1 == a3/((pp^e1-1)*pp^pPow) ) else assert (t1 == a3/(pp^pPow) );
     
     d1 := pp^(pPow); if (e1 != 0) then d1 = d1*(pp^e1-1); --this is our denominator, used
                                                           --for simplifying computations
     a2 := a3 % d1;
     k2 := a3 // d1; --it seems faster to use the fact 
                              --that tau(f^(k+t)) = f^k*tau(f^t).  We'll toss on the multiple 
			      --f^k at the end
	     			  
     local I2;
     --first we compute tau(fk^{a2/(p^e1-1)})
     if (e1 != 0) then (
          I1 = tauAOverPEMinus1QGorAmb(Sk,Jk,hk,ek,fm,a2,e1);
          if (pPow != 0) then (
          	I2 = ethRootSafe(hk, I1, numerator((pp^pPow - 1)/(pp^ek - 1)), pPow)
		)
		else I2 = I1
     )
     else (
	  	I1 = ascendIdeal(Jk, hk, ek);
	  	if (pPow != 0) then (
	  		I2 = ethRootSafeList( (hk, fm), I1, (numerator((pp^pPow - 1)/(pp^ek - 1)), a2), pPow)
	  	)
	  	else I2 = I1
     );

	  
     sub(I2, Rk)*ideal(fk^k2)
)

--Computes tau(Rk,fk^tk), assuming Gorenstein rings
tauGor = (Rg,fg,tg) -> tauQGor (Rg,1,fg,tg)


----------------------------------------------------------------
--************************************************************--
--Test ideals for non-principal ideals                        --
--************************************************************--
----------------------------------------------------------------

flattenedReesAlgebra = (I1) -> (--takes an ideal, forms the rees algebra, and returns the rees algebra in two ways, first with flattened variables and the second without
	S1 := reesAlgebra I1;
	J1 := ideal S1;
	tempMonoid := S1.FlatMonoid;
	k1 := coefficientRing (ring I1);
	S2 := k1 tempMonoid;
	
	J2 := sub(J1, S2);
	
	(S2/J2, S1)
)

needsPackage "BGG"; --we'll be pushing forward...

needsPackage "Divisor";

tauNonPrincipalAOverPEPoly = {Verbose=> false}>> o -> (I1, a1, e1) -> ( -- computes \tau(I^{a/p^e}) for I an ideal in a polynomial ring
	if ( not(codim(I1) > 1)) then error "We can only handle ideals of codimension > 1 at this time.";
	
	--this function currently doesn't take advantage of Skoda's theorem, this will need to be done

	reesList := flattenedReesAlgebra I1;
	A1 := reesList#0; --this one has flattened variables
	A2 := reesList#1;
 	irrIdeal := sub(ideal(first entries vars A1), A1);
 	singLocus := ideal singularLocus (A1);
 	
 	IRees := sub(I1, A2);
 	
 	canList := canonicalIdeal(A1, FullMap=>true);
 	canIdeal := canList#0;
 	canMap := canList#1;
 	
 	paraTest := paraTestModuleAmbient(A1, canIdeal); 
 		
 	newMap := map(A1^1/(paraTest#0), canList#2, matrix(canMap));
 	newKer := (ker newMap)**A2; --this is the parameter test submodule of the canonical module  

	flag := false;
	i1 := e1;
	R1 := ring I1;
	p1 := char R1;
	ascend := I1; --dummy variables for checking whether we are done
	descend := ideal(sub(1, R1)); --dummy variables for checking whether we are done
	
	while (flag == false) do (
		ascend = fancyEthRoot(I1, a1*p1^(i1-e1), i1);
		if (o.Verbose == true) then (print  "Ascending ideal"; print ascend);
		
		flag = isSubset(descend, ascend);
		if (o.Verbose == true) then (print "flag"; print flag);
		if (flag == false) then (
			
			myDirectImage := HH_0(directImageComplex(IRees^(a1*p1^(i1-e1))*newKer, Regularity=>(10+a1))); 	
 	
		 	directIdeal := moduleToIdeal(myDirectImage, R1);
 			if ( codim(directIdeal)==1) then error "This function produced a codimension 1 ideal.";
 	
 			descend = ethRoot(directIdeal, i1);
 			if (o.Verbose == true) then (print  "Descending ideal"; print descend)
		);
		

		flag = isSubset(descend, ascend);
		
		--the following should be removed eventually, it is only here for debug purposes
		if ((flag == true) and (isSubset(ascend, descend)==false)) then error "Major error detected";
		i1 = i1+1;
		if (o.Verbose==true) then (print "Loop complete, continue?"; print (not flag) );
	);
	
	ascend
)

----------------------------------------------------------------
--************************************************************--
--Functions for computing sigma                               --
--************************************************************--
----------------------------------------------------------------


--Computes Non-Sharply-F-Pure ideals over polynomial rings for (R, fm^{a/(p^{e1}-1)}), 
--at least defined as in Fujino-Schwede-Takagi.
sigmaAOverPEMinus1Poly ={HSL=> false}>> o -> (fm, a1, e1) -> ( 
     Rm := ring fm;
     pp := char Rm;
     m1 := 0;
	e2 := e1;
	a2 := a1;
	--if e1 = 0, we treat (p^e-1) as 1.  
     if (e2 == 0) then (e2 = 1; a2 = a1*(pp-1));
     if (a2 > pp^e2-1) then (m1 = floor((a2-1)/(pp^e2-1)); a2=((a2-1)%(pp^e2-1)) + 1 );
     --fpow := fm^a2;
     IN := eR(ideal(1_Rm),e2); -- this is going to be the new value.
     -- the previous commands should use the fast power raising when Emily finishes it
     IP := ideal(0_Rm); -- this is going to be the old value.
     count := 0;
     
     --our initial value is something containing sigma.  This stops after finitely many steps.  
     while (IN != IP) do(
		IP = IN;
	  	IN = ethRootSafe(fm,IP,a2,e2); -- eR(ideal(fpow)*IP,e2);
	  	count = count + 1
     );

     --return the final ideal and the HSL number of this function
     if (o.HSL == true) then {IP*ideal(fm^m1),count} else IP*ideal(fm^m1)
)

--Computes Non-Sharply-F-pure ideals for non-polynomial rings with respect to no pair.
sigmaQGor ={HSL=> false}>> o -> (Rm, gg) -> (
     Sm := ambient Rm; --the polynomial ring that Rm is a quotient of
     hk := findQGorGen(Rm, gg);
     
     IN := ideal(1_Sm);
     count := 0;
     IP := ideal(0_Sm);
     
     while (IN != IP) do(
     	IP = IN;
     	IN = eR(ideal(hk)*IP,gg);
     	count = count + 1
     );
     
     --return the ideal and HSL
     if (o.HSL == true) then {sub(IP,Rm), count} else sub(IP, Rm)
)

--Computes Non-Sharply-F-Pure ideals for non-polynomial rings
--gg is the Gorenstein index
sigmaAOverPEMinus1QGor  ={HSL=> false}>> o -> (fk, a1, e1, gg) -> (
     Rm := ring fk;
     Sm := ambient Rm; --the polynomial ring that Rm is a quotient of
     pp := char Rm;
     ek := lcm(e1,gg); --we need to raise things to common powers
     hk := findQGorGen(Rm, gg); --it will be faster to find the Q-Gorenstein generator
     hl := hk^(sub((pp^ek - 1)/(pp^gg - 1), ZZ) ); --
	fm := lift(fk, Sm); --lift fk to the ambient ring
	fpow := fm^(a1*sub( (pp^ek - 1)/(pp^e1-1), ZZ) );


	IN := sigmaAOverPEMinus1Poly(hk,1,gg);
	count := 0;
	IP := ideal(0_Sm);

	while (IN != IP) do(
		IP = IN;
		IN = eR(ideal(fpow*hl)*IP, e1);
		count = count + 1
	);
	
     --return the final ideal
     if (o.HSL == true) then {sub(IP,Rm), count} else sub(IP,Rm)
	
)



----------------------------------------------------------------
--************************************************************--
--Functions for checking whether a ring/pair is F-pure/regular--
--************************************************************--
----------------------------------------------------------------

-- Given an ideal I of polynomial ring R
-- this uses Fedder's Criterion to check if R/I is F-pure
-- Recall that this involves checking if I^[p]:I is in m^[p]
-- Note:  We first check if I is generated by a regular sequence.

isFPure = I1->(
    maxIdeal:= monomialIdeal(first entries vars ring I1);
    local answer;
    local cond;
    p1:=char ring I1;
    if codim(I1)==numgens(I1) then(
	L:=flatten entries gens I1;
	cond = isSubset(ideal(product(#L, l-> fastExp(L#l, p1-1))),frobeniusPower(maxIdeal,1));
	if(cond==false) then answer=true else answer=false;
	)
    else(
	cond = isSubset((frobeniusPower(I1,1)):I1,frobeniusPower(maxIdeal,1));
	if(cond==false) then answer=true else answer=false;
	);
    answer
)

isFRegularPoly = method();

--Determines if a pair (R, f^t) is F-regular at a prime
--ideal Q in R, R is a polynomial ring  
isFRegularPoly (RingElement, QQ, Ideal) := (f1, t1, Q1) -> (
     not isSubset(tauPoly(f1,t1), Q1)
)

--Determines if a pair (R, f^t) is F-regular, R a polynomial ring
isFRegularPoly (RingElement, QQ) := (f1, t1) -> (
     isSubset(ideal(1_(ring f1)), tauPoly(f1,t1))
)

--Checks whether (R, f1^(a1/(p^e1-1)) is sharply F-pure at the prime ideal m1
isSharplyFPurePoly = (f1, a1, e1,m1) -> (
     if (isPrime m1 == false) then error "isSharplyFPurePoly: expected a prime ideal.";
     not (isSubset(ideal(f1^a1), frobeniusPower(m1,e1)))
)

--Checks whether a Q-Gorenstein pair is strongly F-regular 
isFRegularQGor = method();

--Computes whether (R, f1^t1) is F-regular, assuming the index of R divides p^e1-1
isFRegularQGor (ZZ, RingElement, QQ) := (e1,f1, t1) ->(
     R := ring f1;
     isSubset(ideal(1_R), tauQGor((ring f1),e1,f1,t1))
)

--Computes whether (R, f1^t1) is F-regular at Q1, assuming the index of R divides p^e1-1
isFRegularQGor (ZZ, RingElement, QQ, Ideal) := (e1,f1, t1, Q1) ->(
     not isSubset(tauQGor((ring f1),e1,f1,t1), Q1)
)

--Assuming no pair
isFRegularQGor (Ring,ZZ) := (R,e1) ->(
     isSubset(ideal(1_R), tauQGor(R,e1,1_R,1/1 ) )
)

--Assuming no pair checking at Q1
isFRegularQGor (Ring,ZZ,Ideal) := (R,e1,Q1) ->(
     not isSubset(tauQGor(R,e1,1_R,1/1 ),Q1 )
)


----------------------------------------------------------------
--************************************************************--
--Auxiliary functions for F-signature and Fpt computations.   --
--************************************************************--
----------------------------------------------------------------
 
--- Computes the F-signature for a specific value a1/p^e1
--- Input:
---	e1 - some positive integer
---	a1 - some positive integer between 0 and p^e
---	f1 - some polynomial in two or three variables in a ring R of PRIME characteristic
--- Output:
---	returns value of the F-signature of the pair (R, f1^{a1/p^e1})
--- Code is based on work of Eric Canton
fSig = (f1, a1, e1) -> (
     R1:=ring f1;
     pp:= char ring f1;     
     1-(1/pp^(dim(R1)*e1))*
          degree( (ideal(apply(first entries vars R1, i->i^(pp^e1)))+ideal(fastExp(f1,a1) ))) 
)  

--Calculates the x-int of the secant line between two guesses for the fpt
--Input:
--     t - some positive rational number
--     b - the f-signature of (R,f^{t/p^e})
--     e - some positive integer
--     t1- another rational number > t
--     f - some polynomial in two or three variables in a ring of PRIME characteristic
--
-- Output:
--	fSig applied to (f,t1,e)
--	x-intercept of the line passing through (t,b) and (t1,fSig(f,t1,e))
threshInt = (f,e,t,b,t1)-> (
     b1:=fSig(f,t1,e);
{b1,xInt(t,b,t1/(char ring f)^e,b1)}
)


--********************************************
--Some functions for the purpose of checking whether a map of rings is a splitting.  It also computes images of (field) trace.
--********************************************

needsPackage "PushForward"; 

--checks whether f1 : R1 -> S1 splits as a map of R1-modules
isMapSplit = (f1) -> (
	J1 := imageOfRelativeCanonical(f1);
	val := false;
	if (1 % J1 == 0) then val = true;
	
	val
)

--computes the image of Hom_R1(S1, R1) -> R1.
imageOfRelativeCanonical = (f1) -> (
	outList := pushFwd(f1);
--	myGenerators := first entries (outList#1);	
	target1 := (outList#2)(sub(1, target f1));
	
	h1 := map(outList#0, (source f1)^1, target1);
	
	d1 := Hom(h1, (source f1)^1);
	
	trim ideal( first entries gens image d1)
)

--computes the image of trace : S \to R if S is a finite R-module.
imageOfTrace = (f1) -> (
	print "Warning, this only works right now if S is a free module.  We should try to fix it...";
	outList := pushFwd(f1);
	myGenerators := first entries (outList#1);	
	i := 0;
	traceList := {};
	newMap := null;
	newMatrix := null;
	S1 := target f1;
	
	while (i < #myGenerators) do (
		newMap = map(S1^1, S1^1, {{myGenerators#i}});
		newMatrix = pushFwd(newMap, f1);
		traceList = append(traceList, trace newMatrix);
		i = i+1;
	);
	
	trim ideal traceList
)

--computes the relative e-iterated Frobenius over the base ring (the absolute Frobenius in the case 
frobenius = (R1, e1) -> (
	p1 := char R1;
	genList := first entries gens R1;
	fPowerList := apply(genList, z->z^p1 );
	map(R1, R1, fPowerList);
)


--****************************************************--
--*****************Documentation**********************--
--****************************************************--

beginDocumentation()

doc ///
   Key
      PosChar 
   Headline
      A package for calculations in positive characteristic 
   Description
      Text    
         This will do a lot of cool stuff someday. 
///

doc ///
     Key
     	aPower
     Headline
        Finds the largest power of p dividing x.
     Usage
     	 aPower(x,p)
     Inputs 
		x:ZZ
		p:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the largest exponent e such that p^e divides x.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- START: Transfered to EthRootDoc
doc ///
     Key
     	ascendIdeal
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdeal(J, h, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.
///

doc ///
     Key
     	ascendIdealSafe
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdealSafe(J, h, a, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	a:ZZ
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h^a.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.  It differs from ascendIdeal in that it minimizes the exponents that h is raised to, this can make it faster or slower depending on the circumstances.
///
--- END: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	basePExp 
     Headline
        Base p Expansion of an integer N
     Usage
     	  basePExp(N,p) 
     Inputs
         N:ZZ
	 p:ZZ
     Outputs
        :List
     Description
     	Text
	     Given an integer N and a prime p, outputs the digits of the base p expansion of N in base p.
///

doc ///
     Key
     	basePExpMaxE
     Headline
        Computes the base-p expansion of N from digits zero to e-1.
     Usage
     	 basePExpMaxE(N,p,e)
     Inputs
     	 N:ZZ
	 	 p:ZZ
	 	 e:ZZ
     Outputs
         :List
     Description
	Text
	     This computes the base p expansion of N, from digits 0 to e-1.  The digits are given in a list, and come with leading zeros.  If fewer than e digits are required, the list is padded with zeros.  If more digits are required, the final digit lists them.  Little endian is first.  For example, if p=5 and N = 16, the basePExpMaxE(16,5,4) will return {1,3,0,0} (1 one, 3 fives, 0 twentyfives, 0 onehundred twentyfives).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	binomialFPT
     Headline
        Computes the F-pure threshold of a binomial polynomial.
     Usage
     	 binomialFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a binomial in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	carryTest
     Headline
        Finds the number of digits we must check to see whether x and y add without carrying.
     Usage
     	 carryTest(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  This function finds the number of digit places we must check to see if x and y add without carrying.
///

doc ///
     Key
     	denom
     	(denom,ZZ)
     	(denom,QQ)
     Headline
        Returns the denominator of a rational number.
     Usage
     	 denom(x)
     	 denom(y)
     Inputs 
		x:QQ
		y:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the denominator of a rational number or integer (in the latter case it returns 1).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	diagonalFPT
     Headline
        Computes the F-pure threshold of a diagonal polynomial.
     Usage
     	 diagonalFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a diagonal hypersurface in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	digit
	(digit,ZZ,QQ,ZZ)
	(digit,ZZ,List,ZZ)
     Headline
        Gives the e-th digit of the base p expansion 
     Usage
     	 d=digit(e,x,p), D=digit(e,X,p)
     Inputs
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   consisting of rational numbers
     Outputs
        d:ZZ
	    which is the e-th digit of the non-terminating base p expansion of x
	D:List
	    which contains the e-th digits of the entries of the list X
     Description
	Text
	     Gives the e-th digit, to the right of the decimal point, of the non-terminating base p expansion of x in [0,1]; threads over lists of rational numbers. 
///

doc ///
     Key
     	divideFraction
     Headline
        Converts a rational number into something of the form (a/(p^b p^(c-1)).
     Usage
     	 divideFraction(t, p)
     Inputs 
		t:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Given a rational number t and prime p, this function finds a list of integers {a,b,c} such that t= (a/(p^b p^(c-1)).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 estFPT
     Headline
         Atempts to compute the F-pure threshold, where e is the max depth to search in.  
     Usage
     	  estFPT(f,e,finalCheck=>V,Verbose=>W)
     Inputs
     	 f:RingElement
         e:ZZ
	 V:Boolean
	 W:Boolean
     Outputs
        L:List
	Q:QQ
     Description
     	  Text 
	      This tries to find an exact value for the fpt.  If it can, it returns that value.  Otherwise it should return a range of possible values (eventually).  It first checks to see if the ring is binonmial or diagonal.  In either case it uses methods of D. Hernandez.  Next it tries to estimate the range of the FPT using nu's.  Finally, it tries to use this to deduce the actual FPT via taking advantage of convexity of the F-signature function and a secant line argument.  finalCheck is a Boolean with default value True that determines whether the last isFRegularPoly is run (it is possibly very slow).  If FinalCheck is false, then a last time consuming check won't be tried.  If it is true, it will be.  Verbose set to true displays verbose output.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- START: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doc ///
     Key
     	 ethRoot
     Headline
        Computes $I^{[1/p^e]}$ in a polynomial ring over a perfect field
     Usage
     	  ethRoot(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     In a polynomial ring k[x1, ..., xn], I^{[1/p^e]} is the smallest ideal J such that J^{[p^e]} = FrobeniusPower(J,e) \supseteq I.  This function computes it.
///

doc ///
     Key
     	ethRootSafe
     Headline
        Computes (f^a*I)^{[1/p^e]} in such a way that we don not blow exponent buffers.
     Usage
     	 ethRootSafe(f, I, a, e)
     Inputs
     	 f:RingElement
	 I:Ideal
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the 1/p^e-th root of (f^a*I).  It does it while trying to minimize the power that f gets raised to (in case a is a large number).  This can either be faster or slower than ethRoot.
///
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	fastExp 
     Headline
        Computes powers of elements in rings of characteristic p>0 quickly.
     Usage
     	  fastExp(f,N) 
     Inputs
     	 f:RingElement
         N:ZZ
     Outputs
        :RingElement
     Description
	Text
	     In prime characteristic p > 0, raising a sum (a+b) to a power is more quickly done by simply computing a^p and b^p and adding them.  The basic strategy is to break up the exponent into its base p expansion, and then use the exponent rules.  For example, (x+y)^(3*p^2 + 5*p+2) = ((x+y)^3)^(p^2)*((x+y)^5)^p*(x+y)^2.
///

doc ///
     Key
     	findQGorGen
     Headline
        If R = S/I where S is a polynomial ring, returns the ring element with I^{[p^e]} : I = (f) + I^{[p^e]}.
     Usage
     	 findQGorGen(R, e)
     Inputs
     	 R:Ring
     Outputs
         :RingElement
     Description
	Text
	     If R is Q-Gorenstein with index not divisible by p, then I^{[p^e]} : I = (f) + I^{[p^e]}.  For some e.  This function tries to find the f.  If the argument e is left out then e is assumed to be 1.
///

doc ///
     Key
     	firstCarry
     Headline
        Finds the first spot where (the eth digit of x) + (the eth digit of y) >= p.
     Usage
     	 firstCarry(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  Finds the first place where (the eth digit of x) + (the eth digit of y) >= p, in other words where the numbers add with carrying.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 FPTApproxList
	 (FPTApproxList,Ideal,ZZ)
	 (FPTApproxList,RingElement,ZZ)
     Headline
        Gives a list of nu_I(p^d)/p^d for d=1,...,e.
     Usage
     	  FPTApproxList(I,e)
	  FPTApproxList(f,e) 
     Inputs
     	 I:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I(p^d)/p^d for d = 1, ..., e.  The {nu_I(p^d)/p^d} converge to the F-pure threshold.	     
///

doc ///
     Key
     	FPT2VarHomog
	(FPT2VarHomog,RingElement)
	(FPT2VarHomog,List,List)
     Headline
        F-pure threshold of a form in two variables
     Usage
     	  fpt=FPT2VarHomog(G), fpt=FPT2VarHomog(factors,multiplicities)
     Inputs 
	factors:List
	    which contains the linear factors of a form G in two variables 
	multiplicities:List
	    which contains the multiplicities of those linear factors in G
	G:RingElement
	    a form in two variables
     Outputs
        fpt:QQ
     Description
	Text
	    FPT2VarHomog computes the F-pure threshold of a homogeneous polynomial G
	    	in two variables. 
	    The polynomial G can be entered directly, or if the user knows a factorization
	    	G=L1^(a1)...Ln^(an) into linear forms, that can be used for improved 
		performance: FPT2VarHomog({L1,...,Ln},{a1,...,an}).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 frobeniusPower
     Headline
        The following raises an ideal to the $p^e$th power.
     Usage
     	  frobeniusPower(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     If I = ideal(x1, ..., xn), then frobeniusPower(I,e) outputs ideal(x1^(p^e), ..., xn^(p^e)) where p is the characteristic of the ring.
///

doc ///

     Key
     	 fSig
     Headline
        Computes the F-signature for a specific value $a/p^e$.
     Usage
     	  fSig(f,a,e)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
     Outputs
        :QQ
     Description
	Text
	     This computes the F-signature $s(R, f^{a/p^e})$ if R is a polynomial ring over a perfect field.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 FTApproxList
	 (FTApproxList,Ideal,Ideal, ZZ)
	 (FTApproxList,RingElement,Ideal,ZZ)
     Headline
        Gives a list of nu_I^J(p^d)/p^d for d=1,...,e.
     Usage
     	  FTApproxList(I,J,e)
	  FTApproxList(f,J,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I^J(p^d)/p^d for d = 1, ..., e.  The {nu_I^J(p^d)/p^d} converge to the F-threshold.	     
///

doc ///
     Key
     	 FTHatApproxList
	 (FTHatApproxList,Ideal,Ideal, ZZ)
	 (FTHatApproxList,RingElement,Ideal,ZZ)
     Headline
        Gives a list of nuHat_I^J(p^d)/p^d for d=1,...,e.
     Usage
     	  FTHatApproxList(I,J,e)
	  FTHatApproxList(f,J,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nuHat_I^J(p^d)/p^d for d = 1, ..., e.  The {nuHat_I^J(p^d)/p^d} converge to the FHat-threshold.	     
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	genFrobeniusPower 
     Headline
        Computes the generalized Frobenius power of an ideal
     Usage
     	  frobeniusPower(I1,e1)  
     Inputs
         	I1:Ideal
	 	e1:ZZ
     Outputs
        :Ideal
     Description
     	Text
	     Computes I^[N] for an ideal I and an integer N, where I^[N] is defined as follows. If N's base P-expansion is N=n_0+n_1P+...+n_eP^e then I^[N]=I^(n_0)*(I^(n_1))^[P]*...*(I^(n_e))^[P^e]. When P is prime I^[P^e] is the usual Frobenius power.
 ///
 
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	guessFPT 
     Headline
        Tries to guess the FPT in a really naive way (this should be improved).
     Usage
     	  guessFPT(f,e,d) 
     Inputs
     	 f:RingElement
         e:ZZ
	 d:ZZ
     Outputs
        :List
     Description
	Text
	     This tries to guess the FPT.  In particular, it computes the number nu such that nu/(p^e - 1) <= FPT < (nu+1)/p^e.  It then outputs a list of all rational numbers with denominators less than or equal to d, which lie in that range.  WARNING:  There are several improvements which should be made to this function to rule out many of the possibilies.
///

doc ///
     Key
     	isBinomial 
     Headline
        Checks whether a polynomial is binomial.
     Usage
     	 isBinomial(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a binomial, otherwise returns false.
///

doc ///
     Key
     	isDiagonal 
     Headline
        Checks whether a polynomial is diagonal.
     Usage
     	 isDiagonal(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a diagonal, otherwise returns false.  Recall f is called diagonal if it is of the form x_1^(a_1)+...+x_n^(a_n) up to renumbering of the variables.
///

doc ///
     Key
     	isFJumpingNumberPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFJumpingNumberPoly(f,t,Verbose=>V)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is an F-jumping number, otherwise it returns false.
///

doc ///
     Key
     	isFPTPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFPTPoly(f,t,Verbose=>V,Origin=>W)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
		W:Origin
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is the FPT, otherwise it returns false.  If Origin is true, it only checks it at ideal(vars ring f).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 isFPure 
     Headline
         Tests for a given ideal I, if R/I is F-pure. 
     Usage
     	 isFPure(I)
     Inputs
     	 I:Ideal
     Outputs
         :Boolean
     Description
	Text 
	    In the case where I is a complete intersection,this function applies Fedder's Criterion.
	    Otherwise, checks if I^[p]:I is contained in m^[p]. 
 	   
///

doc ///
     Key
     	 isFRegularPoly
     Headline
        Determines if a pair $(R, f^t)$ is F-regular when R is a polynomial ring. 
     Usage
     	  isFRegularPoly
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Boolean
     Description
	Text
	     This computes the test ideal.  The ring is F-regular if the test ideal is the whole ring, in which case this function returns true.  Otherwise, this function returns false.

///

doc ///
     Key
     	isFRegularQGor
	 (isFRegularQGor, ZZ, RingElement, QQ)
	 (isFRegularQGor, ZZ, RingElement, QQ, Ideal)
	 (isFRegularQGor, Ring, ZZ)
	 (isFRegularQGor, Ring, ZZ, Ideal)
     Headline
        Checks whether a ring or a pair is Q-Gorenstein.
     Usage
     	 isFRegularQGor(e,f,t)
     	 isFRegularQGor(e,f,t,Q)
     	 isFRegularQGor(R,e)
     	 isFRegularQGor(R,e,Q)
     Inputs 
		R:Ring
	 f:RingElement
	 e:ZZ
	 t:QQ
	 Q:Ideal
     Outputs
         :Boolean
     Description
	Text
	     Checks whether R, or the pair (R, f^t),  is strongly F-regular at Q (respectively the origin).  It assumes the Q-Gorenstein index divides (p^e - 1).
///

doc ///
     Key
     	 isSharplyFPurePoly
     Headline
        Checks whether (R, f^(a/(p^e - 1))) is F-pure at the prime ideal m.
     Usage
     	 isSharplyFPurePoly(f,a,e,m)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
	 m:Ideal
     Outputs
         :Boolean
     Description
	Text
	     This checks whether (R, f^(a/(p^e-1))) is F-pure at the prime ideal m at least in the case that R is a polynomial ring.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 nu
	 (nu,Ideal,Ideal,ZZ)
	 (nu,Ideal,ZZ)
	 (nu, RingElement,Ideal, ZZ)
	 (nu, RingElement, ZZ)
     Headline
        Gives $\(nu_I)^J(p^e)$ or $\(nu_f)^J(p^e)$
     Usage
     	  nu(I,J,e)
	  nu(I,e)
	  nu(f,J,e)
	  nu(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :ZZ
     Description
	Text
	    Given an ideal I in a polynomial ring k[x1, ..., xn], this function outputs the maximal integer nu such that I^nu is not in ideal J^[p^e].  If the input is (Ideal,ZZ) then the function computes the maximal integer nu such that I^nu in not in (x_1, ...,x_n)^[p^e]. If a RingElement is passed, it computes nu of the principal ideal generated by this element.This is used frequently to compute the F-pure threshold.
///

doc ///
     Key
     	 nuList
	 (nuList, Ideal,Ideal,ZZ)
	 (nuList, Ideal, ZZ)
	 (nuList, RingElement, Ideal, ZZ)
	 (nuList, RingElement, ZZ)
     Headline
        Lists $\(nu_I)^J(p^d)$ for d = 1,...,e.
     Usage
     	  nuList(I,J,e)
	  nuList(I,e)
	  nuList(f,J,e)
	  nuList(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :List
     Description
	Text
	     Given an ideal I in a polynomial ring k[x1,...,xn], this function computes nu(I,d) for d = 1,...,e. If a RingElement is passed, it computes nu of the principal ideal generated by this element for d=1,...,e
///

doc ///
     Key
     	 nuHat
	 (nuHat,Ideal,Ideal,ZZ)
	 (nuHat,Ideal,ZZ)
	 (nuHat, RingElement,Ideal, ZZ)
	 (nuHat, RingElement, ZZ)
     Headline
        Gives $\hat(nu_I)^J(p^e)$ or $\hat(nu_f)^J(p^e)$
     Usage
     	  nuHat(I,J,e)
	  nuHat(I,e)
	  nuHat(f,J,e)
	  nuHat(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :ZZ
     Description
	Text
	    Given an ideal I in a polynomial ring k[x1, ..., xn], this function outputs the maximal integer nu such that I^[nu] is not in ideal J^[p^e].  If the input is (Ideal,ZZ) then the function computes the maximal integer nu such that I^[nu] in not in (x_1, ...,x_n)^[p^e]. If a RingElement is passed, it computes nuHat of the principal ideal generated by this element.This is used frequently to compute the generalized F-pure threshold.
///

doc ///
     Key
     	 nuHatList
	 (nuHatList, Ideal,Ideal,ZZ)
	 (nuHatList, Ideal, ZZ)
	 (nuHatList, RingElement, Ideal, ZZ)
	 (nuHatList, RingElement, ZZ)
     Headline
        Lists $\hat(nu_I)^J(p^d)$ for d = 1,...,e.
     Usage
     	  nuHatList(I,J,e)
	  nuHatList(I,e)
	  nuHatList(f,J,e)
	  nuHatList(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :List
     Description
	Text
	     Given an ideal I in a polynomial ring k[x1,...,xn], this function computes nuHat(I,d) for d = 1,...,e. If a RingElement is passed, it computes nuHat of the principal ideal generated by this element for d=1,...,e
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	sigmaAOverPEMinus1Poly
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}) when R is a polynomial ring.
     Usage
     	 sigmaAOverPEMinus1Poly (f, a, e, HSL=>W)
     Inputs 
		f:RingElement
	a:ZZ
	e:ZZ
	W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by f^a.  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	sigmaAOverPEMinus1QGor
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}).
     Usage
     	 sigmaAOverPEMinus1QGor(f, a, e, g,HSL=>W)
     Inputs 
		f:RingElement
		a:ZZ
		e:ZZ
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace of R by f^a (we assume that the Q-Gorenstein index of R divides p^g-1).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number of sigma(ring f) relative to f.
///

doc ///
     Key
     	sigmaQGor
     Headline
        Computes the non-sharply F-pure ideal of R, where R is Q-Gorenstein with index dividing (p^g-1).
     Usage
     	 sigmaQGor(R, g,HSL=>W)
     Inputs 
		R:Ring
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the  g-th Frobenius trace of R (we assume that g is the Q-Gorenstein index of R).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	tauAOverPEMinus1Poly
     Headline
        Computes the test ideal of f^(a/(p^e-1)) if f is in a polynomial ring.
     Usage
     	 tauAOverPEMinus1Poly(f, a, e)
     Inputs
     	 f:RingElement
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the test ideal tau(f^(a/(p^e-1)) ).  The basic idea first appeared in a paper of Mordechai Katzman.
///

doc ///
     Key
     	tauGor
     Headline
        Computes tau(R,f^t) for a Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauGor(R,f,t)
     Inputs
     	 R:Ring
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Gorenstein ring.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing the images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	tauGorAmb
     Headline
        Computes tau(R) for a Gorenstein ring.
     Usage
     	 tauGorAmb(R)
     Inputs
     	 R:Ring
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a quasi-Gorenstein ring R.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S. 
///

doc ///
     Key
     	 tauPoly
     Headline
        Computes the test ideal of $(R, f^t)$.
     Usage
     	  tauPoly(f,t) 
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Ideal
     Description
	Text
	     This computes the test ideal of (R, f^t) when R is a polynomial ring over a perfect field.  It is done as follows.  If t = a/(p^e - 1) then tau(R, f^t) is computed as a sum of (f^{lceil t rceil}*f^{lceil t(p^e-1) rceil})^{[1/p^e]} until the sum stabilizes.  For the more general case, we use the formula tau(R, f^t)^{[1/p^d]} = tau(R, f^{t/p^d}).
///

doc ///
     Key
     	tauQGorAmb
     Headline
        Computes tau(R) for a Q-Gorenstein ring with index not dividing p^e - 1.
     Usage
     	 tauQGorAmb(R,e)
     Inputs
     	 R:Ring
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a Q-Gorenstein ring R with index dividing p^e - 1.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S.  e is the index of the canonical divisor on R.
///

doc ///
     Key
     	tauQGor
     Headline
        Computes tau(R,f^t) for a Q-Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauQGor(R,e,f,t)
     Inputs
     	 R:Ring
	 e:ZZ
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Q-Gorenstein ring such that the index divides p^e -1.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	truncation
	(truncation,ZZ,QQ,ZZ)
	(truncation,ZZ,List,ZZ)
     Headline
        Truncations of base p expansions of rational numbers
     Usage
     	 t=truncation(e,x,p), T=truncation(e,X,p)
     Inputs 
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   which contains rational numbers
     Outputs
        t:QQ
	    which is the e-th truncation of the non-terminating base p expansion of x
	T:List
	    which contains the e-th truncations of the entries of the list X
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of a nonnegative rational number x, as a fraction; threads over lists of rational numbers.
///

doc ///
     Key
     	truncationBaseP
     Headline
        Gives the first e digits of the non-terminating base p expansion of x.
     Usage
     	 truncationBaseP(e,x,p)
     Inputs 
		e:ZZ
	x:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of x in [0,1], as a list.
///

end
