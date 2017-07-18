--****************************************************
--****************************************************
--This file contains functions related to test ideals.
--****************************************************
--****************************************************
 
-- This function computes the element f in the ambient ring S of R=S/I such that
-- I^{[p^e]}:I = (f) + I^{[p^e]}.
-- If there is no such unique element, the function returns an error.
if  (not (class Fsing === Package)) then (
    needs "BasicFunctions.m2";
    needs "EthRoots.m2";
    needs "frobeniusPowers.m2";
    needs "parameterTestIdeal.m2";
);

needsPackage "Divisor";

--*********************
--Preliminary functions
--*********************

findQGorGen=method()

findQGorGen ( ZZ, Ring ) := ( e, R ) -> 
(
     S := ambient R; -- the ambient ring
     I := ideal R; -- the defining ideal
     Ie := frobenius( e, I );     
     J := trim ( Ie : I ); --compute the colon
     J = trim sub( J, S/Ie ); -- extend colon ideal to S/Ie
     L := J_*; -- grab generators
     if ( #L != 1 ) then 
	  error "findQGorGen: this ring does not appear to be (Q-)Gorenstein, or you might need to work on a smaller chart. Or the index may not divide p^e-1 for the e you have selected.";
     lift( L#0, S )
)

findQGorGen ( Ring ) := R -> findQGorGen( 1, R )

--Finds a test element of a ring R = k[x, y, ...]/I (or at least an ideal 
--containing a nonzero test element).  It views it as an element of the ambient ring
--of R.  It returns an ideal with some of these elements in it.
--One could make this faster by not computing the entire Jacobian / singular locus
--instead, if we just find one element of the Jacobian not in I, then that would also work
--and perhaps be substantially faster
--it assumes that R is a reduced ring.
findTestElementAmbient = method()
randomSubset = method()

findTestElementAmbient(Ring) := (R1) ->
(
	I1 := ideal R1;
	n1 := #gens R1 - dim R1;
	M1 := jacobian I1;
	r1 := rank target M1;
	c1 := rank source M1;
	testEle := sub(0,ambient R1);
	primesList := minimalPrimes I1;
    curMinor := ideal(sub(0, ambient R1));
	while(any(primesList, II->isSubset(ideal(testEle), II))) do(
	    curMinor = first entries gens  minors(n1,M1, First =>{randomSubset(r1,n1),randomSubset(c1,n1)}, Limit =>1);
	    if (#(curMinor) > 0) then (
            testEle = testEle + (random(coefficientRing R1))*(first curMinor);
        );
	);
	testEle
);

randomSubset(ZZ,ZZ) := (m,n) ->
(
	--L = for i from 0 to m-1 list i;
	L := toList(0..m-1);
	for i from 0 to m-n-1 do (L = delete(L#(random(0,m-1-i)),L));
	L
)

--****************************
--****************************
--**New test ideal functions**
--****************************
--****************************

--gets a nonzero generator of an ideal.
getNonzeroGenerator := (I2) -> (
    i := 0;
    flag := false;
    genList := first entries gens I2;
    localZero := sub(0, ring I2);
    while ((i < #genList) and (flag == false)) do (
        if (genList#i != localZero) then (            
            flag = true;
        );
        i = i + 1;
    );
    if (flag == true) then (
        genList#(i-1)
    )
    else (
        null
    )
);

--the following function should go elsewhere, it checks whether a given ideal is locally principal (really, invertible).  If it is locally principal, it returns the inverse ideal.
isLocallyPrincipalIdeal := (I2) -> (
    localGen := getNonzeroGenerator(I2);
    if (localGen === null) then (
        return {false, sub(0, ring I2)};
    );
    inverseIdeal := (ideal(localGen) : I2);
    idealProduct := inverseIdeal*I2;
    isLocPrinc := (reflexify(idealProduct) == idealProduct);
    if (isLocPrinc == true) then (
        return {true, inverseIdeal};
    )
    else (
        return {false, sub(0, ring I2)};
    );
    
);

--the following is the new function for computing test ideals written by Karl.

testIdeal = method(Options => {MaxCartierIndex => 100});

testIdeal(Ring) := o->(R1) -> (
    canIdeal := canonicalIdeal(R1);
    pp := char R1;
    cartIndex := 0;
    fflag := false;
    computedFlag := false;
    curIdeal := ideal(sub(0, R1));
    locPrincList := null;
    computedTau := ideal(sub(0, R1));
    while ( (fflag == false) and (cartIndex < o.MaxCartierIndex) ) do (
        cartIndex = cartIndex + 1;
        curIdeal = reflexivePower(cartIndex, canIdeal);
        locPrincList = isLocallyPrincipalIdeal(curIdeal);
        if (locPrincList#0 == true) then (            
            fflag = true;
        );
    );
    if (fflag == false) then error "testIdeal: Ring does not appear to be Q-Gorenstein, perhaps increase the option MaxCartierIndex";
    if ((pp-1)%cartIndex == 0) then ( 
        J1 := findTestElementAmbient( R1 );
        h1 := sub(0, ambient R1);
        try (h1 = findQGorGen( 1, R1)) then (
            computedTau = ascendIdeal(1, h1, sub(ideal(J1), R1) );
            computedFlag = true;
        )
        else (
            computedFlag = false;
        );
    );
    if (computedFlag == false) then (
        gg := first first entries gens trim canIdeal;
        dualCanIdeal := (ideal(gg) : canIdeal);
        nMinusKX := reflexivePower(cartIndex, dualCanIdeal);
        gensList := first entries gens trim nMinusKX;
        
        runningIdeal := ideal(sub(0, R1));
        omegaAmb := sub(canIdeal, ambient R1) + ideal(R1);
    	u1 := (findusOfIdeal(ideal R1, omegaAmb));
    
--    print gensList;
--    1/0;
        for x in gensList do (
            runningIdeal = runningIdeal + (testModule(1/cartIndex, sub(x, R1), canIdeal, u1))#0;        
        );
--    1/0;
    
        newDenom := reflexify(canIdeal*dualCanIdeal);
        computedTau = (runningIdeal*R1) : newDenom;
    );
    computedTau
)

testIdeal(QQ, RingElement, Ring) := o->(t1, f1, R1) -> (
    --this computes \tau(R, f^t)
    testIdeal({t1}, {f1}, R1)
);

testIdeal(ZZ, RingElement, Ring) := o->(t1, f1, R1) -> (
    --this computes \tau(R, f^t)
    testIdeal({t1/1}, {f1}, R1)
);

testIdeal(QQ, RingElement) := o->(t1, f1) -> (
    testIdeal({t1}, {f1}, ring f1)
);

testIdeal(ZZ, RingElement) := o->(t1, f1) -> (
    testIdeal({t1/1}, {f1}, ring f1)
);

testIdeal(List, List) := o->(tList, fList) ->(
    testIdeal(tList, fList, ring (fList#0))
);

testIdeal(List, List, Ring) := o->(tList, fList, R1) ->(
    canIdeal := canonicalIdeal(R1);
    pp := char R1;
    cartIndex := 0;
    fflag := false;
    computedFlag := false;
    curIdeal := ideal(sub(0, R1));
    locPrincList := null;
    computedTau := ideal(sub(0, R1));
    while ( (fflag == false) and (cartIndex < o.MaxCartierIndex) ) do (
        cartIndex = cartIndex + 1;
        curIdeal = reflexivePower(cartIndex, canIdeal);
        locPrincList = isLocallyPrincipalIdeal(curIdeal);
        if (locPrincList#0 == true) then (            
            fflag = true;
        );
    );
    if (fflag == false) then error "testIdeal: Ring does not appear to be Q-Gorenstein, perhaps increase the option MaxCartierIndex";
    if ((pp-1)%cartIndex == 0) then ( 
        J1 := findTestElementAmbient( R1 );
        h1 := sub(0, ambient R1);
        try (h1 = findQGorGen( 1, R1)) then (
            --do stuff
            computedTau = testModule(tList, fList, ideal(sub(1, R1)), {h1});
        ) else (
            computedFlag = false;
        );
    );
    if (computedFlag == false) then (
        gg := first first entries gens trim canIdeal;
        dualCanIdeal := (ideal(gg) : canIdeal);
        nMinusKX := reflexivePower(cartIndex, dualCanIdeal);
        gensList := first entries gens trim nMinusKX;
        
        runningIdeal := ideal(sub(0, R1));
        omegaAmb := sub(canIdeal, ambient R1) + ideal(R1);
    	u1 := (findusOfIdeal(ideal R1, omegaAmb));
    	
        t2 := append(tList, 1/cartIndex);
        f2 := fList;
    
        for x in gensList do (
            f2 = append(fList, x);
            runningIdeal = runningIdeal + (testModule(t2, f2, canIdeal, u1))#0;        
        );
    
        newDenom := reflexify(canIdeal*dualCanIdeal);
        computedTau = (runningIdeal*R1) : newDenom;
    );
    computedTau        
);

--We can now check F-regularity

isFregular = method(Options => {MaxCartierIndex => 100});

isFregular(Ring) := R1 -> (
    tau := testIdeal(R1);
    if (isSubset(ideal(sub(1, R1)), tau)) then true else false
);

isFregular(QQ, RingElement) := (tt, ff) -> (
    R1 := ring ff;
    tau := testIdeal(tt, ff);
    if (isSubset(ideal(sub(1, R1)), tau)) then true else false
);

isFregular(ZZ, RingElement) := (tt, ff) -> (
    R1 := ring ff;
    tau := testIdeal(tt, ff);
    if (isSubset(ideal(sub(1, R1)), tau)) then true else false
);

isFregular(List, List) := (ttList, ffList) -> (
    R1 := ring (ffList#0);
    tau := testIdeal(ttList, ffList);
    if (isSubset(ideal(sub(1, R1)), tau)) then true else false
);

--*************
--all the below versions should eventually be deprecated
--*************

--Outputs the test ideal of a (Q-)Gorenstein ring (with no pair or exponent)
--e is the number such that the index divides (p^e - 1)
--It actually spits ourt the appropriate stable/fixed ideal inside the ambient ring

--Computes the test ideal of (R, f^(a/(p^e - 1)))
--when R is a polynomial ring.  This is based upon ideas of Moty Katzman.
tauAOverPEMinus1Poly = ( f, a, e ) -> 
(
     R := ring f;
     p := char R;
     b := a % (p^e - 1);
     k := a // (p^e - 1); --it seems faster to use the fact that tau(f^(1+k)) = f*tau(f^k) 
     I := ethRoot( e, a, f, ideal(f) );     
     I = ascendIdeal( e, b, f, I );
     I*ideal(f^k)
)

--Computes the test ideal of (R, f^t), where t is rational and R is a polynomial ring. 
tauPoly = ( f, t ) -> 
(
     R := ring f; 
     p := char R;
     (a,b,c) := toSequence( divideFraction( p, t) ); --this breaks up t into the pieces we need
     local I;
     --first we compute tau(f^{a/(p^c-1)})
     if (c != 0) then 
     (
     	I = tauAOverPEMinus1Poly( f, a, c);
     	if (b != 0) then I = ethRoot( b, I)
     )
     else 
     (
     	if (b != 0) then I = ethRoot( b, a, f )
     	else I = ideal( f^a )
     );
     I
)

--This is an internal function
--It is used to compute the test ideals of pairs (R, fm^(a1/p^e1-1)) where
--R = Sk/Ik.
--Inputs are Jk, a nonzero ideal contained in the test ideal
--hk, the multiple used to form phi of the ambient ring.  ek is the power associated with hk
--a1 and e1 and fm are as above

tauAOverPEMinus1QGorAmb = (Sk, Jk, hk, ek, fm, a1, e1) -> 
(
     pp := char Sk;
     et := lcm(ek, e1);
     
     ak1 := numerator ((pp^et - 1)/(pp^ek - 1)); --an exponent for hk
     a3 := numerator (a1*(pp^et - 1)/(pp^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
                                               
	a2 := a3 % (pp^et - 1);
     k2 := a3 // (pp^et - 1); --it seems faster to use the fact that we can do simple Skoda for tau
     
     Jl := ascendIdeal(ek, 1, hk, Jk );
     expList := (a2, numerator ((pp^et - 1)/(pp^ek - 1)));
        --          assert false;                             
     Iasc := ascendIdeal( et, expList, (fm, hk), Jk*ideal(fm)^(ceiling(a3/(pp^et - 1))) );
--     assert false;
     
     Iasc*ideal(fm^k2)
)


--Computes the test ideal of (Rk, fk^t1).  Here we assume the index of the canonical
--divides (p^ek - 1)
tauQGor = (Rk, ek, fk, t1) -> 
(
     Sk := ambient Rk;
     pp := char Sk;
     L1 := divideFraction(pp,t1); --this breaks up t1 into the pieces we need
     hk := findQGorGen(ek, Rk); --the term in the test ideal
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
          	I2 = ethRootRingElements(pPow, numerator((pp^pPow - 1)/(pp^ek - 1)), hk, I1 )
		)
		else I2 = I1
     )
     else (
	  	I1 = ascendIdeal(ek, hk, Jk);
	  	if (pPow != 0) then (
	  		I2 = ethRootRingElements( pPow, {numerator((pp^pPow - 1)/(pp^ek - 1)), a2}, {hk, fm}, I1 )
	  	)
	  	else I2 = I1
     );

	  
     sub(I2, Rk)*ideal(fk^k2)
)

--Computes tau(Rk,fk^tk), assuming Gorenstein rings
tauGor = ( R, f, t ) -> tauQGor( R, 1, f, t )


----------------------------------------------------------------
--************************************************************--
--Test ideals for non-principal ideals                        --
--************************************************************--
----------------------------------------------------------------

--takes an ideal, forms the rees algebra, and returns the rees algebra in two ways, first with flattened variables and the second without
flattenedReesAlgebra = I -> 
(
	S1 := reesAlgebra I;
	J1 := ideal S1;
	tempMonoid := S1.FlatMonoid;
	k := coefficientRing (ring I);
	S2 := k tempMonoid;
	J2 := sub(J1, S2);	
	(S2/J2, S1)
)

needsPackage "BGG"; --we'll be pushing forward...

needsPackage "Divisor";

tauNonPrincipalAOverPEPoly = {Verbose=> false}>> o -> (I1, a1, e1) -> ( -- computes \tau(I^{a/p^e}) for I an ideal in a polynomial ring
	if ( not(codim(I1) > 1)) then error "We can only handle ideals of codimension > 1 at this time.";

	R1 := ring I1;
	p1 := char R1;

    -- Skoda's theorem
    n:= numgens I1;
    mSkoda := 0;
    if (a1 > p1^e1*n) then (
        mSkoda = floor(a1/p1^e1) + 1 - n;
        a1 = a1 - mSkoda *p1^e1;
    );


	reesList := flattenedReesAlgebra I1;
	A1 := reesList#0; --this one has flattened variables
	A2 := reesList#1;
 	irrIdeal := sub(ideal(first entries vars A1), A1);
 	singLocus := ideal singularLocus (A1);
 	
 	IRees := sub(I1, A2);
 	
 	canList := canonicalIdeal(A1, FullMap => true);
    --canonicalIdeal doesn't return a list though ?

 	canIdeal := canList#0;
 	canMap := canList#1;
 	
 	paraTest := testModule(A1, canIdeal); 
 		
 	newMap := map(A1^1/(paraTest#0), canList#2, matrix(canMap));
 	newKer := (ker newMap)**A2; --this is the parameter test submodule of the canonical module  

	flag := false;
	i1 := e1;
		ascend := I1; --dummy variables for checking whether we are done
	descend := ideal(sub(1, R1)); --dummy variables for checking whether we are done
	
	while (flag == false) do (
		if (o.Verbose == true) then (print  "Ascending ideal");
		ascend = ethRoot(i1, a1*p1^(i1-e1), I1);
		if (o.Verbose == true) then (print ascend);
		
		flag = isSubset(descend, ascend);
		if (o.Verbose == true) then (print "flag"; print flag);
		if (flag == false) then (
			
			myDirectImage := HH_0(directImageComplex(IRees^(a1*p1^(i1-e1))*newKer, Regularity=>(10+a1))); 	
 	
		 	directIdeal := embedAsIdeal(myDirectImage**R1); --probably the tensor product is unnecessary
 			if ( codim(directIdeal)==1) then error "This function produced a codimension 1 ideal.";
 	
 			descend = ethRoot(i1,directIdeal);
 			if (o.Verbose == true) then (print  "Descending ideal"; print descend)
		);
		

		flag = isSubset(descend, ascend);
		
		--the following should be removed eventually, it is only here for debug purposes
		if ((flag == true) and (isSubset(ascend, descend)==false)) then error "Major error detected";
		i1 = i1+1;
		if (o.Verbose==true) then (print "Loop complete, continue?"; print (not flag) );
	);
	
	ascend*I1^mSkoda
)


