----------------------------------------------------------------
--************************************************************--
--Functions for computing parameter test modules and ideals   --
--************************************************************--
----------------------------------------------------------------


--This function computes the parameter test module of a ring, it returns it as a submodule of a canonical ideal.
--this is a slightly modified function originally written by Moty Katzman for "Parameter test ideals of Cohen Macaulay rings"
--it returns the lift of the canonical module to the ambient ring
needsPackage "Divisor";

canonicalIdeal = method();

canonicalIdeal(Ring) := (R1) -> (
    S1 := ambient R1;
	I1 := ideal R1;
	dR := dim R1;
	dS := dim S1;
	varList := first entries vars S1;
	degList := {};
	if (#(degree(varList#0)) == 1) then (
		degList = apply(varList, q -> (degree(q))#0); )
	else (
		degList = apply(varList, q -> (degree(q))); );
	M1 := (Ext^(dS - dR)(S1^1/I1, S1^{-(sum degList)}))**R1;
	moduleToIdeal(M1)
)


--the following function computes the u of a canonical ideal in a polynomial ring
--it uses previous work of Katzman
finduOfIdeal = method();

finduOfIdeal(Ideal, Ideal) := (defIdeal, canIdeal) -> (
	Ip := frobeniusPower(1, defIdeal);
	tempIdeal := intersect( (frobeniusPower(1,canIdeal)) : canIdeal, Ip : defIdeal );
	
	M1 := compress ((gens tempIdeal)%(gens Ip));
	first first entries M1
);

--computes the parameter test submodule of a given ring.  It outputs the parameter test module (as an ideal), it then outputs the canonical module (as an ideal), and finally it outputs the term u used as the action on the ideal
paraTestModuleAmbient = method();

paraTestModuleAmbient (Ring) := (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal := sub(canonicalIdeal(R1), S1) + I1;
	
	J1 := (findTestElementAmbient(R1));
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(I1, canIdeal); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(1, u1, tau0);
	
	(sub(tauOut, R1), sub(canIdeal, R1), 
	u1)
)

paraTestModuleAmbient (Ring, Ideal) := (R1, canIdeal) -> (--it expects the canonical ideal to be lifted to the ambient ring
	S1 := ambient R1;
	I1 := ideal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(I1, canIdeal); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(1, u1, tau0);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

--computes the parameter test ideal of an ambient ring
paraTestIdealAmbient = (R1) -> (
	tempList := paraTestModuleAmbient(R1);
	(tempList#0) : (tempList#1)
)
paraTestModule = method(Options=>{AscentCount=>false})
--this computes the parameter test module \tau(R, f^t).  It does not assume that R is a polynomial ring.
paraTestModule(QQ, RingElement) := o -> (t1, fk) -> ( --maintained by Karl
	R1 := ring fk;
	S1 := ambient R1;
	f1 := sub(fk, S1);
	I1 := ideal R1;
	pp := char R1;
	funList := divideFraction(pp, t1);
	
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	
--	tempList := paraTestModuleAmbient(R1);
--	tauAmb := sub(tempList#0, S1);
--	omegaAmb := sub(tempList#1, S1);
--	u1 := tempList#2;

	omegaAmb := sub(canonicalIdeal(R1), S1) + I1;
	J1 := findTestElementAmbient(R1)*omegaAmb;
	u1 := finduOfIdeal(I1, omegaAmb);

	uPower := 1;
	if (cc != 0) then
		uPower = floor((pp^cc-1)/(pp-1));
	firstTau := J1;
	local tempList;
	ascendingCount := 0;
--	assert false;
	if (cc != 0) then	
--??? REORDER PARAMETERS
		if (o.AscentCount == false) then (firstTau = ascendIdeal( (aa, uPower), cc, (f1, u1), J1*ideal(f1^(pp^bb*ceiling(t1))) ))
		else (tempList = ascendIdeal(  cc, (aa, uPower), (f1, u1), J1*ideal(f1^(pp^bb*ceiling(t1))), AscentCount=>true);
			firstTau = tempList#0;
			ascendingCount = tempList#1;
		)
--		firstTau = ascendIdeal(cc, f1^aa*u1^(uPower), J1*ideal(f1^(aa)))

	else 
--		firstTau = ascendIdeal(1, u1^(uPower), J1)*ideal(f1^aa);
		firstTau = ascendIdeal( 1, uPower, u1, J1);
			
	secondTau := firstTau;
	if (bb != 0) then
		secondTau = ethRootRingElements(bb, floor((pp^bb-1)/(pp-1)), u1, firstTau); --??? REORDER PARAMETERS

	if (o.AscentCount == false) then (sub(secondTau, R1), omegaAmb, u1) else (sub(secondTau, R1), omegaAmb, u1, ascendingCount)
)

--this computes the parameter test module \tau(R, f^t).  It does not assume that R is a polynomial ring.
paraTestModule(QQ, RingElement, Ideal, RingElement) := o -> (t1, fk, omegaAmb, u1) -> ( --maintained by Karl
	R1 := ring fk;
	S1 := ambient R1;
	f1 := sub(fk, S1);
	I1 := ideal R1;
	pp := char R1;
	funList := divideFraction(pp, t1);
	
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	
--	tempList := paraTestModuleAmbient(R1);
--	tauAmb := sub(tempList#0, S1);
--	omegaAmb := sub(tempList#1, S1);
--	u1 := tempList#2;

	J1 := findTestElementAmbient(R1)*omegaAmb;

	uPower := 1;
	if (cc != 0) then
		uPower = floor((pp^cc-1)/(pp-1));
	firstTau := J1;
	local tempList;
	ascendingCount := 0;
--	assert false;
	if (cc != 0) then	
--??? REORDER PARAMETERS
		if (o.AscentCount == false) then (firstTau = ascendIdeal(cc, {aa, uPower}, (f1, u1), J1*ideal(f1^(pp^bb*ceiling(t1))) ))
		else (tempList = ascendIdeal(  cc, (aa, uPower), (f1, u1), J1*ideal(f1^(pp^bb*ceiling(t1))), AscentCount=>true);
			firstTau = tempList#0;
			ascendingCount = tempList#1;
		)
--		firstTau = ascendIdeal(cc, f1^aa*u1^(uPower), J1*ideal(f1^(aa)))

	else 
--		firstTau = ascendIdeal(1, u1^(uPower), J1)*ideal(f1^aa);
		firstTau = ascendIdeal( 1, {uPower},  {u1}, J1);
			
	secondTau := firstTau;
	if (bb != 0) then
		secondTau = ethRootRingElements(bb, floor((pp^bb-1)/(pp-1)), u1, firstTau); --??? REORDER PARAMETERS

	if (o.AscentCount == false) then (sub(secondTau, R1), omegaAmb, u1) else (sub(secondTau, R1), omegaAmb, u1, ascendingCount)
)


--****************************************************
--*****Karl is starting a rewrite of some of this*****
--****************************************************

--this function finds the generators of the intersection of 
--J^{[p]} : J and I^{[p]} : I where I is the defining ideal and J is the canonical
--ideal lifted to the ambient ring (in a maximal way).
findusOfIdeal = (defIdeal, canIdeal) -> (
	Ip := frobeniusPower(1, defIdeal);
	tempIdeal := intersect( (frobeniusPower(1,canIdeal)) : canIdeal, Ip : defIdeal );
	
	M1 := compress ((gens tempIdeal)%(gens Ip));
	first entries M1
)

testModule = method(); --a rewritten function to construct the (parameter) test (sub)module of a given ring.  
                       --it returns two ideals and an element.  
                       --The first ideal is an ideal isomorphic to the test module and the
                       --and the second is an ideal isomorphic to the canonical module, in which the parameter
                       --resides.  The locus where these two ideals differ (plus the non-CM locus) is the
                       --locus where the ring does not have rational singularities.
                       --the final element is the element of the ambient polynomial ring which is used to induce
                       --the canonical trace map

testModule(Ring) := (R1) -> (
    S1 := ambient R1;
	I1 := ideal R1;
    J1 := sub(canonicalIdeal(R1), S1);
    C1 := findTestElementAmbient(R1);
    
    u1 := findusOfIdeal(I1, J1+I1);
    tau := I1;
    if (#u1 > 1) then(
        print "Multiple trace map for omega generators.  Using them all.";
        j := 0;
        while (j < #u1) do (
            tau = tau + ascendIdeal(1, u1, C1*J1*R1);
            j = j+1;
        );
    )
    else (
        u1 = u1#0;
        tau = ascendIdeal(1, u1, C1*J1*R1);
    );

    (sub(tau, R1), sub(J1, R1), u1)
);


