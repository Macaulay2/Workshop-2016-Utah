----------------------------------------------------------------
--************************************************************--
--Functions for computing parameter test modules and ideals   --
--************************************************************--
----------------------------------------------------------------


--This function computes the parameter test module of a ring, it returns it as a submodule of a canonical ideal.
--this is a slightly modified function originally written by Moty Katzman for "Parameter test ideals of Cohen Macaulay rings"
--it returns the lift of the canonical module to the ambient ring

canonicalIdeal ={FullMap=> false} >> o -> (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	d1 := (dim S1) - (dim R1);
	local answer2;
	
	degShift := sum degrees S1;
	myExt := prune( Ext^d1(S1^1/I1, S1^{-degShift}));
	canModuleMatrix := relations(myExt);
	
	answer:=0;
	s1:=syz transpose substitute(canModuleMatrix,R1);
	s2:=entries transpose s1;
	--use S1;
	apply(s2, t->
	{
		s3:=substitute(syz gens ideal t,S1);
---		print(s3%canModuleMatrix);
		if ((s3%canModuleMatrix)==0) then
		{
			answer2 = t;
			answer=substitute(mingens ideal t,S1);
			break;
		};
	});
	
	
	
	if (o.FullMap == true) then (ideal answer, map(R1^1, myExt**R1, matrix {answer2}), (myExt**S1^{-degShift})**R1) else ideal answer
)

--moduleToIdeal = (M1, R1) -> (--turns a module to an ideal of a ring, it returns the lift of the ideal to the ambient ring
--	S1 := ambient R1;
---	myMatrix := substitute(relations prune M1, S1);
--	
--	answer:=0;
--	s1:=syz transpose substitute(myMatrix,R1);
--	s2:=entries transpose s1;
--	
--	apply(s2, t->
--	{
--		s3:=substitute(syz gens ideal t,S1);
---		print(s3%canModuleMatrix);
--		if ((s3%myMatrix)==0) then
--		{
--			answer=substitute(mingens ideal t,S1);
--			break;
--		};
--	});
--	ideal answer	
--)

--the following function computes the u of a canonical ideal in a polynomial ring
--it uses previous work of Katzman
finduOfIdeal = (canIdeal, defIdeal) -> (
	Ip := frobeniusPower(1, defIdeal);
	tempIdeal := intersect( (frobeniusPower(1,canIdeal)) : canIdeal, Ip : defIdeal );
	
	M1 := compress ((gens tempIdeal)%(gens Ip));
	first first entries M1
)

--computes the parameter test submodule of a given ring.  It outputs the parameter test module (as an ideal), it then outputs the canonical module (as an ideal), and finally it outputs the term u used as the action on the ideal
paraTestModuleAmbient = method();

paraTestModuleAmbient (Ring) := (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal := canonicalIdeal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(1, u1 tau0);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

paraTestModuleAmbient (Ring, Ideal) := (R1, canIdeal) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(1, u1, tau0);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

paraTestModuleAmbient (Ring, Ideal) := (R1, canIdeal) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(1, u1, tau0);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

--computes the parameter test ideal of an ambient ring
paraTestIdealAmbient = (R1) -> (
	tempList := paraTestModuleAmbient(R1);
	(tempList#0) : (tempList#1)
)

--this computes the parameter test module \tau(R, f^t).  It does not assume that R is a polynomial ring.
paraTestModule ={AscentCount=>false} >> o -> (fk, t1) -> ( --maintained by Karl
	R1 := ring fk;
	S1 := ambient R1;
	f1 := sub(fk, S1);
	I1 := ideal R1;
	pp := char R1;
	funList := divideFraction(t1, pp);
	
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	
--	tempList := paraTestModuleAmbient(R1);
--	tauAmb := sub(tempList#0, S1);
--	omegaAmb := sub(tempList#1, S1);
--	u1 := tempList#2;

	omegaAmb := canonicalIdeal(R1);
	J1 := findTestElementAmbient(R1)*omegaAmb;
	u1 := finduOfIdeal(omegaAmb, I1);

	uPower := 1;
	if (cc != 0) then
		uPower = floor((pp^cc-1)/(pp-1));
	firstTau := J1;
	local tempList;
	ascendingCount := 0;
--	assert false;
	if (cc != 0) then	
--??? REORDER PARAMETERS
		if (o.AscentCount == false) then (firstTau = ascendIdealSafeList( J1*ideal(f1^(pp^bb*ceiling(t1))), (f1, u1), (aa, uPower), cc))
		else (tempList = ascendIdealSafeList( J1*ideal(f1^(pp^bb*ceiling(t1))), (f1, u1), (aa, uPower), cc, AscentCount=>true);
			firstTau = tempList#0;
			ascendingCount = tempList#1;
		)
--		firstTau = ascendIdeal(cc, f1^aa*u1^(uPower), J1*ideal(f1^(aa)))
		--I should write an ascendIdealSafe that works for multiple elements raised to powers...	
	else 
--		firstTau = ascendIdeal(1, u1^(uPower), J1)*ideal(f1^aa);
		firstTau = ascendIdealSafe(J1, u1, uPower, 1); --??? REORDER PARAMETERS
			
	secondTau := firstTau;
	if (bb != 0) then
		secondTau = ethRootRingElements(u1, firstTau, floor((pp^bb-1)/(pp-1)) , bb); --??? REORDER PARAMETERS

	if (o.AscentCount == false) then (sub(secondTau, R1), omegaAmb, u1) else (sub(secondTau, R1), omegaAmb, u1, ascendingCount)
)

