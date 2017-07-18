----------------------------------------------------------------
--************************************************************--
--Functions for computing compatibly split ideals             --
--************************************************************--
----------------------------------------------------------------

-----------------------------------------------------------------------


--- Start of MK ---------------------------------------------------------------------------------------------------

-- FIND IDEALS COMPATIBLE WITH A GIVEN NEAR-SPLITTING
-- This is an implementation of the algorithm described in
-- Moty Katzman and Karl Schwede's paper 
-- "An algorithm for computing compatibly Frobenius split subvarieties"
-- J. Symbolic Comput. 47 (2012), no. 8, 996\961008. 

----------------------------------------------------------------------------------------


--- Input:
---   	an element u of the polynomial ring R OVER A PRIME FIELD.
--- Output:
---	A list of all prime ideals P such that
---	(a) u P \subseteq P^{[p]}, and
---	(b) the action of uT on the the annihilator of P on the injective hull of the residue field of R 
---	is not the zero Frobenius map.


findAllCompatibleIdeals = (u) ->(
	L:={}; R:=ring u; p:=char R;
	P:=ideal(0_R);
	J:=ethRoot(1,ideal(u));
	t:=1_R % (gens J);
	if (t != 0_R) then print("*** WARNING *** Frobenius action has nilpotent elements");
	findAllCompatibleIdealsInnards (u,L,P)
)



findAllCompatibleIdealsInnards = (u,L,P) ->(
	R:=ring u;
	p:=char R;
	local tau;
	local Plist;
	P1:=frobenius(P);
	C1:=ideal((singularLocus(P)).relations);
	---tau=ideal mingens star(C1,u,1) ; ---OLD VERSION
	tau=ideal mingens minimalCompatible (1,u,C1);
	Plist=minimalPrimes tau;
	local Q;
	local T;
	apply(Plist, Q->
	{
		f:= any(L,T -> T == Q);
---print(L,Q,f);
		if (not f) then
		{
			L=append(L,Q);
			L=unique(L | findAllCompatibleIdealsInnards(u,L,Q));
		};
	});
---
	C2:=(P1+ideal(u)):(P1:P);
---	JB:=C1*C2; ---MK
---print(mingens P, mingens JB);
---tau=ideal mingens star(C2,u,1) ;  --- OLD VERSION
	tau=ideal mingens minimalCompatible (1,u,C2);
	Plist=minimalPrimes tau;
	local Q;
	local T;
	apply(Plist, Q->
	{
		f:= any(L,T -> T == Q);
	---print(L,Q,f);
		if (not f) then
		{
			L=append(L,Q);
			L=unique(L | findAllCompatibleIdealsInnards(u,L,Q));
		};
	});
	---
	L
)
