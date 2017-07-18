---compute a generating morphism for H^(dim R -i)_I(R)
---the output is (A,U) where U:coker(A) -> F(coker A) is the generating morphism
generatingMorphism= (i,I) ->(
	local F1; local K; local C;
	local F1p; local Kp; local Cp;
	R:=ring(I);
	Ip:=frobenius( I );
	M:=coker I;
	Mp:=coker Ip;
	resM:=res M;
	resMp:=res Mp;
	f:=inducedMap(M,Mp);
	resf:=res f;
---Compute the generating map Ext^i(R/I,R) -> F Ext^i(R/I,R) [[The top i doesnt work: the zero matrix is not given]]
	G:=resf#i; G=transpose(G);
	F0:=(resM.dd)#(i); F0=transpose(F0);
	if (resM.dd)#?(i+1) then
	{
		F1=(resM.dd)#(i+1); F1=transpose(F1);
		K=ker F1;
	} else
	{
		K=target(F0);
	};
	temp1:=substitute(gens K,R);
	if (temp1==0) then C=coker(F0) else C=subquotient(substitute(gens K,R),F0);
---C=subquotient(substitute(gens K,R),F0);
	C1:=prune(C);
	h:=C1.cache.pruningMap;
	generatingMorphism0:=G*gens(K)*matrix(entries h);
	F0p:=(resMp.dd)#(i); F0p=transpose(F0p);
	if (resMp.dd)#?(i+1) then
	{
		F1p=(resMp.dd)#(i+1); F1p=transpose(F1p);
		Kp=ker F1p;
	} else
	{
		Kp=target(F0p);
	};
	temp1=substitute(gens Kp,R);
	if (temp1==0) then Cp=coker(F0p) else Cp=subquotient(substitute(gens Kp,R),F0p);
--Cp=subquotient(gens Kp,F0p);
	C1p:=prune(Cp);
	hp:=C1p.cache.pruningMap;
	A0:=gens(Kp)*matrix(entries hp); A:=A0| F0p;
	gbA:=gb(A, ChangeMatrix => true) ;
	B:=generatingMorphism0// A;
--- Now generatingMorphism0=A*B
	k:=rank source A0;
	(relations(C1), submatrix(B,toList(0..(k-1)),))
)



-- Produce a sequence of maps Ext^i(R/I,R) ->  Ext^i(R/I^{[p]},R) induced
-- by the surjections R/I^{[p]} -> R/I
-- for i=1..pdim(coker I)
-- The output consists of a sequence of pairs (A,B) where the induced maps are
-- B: coker A -> coker A^{[p]}
findGeneratingMorphisms = (I) ->
(
	local i;
	Ip:=frobenius( I );
	M:=coker I;
	Mp:=coker Ip;
	resM:=res M;
	resMp:=res Mp;
	f:=inducedMap(M,Mp);
	resf:=res f;
	resLength:=length(resM);
	answer:=();
	apply(1..resLength, i->
	{
		G:=resf#i; G=transpose(G);
		F1:=(resM.dd)_(i+1); F1=transpose(F1);
		F0:=(resM.dd)_(i); F0=transpose(F0);
		K:=ker F1;
		C:=subquotient(gens K,F0);
		C1:=prune(C);
		h:=C1.cache.pruningMap;
--
		generatingMorphism0:=G*gens(K)*matrix(entries h);
		F1p:=(resMp.dd)_(i+1); F1p=transpose(F1p);
		F0p:=(resMp.dd)_(i); F0p=transpose(F0p);
		Kp:=ker F1p; 
		Cp:=subquotient(gens Kp,F0p);
		C1p:=prune(Cp);
		hp:=C1p.cache.pruningMap;
--
		A0:=gens(Kp)*matrix(entries hp); 
		A:=A0| F0p;
		gbA:=gb(A, ChangeMatrix => true) ;
		B:=generatingMorphism0// A;
--- Now generatingMorphism0=A*B
		k:=rank source A0;
		generatingMorphism:=submatrix(B,toList(0..(k-1)),);
		answer=append(answer, (C1,generatingMorphism));
---		print(generatingMorphism);
	});
answer
)

