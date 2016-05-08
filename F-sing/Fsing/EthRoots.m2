--*************************************************
--*************************************************
--This file is used for doing the [1/p^e] operation
--in the sense of Blickle-Mustata-Smith.
--This operation is also called I_e in Katzman or 
--simply the image of 
--M \subseteq \Hom_R(R^{1/p^e}, R) -> R
--under evaluation at 1.
--*************************************************
--*************************************************

ethRoot = method(); --- MK


--Computes I^{[1/p^e]}, we must be over a perfect field. and working with a polynomial ring
--This is a slightly stripped down function due to Moty Katzman, with some changes to avoid the
--use(Rm) which is commented out below
--The real meat of the function is ethRootInternal, this function just looks for a certain error and calls 
--the other function depending on that error.
ethRoot(Ideal,ZZ) := (Im,e) -> (
     J := Im;
     success := false;
     count := 0;
     try J = ethRootInternal(J,e) then success = true else (
--     print "blew a buffer";
	 while(count < e) do (	 	
	      J = ethRootInternal(J, 1);
	      count = count + 1
	 )
     );
     J
)

ethRoot(MonomialIdeal, ZZ ) := ( I, e ) ->
(
     R := ring I;
     p := char R;
     G := I_*;
     if #G == 0 then ideal( 0_R ) else ideal( apply( G, j -> R_((exponents(j))#0//p^e )))
)

--This tries to compute (f^a*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
--We should probably just use it.  It relies on the fact that (f^(ap+b))^{[1/p^2]} = (f^a(f^b)^{[1/p]})^{[1/p]}.
ethRootSafe = method(); 

ethRootSafe( RingElement, Ideal, ZZ, ZZ ) := ( f, I, a, e ) -> (
	R1 := ring I;
	p1 := char R1;
	
	aRem := a%(p1^e);
	aQuot := floor(a/p1^e);
	
	expOfA := basePExpMaxE(aRem,p1,e); --this gives "a base p", with the left-most term the smallest "endian".
	
	IN1 := I;
	
	if (e > 0) then (
		IN1 = IN1*ideal(f^(expOfA#0));
		IN1 = ethRoot(IN1, 1);
		i := 1;
	
		while(i < #expOfA) do (
			IN1 = ethRoot( IN1*ideal(f^(expOfA#i)), 1);
			i = i + 1;
		)
	);
	IN1*ideal(f^(aQuot))
)

ethRootSafe( RingElement, ZZ, ZZ ) := ( f, a, e ) -> 
    ethRootSafe( f, ideal( 1_(ring f) ), a, e )

--This tries to compute (f1^a1*f2^a2*...fk^ak*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
ethRootSafeList = method();

ethRootSafeList( List, Ideal, List, ZZ ) := ( elmtList, I1, aList, e1 ) -> (
	   R1 := ring I1;
        p1 := char R1;
        
        aListRem := apply(aList, z1 -> z1%(p1^e1) );
        aListQuot := apply(aList, z1 -> floor(z1/p1^e1) );
        
        expOfaList := apply(aListRem, z1-> basePExpMaxE(z1, p1, e1) );
        
        aPowerList := apply(elmtList, expOfaList, (f1, z1) -> f1^(z1#0));
        
        IN1 := I1*ideal(fold(times, aPowerList));
        if (e1 > 0) then (
                IN1 = ethRoot(IN1, 1);
                i := 1;
                while(i < e1) do (
                        aPowerList = apply(elmtList, expOfaList, (f1, z1) -> f1^(z1#i));
                        IN1 = ethRoot( IN1*ideal(fold(times, aPowerList)), 1);
                        i = i + 1;
                )
        );
        aPowerList = apply(elmtList, aListQuot, (f1, z1) -> f1^z1);
        IN1*ideal(fold(times, aPowerList))
)

ethRootSafeList( List, List, ZZ ) := ( F, a, e ) ->
    ethRootSafeList( F, ideal( 1_( ring( F#0 ) ) ), a, e )
	
ethRoot(RingElement, Ideal, ZZ, ZZ) := (f, I, a, e) -> ethRootSafe (f, I, a, e) ---MK

ethRootInternal = (I,e) -> (
     if (not isIdeal(I)) then (error "ethRoot: Expected first argument to be an ideal.");
     if (not e >= 0) then (error "ethRoot: Expected second argument to be a nonnegative integer.");
     R:=ring(I); --Ambient ring
     if (class R =!= PolynomialRing) then (error "ethRoot: Expected an ideal in a PolynomialRing.");
     p:=char(R); --characteristic
     kk:=coefficientRing(R); --base field
     if ((kk =!= ZZ/p) and (class(kk) =!= GaloisField)) then (error "ethRoot: Expected the coefficient field to be ZZ/p or a GaloisField.");
     var:=R_*; --the variables (henceforth denoted X_i)
     n:=#var; --number of variables
     Y:=local Y;
     newvar := var | toList(Y_1..Y_n);
     S:=kk(monoid[newvar, MonomialOrder=>ProductOrder{n,n},MonomialSize=>64]);
         -- brand new ring, with a variable Y_i for each X_i
     J:=matrix {toList apply(n, i->newvar#(n+i)-newvar#(i)^(p^e))}; 
         -- J = (Y_i-X_i^(p^e)) 
     rules:=toList apply(n, i->newvar#(n+i)=>substitute(var#(i),S)); 
         -- {Y_i =>X_i} 
     G:=first entries compress((gens substitute(I,S)) % J);
     	 -- replaces X_i^(p^e) with Y_i 
     L:=sum(G,t->ideal((coefficients(t,Variables=>var))#1));	 
     L=first entries mingens L;
     L=apply(L, t->substitute(t,rules));
     q:=kk#order;
     if (q > p) then 
     (
	 a:=(gens kk)#0;
	 e0:=floorlog(p,q); 
	 root:=a^(p^(e0-(e%e0)));
     	 L=apply(L,t->substitute(t,a=>root));
	     -- substitute generator of kk with its p^e-th root
     );
     substitute(ideal L,R)
)


--A short version of ethRoot
eR = (I1,e1)-> (ethRoot(I1,e1) )

---------------------------------------------------------------------------------------
--- The following code was written in order to more quickly compute eth roots of (f^n*I)
--- It is used in fancyEthRoot
----------------------------------------------------------------------------------------
--- Find all ORDERED partitions of n with k parts
allPartitions = (n,k)->
(
	PP0:=matrix{ toList(1..k) };
	PP:=mutableMatrix PP0;
	allPartitionsInnards (n,k,PP,{})
)

allPartitionsInnards = (n,k,PP,answer)->
(
	local i;
	if (k==1) then 
	{
		PP_(0,k-1)=n;
		answer=append(answer,first entries (PP));
	}
	else
	{
		for i from 1 to n-(k-1) do
		{
			PP_(0,k-1)=i;
			answer=allPartitionsInnards (n-i,k-1,PP,answer)	;	
		};
	};
	answer
)




--- write n=a*p^e+a_{e-1} p^{e-1} + \dots + a_0 where 0\leq e_j <p 
baseP1 = (n,p,e)->
(
	a:=n//(p^e);
	answer:=1:a;
	m:=n-a*(p^e);
	f:=e-1; 
	while (f>=0) do
	{
		d:=m//(p^f);
		answer=append(answer,d);
		m=m-d*(p^f);
		f=f-1;
	};
	answer
)	


fancyEthRoot = (I,m,e) ->
(
	G:=first entries mingens I;
	k:=#G;
	P:=allPartitions(m,k);
	R:=ring(I);
	p:=char(R);
	answer:=ideal(0_R);
	apply(P, u->
	{
	---print("Partition: ",u);
		a:=ideal(1_R);
		U:=apply(u, v->baseP1(v,p,e));
		for i from 0 to e do
		{
			j:=e-i;
			g:=1_R;
			for l from 0 to k-1 do g=g*(G#l)^((U#l)#j); 
			a=ideal(g)*a;
			if (i<e) then a=ethRoot(a ,1);
---print(g,answer);
		};
		answer=answer+a;
	});
	ideal(mingens(answer))
)

ethRoot (Ideal, ZZ, ZZ) := (I,m,e) -> fancyEthRoot (I,m,e)  --- MK

ethRoot( RingElement, ZZ, ZZ ) := ( f, a, e ) -> ethRoot( ideal( f ), a, e )

--Computes I^{[1/p^e]}, we must be over a perfect field. and working with a polynomial ring
--This is a slightly stripped down function due to Moty Katzman, with some changes to avoid the
--use(Rm) which is commented out below
--The real meat of the function is ethRootInternal, this function just looks for a certain error and calls 
--the other function depending on that error.
ethRoot(Ideal,ZZ) := (Im,e) -> (
     J := Im;
     success := false;
     count := 0;
     try J = ethRootInternal(J,e) then success = true else (
--     print "blew a buffer";
	 while(count < e) do (	 	
	      J = ethRootInternal(J, 1);
	      count = count + 1
	 )
     );
     J
)
