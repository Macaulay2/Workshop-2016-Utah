-- changes in argument order:

-- INTERNAL: 

-- ethRoot: (I,e) -> (e,I), (f, I, a, e) -> (e,a,f,I), (f, I, a, e) -> (e,a,f,I), (Matrix, ZZ) -> (ZZ, Matrix), 
--    	    (Ideal, ZZ, ZZ),  (I,m,e) -> ( e, m, I ) 

-- eR: gone

-- ethRootSafe: (f,I,a,e) -> (e,a,f,I), (f,a,e) -> (e,a,f)
		
-- ethRootSafeList: (fList,I,aList,e) -> (e,aList,fList,I), (fList,aList,e) -> (e,aList,fList)

-- mEthRoot, mEthRootOneElement: (A,e) -> (e,A)

-- fancyEthRoot (I,m,e) -> (e, m, I)

-- EXTERNAL: basePExp

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

-------------------------------------------------------
---------- List of functions in this file -------------
-----------------(as of 2016-05-25)--------------------
-------------------------------------------------------
-- ethRoot
-- getFieldGenRoot 
-- ethRootMonStrat 
-- ethRootSubStrat 
-- ethRootSafe
-- ethRootSafeList 
-- ascendIdeal 
-- ascendIdealSafe 
-- ascendIdealSafeList 
-- getExponents
-- mEthRootOfOneElement
-- mEthRoot 
-- Mstar 
-------------------------------------------------------
-------------------------------------------------------
-------------------------------------------------------

needs "BasicFunctions.m2" -- maybe this should be removed
                          -- when we publish this package

ethRoot = method(Options => {EthRootStrategy => Substitution});
--ethRoot takes two strategy options: Substitution and MonomialBasis
--The second strategy seems to generally be faster for computing I^[1/p^e] when e = 1, especially for polynomials of
--high degree, but slower for larger e. 
-- Dan: I wonder if this is because getFieldGen is not optimized? It's called many times per
-- generator of the ideal in the monomial strategy. Though I see it's also called for the 
-- substitution strategy...

ethRoot ( ZZ, Ideal ) := Ideal => opts -> (e,I) -> (
    if (not e >= 0) then (error "ethRoot: Expected second argument to be a nonnegative integer.");
    R := ring I;
    if (class R =!= PolynomialRing) then (error "ethRoot: Expected an ideal in a PolynomialRing.");
    p := char R;
    k := coefficientRing(R);
    if ((k =!= ZZ/p) and (class(k) =!= GaloisField)) then (error "ethRoot: Expected the coefficient field to be ZZ/p or a GaloisField.");

    q := k#order;
    --Gets the cardinality of the base field.
    G := I_*;
    --Produces a list of the generators of I.
    if #G == 0 then ideal(0_R) 
    else if opts.EthRootStrategy == MonomialBasis then (
	L := sum( apply( G, f -> ethRootMonStrat(e,p,q,k,f,R) ) );
	L = first entries mingens L;
	ideal(L)
	)
    else ethRootSubStrat(e,p,q,k,I,R)  
)

-----------------------------------------------------------------------------

ethRoot ( ZZ, MonomialIdeal ) := Ideal => opts -> (e,I) -> (
     R := ring I;
     p := char R;
     G := I_*;
     if #G == 0 then ideal( 0_R ) else ideal( apply( G, f -> R_((exponents(f))#0//p^e )))
)
-----------------------------------------------------------------------------

ethRoot ( ZZ, ZZ, RingElement, Ideal ) := opts -> ( e, a, f, I ) -> ethRootSafe ( e, a, f, I ) ---MK

-----------------------------------------------------------------------------

ethRoot ( ZZ, ZZ, RingElement ) := opts -> ( e, a, f ) -> ethRootSafe ( e, a, f ) ---MK

-----------------------------------------------------------------------------

ethRoot ( ZZ, ZZ, Ideal ) := opts -> ( e, m, I ) -> fancyEthRoot( e, m, I )  --- MK

-----------------------------------------------------------------------------

ethRoot ( RingElement, ZZ, ZZ ) := opts -> ( f, a, e ) -> ethRoot( ideal( f ), a, e )

-----------------------------------------------------------------------------

ethRoot ( ZZ, Matrix ) := opts -> (e, A) -> mEthRoot (e,A)  --- MK

-----------------------------------------------------------------------------
-- MACHINERY
-----------------------------------------------------------------------------

getFieldGenRoot = (e,p,q,k) -> (
    s := floorLog(p,q);
    -- Gets the exponent s such that q = p^s.
    a := (gens k)#0;
    a^(p^(s-e%s))
    -- Gets the p^e-th root of the cyclic generator a for the field extension k 
    -- over ZZ/p.  If 1,a,..,a^t is a basis for k over ZZ/p and 
    -- c = c_0 + c_1a + .. + c_ta^t in k, then replacing a with its p^e-th root 
    -- in the preceding expansion using substitute(c,a => getFieldGenRoot(e,p,q,k))
    -- yields the p^e-th root of c.
)


-----------------------------------------------------------------------------

ethRootMonStrat = (e,p,q,k,f,R) -> (
    -- e = exponent, p = prime, q = size of coeff field, k = coeff field, 
	-- f = a generator of the ideal in question, R = the ring
	-- to use this strategy to find the p^eth root of an ideal, you need to apply this
	-- function to each generator of the ideal and sum the results. 
	-- maybe this should just return the ideal though? I guess it's an internal
	-- function, so it doesn't matter.
    expDecomp := apply(exponents(f),exponent->{coefficient(R_exponent,f)*R_(exponent //p^e),exponent%p^e});
    --Gets the exponent vectors of each monomial X^u of the polynomial f, and associates to u the two-element list whose
    --first entry is cX^v and second entry is w, where c is the coefficient of X^u in f and u = p^e*v + w. 
    if q > p then (
	substRule := ( (gens k)#0 => getFieldGenRoot(e,p,q,k) );
	expDecomp = apply( expDecomp, pair -> { substitute( pair#0, substRule ), pair#1 } );
    );
    remainders := partition(x-> x#1, expDecomp);
    --Sorts the preceding list of two-element lists into a hash table with keys the remainder w of the exponent vector.
    --The value of each key is a list of two-element lists {cX^v,w} with the same remainder.
    remainders = applyValues(remainders,v->apply(v,w->(w#0)));
    --Forgets the second entry of each two-element list in the preceding hash table.
    remainders = applyValues(remainders,v->sum(v));
    --Adds together all the terms for each key w in the hash table to get the coefficient of the basis monomial X^w
    --for R over R^(p^e).
    return ideal(values(remainders))
)

-----------------------------------------------------------------------------

ethRootSubStrat = (e,p,q,k,I,R) -> (
    n := numgens R;
    Rvars := R_*;
    Y := local Y;
    S := k(monoid[(Rvars | toList(Y_1..Y_n)), MonomialOrder=>ProductOrder{n,n},MonomialSize=>64]);
    --Produces a polynomial ring with twice as many variables as R.  The peculiar notation in the previous two lines
    --is required to ensure that the variables of S are hidden from the user.  In particular, the variables in R_* are
    --still recognized as variables of R and not S, and the code will not break if the variables in R happen to be called
    --Y_i also.  
    Svars := S_*;
    J := ideal(apply(n,i->Svars#(n+i) - Svars#i^(p^e)));
    H := apply((substitute(I,S))_*, f -> f % J);
    --If we denote the variables in R as X_1 .. X_n, then this replaces each occurrence of X_i^(p^e) in the polynomial f
    --with a Y_i.
    L := sum(H, f -> ideal((coefficients(f,Variables => Rvars))#1));
    --Peals off the coefficients of the basis polynomials for R over R^(p^e) as polynomials in the Y_i, and produces the
    --ideal generated by these coefficient polynomials.
    L = first entries mingens L;
    subRelations := apply(n,i->Svars#(n+i) => Svars#i);
    if q > p then subRelations = subRelations|{(gens k)#0 => getFieldGenRoot(e,p,q,k)};
    L = apply(L, g ->substitute(g,subRelations));
    --Pushes the ideal of coefficient polynomials down to R by substituting Y_i => X_i.
    --q := k#order;
    --Gets the size of the base field.
    substitute(ideal L, R)
)

--This tries to compute (f^a*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
--We should probably just use it.  It relies on the fact that (f^(ap+b))^{[1/p^2]} = (f^a(f^b)^{[1/p]})^{[1/p]}.
ethRootSafe = method(); 

ethRootSafe( ZZ, ZZ, RingElement, Ideal ) := ( e, a, f, I ) -> (
	R := ring I;
	p := char R;
	
	aRem := a%(p^e);
	aQuot := floor(a/p^e);
	
	expOfA := basePExp(p,aRem); --this gives digits of aRem base p as a list, left-endian first 
	zeros := apply(e-#expOfA,i->0);
	expOfA = expOfA|zeros;
	
	IN1 := I;
	
	if (e > 0) then (
		IN1 = IN1*ideal(f^(expOfA#0));
		IN1 = ethRoot( 1, IN1 );
		i := 1;
	
		while(i < #expOfA) do (
			IN1 = ethRoot( 1, IN1*ideal(f^(expOfA#i)) );
			i = i + 1;
		)
	);
	IN1*ideal(f^(aQuot))
)

ethRootSafe( ZZ, ZZ, RingElement ) := ( e, a, f ) -> 
    ethRootSafe( e, a, f, ideal( 1_(ring f) ) )

--This tries to compute (f1^a1*f2^a2*...fk^ak*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
ethRootSafeList = method();

ethRootSafeList( ZZ, List, List, Ideal ) := ( e, aList, elmtList, I ) -> (
    R := ring I;
    p := char R;
    
    aListRem := aList % p^e;
    aListQuot := aList // p^e;
    
    -- gives the basePexpansion of each element of aListRem
    -- expOfaList is thus a list of lists.
    expOfaList := apply(aListRem, z -> basePExp( p, z ) );
    
    -- this computes { ... f_i^b_i ... } where b_i = a_i % p
    aPowerList := apply(elmtList, expOfaList, (f, z) -> f^(z#0));
    
    IN1 := I*ideal(product(aPowerList));
    if (e > 0) then (
        IN1 = ethRoot( 1, IN1 );
        i := 1;
        while(i < e) do (
            aPowerList = apply(elmtList, expOfaList, (f, z) -> f^(z#i));
            IN1 = ethRoot( 1, IN1*ideal(product(aPowerList)) );
            i = i + 1;
        )
    );
    aPowerList = apply(elmtList, aListQuot, (f, z) -> f^z);
    IN1*ideal(product(aPowerList))
)

ethRootSafeList( ZZ, List, List ) := ( e, a, F ) ->
    ethRootSafeList( e, a, F, ideal( 1_( ring( F#0 )  ) ) )
	

-- DS: what does fancyEthRoot do?
-- is it (I^m)^[1/p^e]? Then it's just a special case of ethRootSafeList, no?
fancyEthRoot = (e,m,I) ->
(
    if (e < 0) then (error "expected a nonnegative integer e"); 
	G:=first entries mingens I;
	k:=#G;
	P:=allPartitions(m,k); --partitions of m of size k
	                       --this function is used here and only here

	                       --it seems like the point of this function
	                       --is just to find the generators of I^m where 
	                       --each generator of I appears
	R:=ring(I);
	p:=char(R);
	answer:=ideal(0_R);
	apply(P, u-> -- take one of your partitions and call it P
	{
		a:=ideal(1_R);
		U:=apply(u, v->baseP1(v,p,e)); -- so we're replacing each element of  P with its
		                               -- base P expansion (but stopping at p^e)
		for i from 0 to e do
		{
			j:=e-i;
			g:=1_R;
			for l from 0 to k-1 do g=g*(G#l)^((U#l)#j);  -- take each generator of I and raise it to 
			                                             -- the the coefficient of p^j of the base-p expansion
			                                             -- of appropriate element of your partition, a la the method
			                                             -- of EthRootSafe. 
			                                             -- Then multiply these together. 
			a=ideal(g)*a;
			if (i<e) then a=ethRoot( 1, a ); -- so we're iteratively doing answer = (prevanswer*f^a_1*... *f_k^a_k)^[1/p]
		};
		answer=answer+a; --adding up the above results for each partition
	});
	answer
)


--this is a working name for a function which should generalize both ethRootSafe and fancyEthRoot
--it will take 
smartEthRoot = method();

--smartEthRoot(ZZ, List, List) := (e, exponentList, idealList) -> (
----the idealList, can take a list of ideals, a list of elements, or a mix of lists of ideals or elements
--    -- make a list of the number of generators of each ideal in idealList
--    minGensList = apply(idealList, jj -> (if (class jj === Ideal) then #(first entries mingens (jj)) else 1 ));
--
--    -- see what's the biggest power of p smaller than a/m where 'a' is in exponentList
--    -- and m is in minGensList. I.e. we want to find the largest e such that a >= mp^e
--    minGensLog = apply(minGensList, exponentList, (mm, aa) -> (
--        n = floorLog(p, aa//mm);
--        if (n > e) then e else n
--    ));
--
--    tripleList = sort apply(minGensLog, idealList, exponentList, (a,b,c) -> {a,b,c});
--
--
--    R := ring(idealList#0);
--    answer :=  ideal(1_R);
--    p := char(R);
--
--    for i from 0 to length(idealList) - 1 do (
--        answer = answer*(idealList#j)^(exponentList#j - p^minGensLog#j)
--    );
--    
--    j := 0;
--    for i from 0 to e do (
--        if i == tripleList#j#0 then (
--            answer = answer*
--            j = j+1;
--        );
--    );
--
--);

--this function is the same as the above, it just explicitly adds J to the end of the ideal list and 1 to the end of the exponent list
smartEthRoot( ZZ, List, List, Ideal) := (e, exponentList, idealList, J) ->
   smartEthRoot(e, append(exponentList, 1), append(idealList, J));

smartEthRoot( ZZ, ZZ, RingElement, Ideal ) := ( e, a, f, I ) -> 
    smartEthRoot(e, {a, 1}, {f, I});

smartEthRoot( ZZ, ZZ, RingElement ) := ( e, a, f ) -> smartEthRoot(e, {a}, {f}); 

smartEthRoot( ZZ, ZZ, Ideal ) := ( e, m, I ) -> fancyEthRoot( e, {m}, {I} );

-- this a slow but easy to write implementation of smartEthRoot
smartEthRoot= method();
smartEthRoot(ZZ, List, List) := (e, exponentList, idealList) -> (
    --include the following line to set a break point: 
    --error "break here";
    I := null;
    if e == 0 then (
        I = idealList#0^(exponentList#0);
        for j from 1 to length(idealList) - 1 do I = I*(idealList#j)^(exponentList#j);
        return I;
    );

    R := ring(idealList#0);
    p := char(R);
    minGensList = apply(idealList, jj -> (if (class jj === Ideal) then #(first entries mingens (jj)) else 1 ));

    -- find max n such that a - (n-1)p > m*p. This is the number of copies of $I$ we can
    -- move outside the pth root. 

    nsList = apply(exponentList, minGensList, (aa, mm) -> (
       max(0, floor(aa/p - mm + 1)) 
    ));
    I = R;
    for j from 0 to length(idealList) - 1 do I = I*(idealList#j)^(exponentList#j - nsList#j * p);
    I = ethRoot(1, I);
    smartEthRoot(e - 1, append(nsList, 1), append(idealList, I))
);


----------------------------------------------------------------
--************************************************************--
--Functions for computing test ideals, and related objects.   --
--************************************************************--
----------------------------------------------------------------


--Finds the smallest phi-stable ideal containing the given ideal Jk
--in a polynomial ring Sk
--Jk is the given ideal, ek is the power of Frobenius to use, hk is the function to multiply 
--trace by to give phi:  phi(_) = Tr^(ek)(hk._)
--This is based on ideas of Moty Katzman, and his star closure
ascendIdeal = (ek, hk, Jk) -> (
     Sk := ring Jk;
     pp := char Sk;
     IN := Jk;
     IP := ideal(0_Sk);
     --we want to make the largest ideal that is phi-stable, following Moty Katzman's idea
     --we do the following
     while (isSubset(IN, IP) == false) do(
     	  IP = IN;
--	  error "help";
	  IN = ethRoot(ek,ideal(hk)*IP)+IP
     );

     --trim the output
     trim IP
)

--Works like ascendIdeal but tries to minimize the exponents elements are taken to
ascendIdealSafe = ( ak, ek, hk, Jk) -> (
	Sk := ring Jk;
     pp := char Sk;
     IN := Jk;
     IP := ideal(0_Sk);
     --we want to make the largest ideal that is phi-stable, following Moty Katzman's idea
     --we do the following
     while (isSubset(IN, IP) == false) do(
     	  IP = IN;
--	  error "help";
	  	IN = ethRootSafe( ek, ak, hk, IP )+IP
     );

     --trim the output
     trim IP
)




--works just like ascendIdealSafe but also handles lists of hk to powers...
ascendIdealSafeList ={AscentCount=>false} >> o ->  (akList, ek, hkList, Jk) -> (
	Sk := ring Jk;
	pp := char Sk;
	IN := Jk;
	IP := ideal(0_Sk);
	i1 := 0;
	--we ascend the ideal as above
	while (isSubset(IN, IP) == false) do(
		i1 = i1 + 1; --
		IP = IN;
		IN = ethRootSafeList( ek, akList, hkList, IP ) + IP
	);
	
	--trim the output
	if (o.AscentCount == false) then 
		trim IP
	else (trim IP, i1)
)

--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK
-- minimalCompatible is a method which is implemented as:
-- (1) the finding of the smallest ideal J which satisfies uJ\subset J^{[p^e]} 
---    containg a given ideal for a given ring element u,
-- (2) the finding of the smallest submodule V of a free module which satisfies UV\subset V^{[p^e]} 
--     containg a given submodule for a given matrix U.
--minimalCompatible = method();
--minimalCompatible(Ideal,RingElement,ZZ) :=  (Jk, hk, ek) -> ascendIdeal (Jk, hk, ek)
--minimalCompatible(Ideal,RingElement,ZZ,ZZ) :=  (Jk, hk, ak, ek) -> ascendIdealSafe (Jk, hk, ak, ek)
--minimalCompatible(Matrix,Matrix,ZZ) := (A,U,e) -> Mstar (A,U,e)
minimalCompatible = method();
minimalCompatible(ZZ,RingElement,Ideal) :=  (ek, hk, Jk) -> ascendIdeal (ek, hk, Jk)
minimalCompatible(ZZ,ZZ,RingElement,Ideal) :=  ( ak, ek,hk,Jk) -> ascendIdealSafe (ak, ek,hk,Jk)
minimalCompatible(ZZ,Matrix,Matrix) := (e,A,U) -> Mstar (e,A,U)

--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK



-----------------------------------------------------------------------------
--- Extend the Frobenius p^e th roots and star operations to submodules of
--- free modules (over polynomial rings with *prime* coeeficient field)
--- This implements the methods described in 
--- Moty Katzman and Wenliang Zhang's paper
--- "Annihilators of Artinian modules compatible with a Frobenius map"
--- Journal of Symbolic computation, 2014

-----------------------------------------------------------------------------


getExponents=(f)->(
answer:={};
t:=terms(first first entries f);
apply(t, i->
{
	exps:=first exponents(i);
	c:=(coefficients(i))#1;
	c=first first entries c;
	answer=append(answer,(c,exps));
});
answer
)

mEthRootOfOneElement= (e,v) ->(
	local i; local j;
	local d;
	local w;
	local m;
	local answer;
	R:=ring(v); p:=char R; q:=p^e;
	F:=coefficientRing(R);
	n:=rank source vars(R);
	V:=ideal vars(R);
	vv:=first entries vars(R);
	T:=new MutableHashTable;
	alpha:=rank target matrix(v);
	B:={};
	for i from 1 to alpha do
	{
		vi:=v^{i-1};
---print("i=",i);
---print("vi=",vi);
		C:=getExponents(vi);
---print(C);
		apply(C, c->
		{
			lambda:=c#0;
			beta:=c#1;
			gamma:=apply(beta, j-> (j%q));
			B=append(B,gamma);
			key:=(i,gamma);
---print(beta, #beta,vv);
			data:=apply(1..(#beta), j-> vv_(j-1)^((beta#(j-1))//q));
			data=lambda*product(toList data);
---print(beta, key, data);
			if (T#?key) then
			{
				T#key=(T#key)+data;
			}
			else
			{
				T#key=data;
			};
		});
	};
	B=unique(B);
	TT:=new MutableHashTable;
	apply(B, b->
	{
		ww:={};
		for i from 1 to alpha do if T#?(i,b) then ww=append(ww,T#(i,b)) else ww=append(ww,0_R);
		ww=transpose matrix {ww};
		TT#b=ww;
	});
	KEYS:=keys(TT);
	answer=TT#(KEYS#0);
	for i from 1 to (#KEYS)-1 do answer=answer | TT#(KEYS#i);
	answer
)

mEthRoot = (e,A) ->(
	local i;
	local answer;
	answer1:=apply(1..(rank source A), i->mEthRootOfOneElement (e, A_{i-1}));
	if (#answer1==0) then 
	{
		answer=A;
	}	
	else
	{
		answer=answer1#0;
		apply(2..(#answer1), i->answer=answer | answer1#(i-1));
	};
	mingens( image answer )
)	








--- Mstar is the implementaion of the star closure operation desribed in 
--- M Katzman's "Parameter test ideals of Cohen Macaulay rings" 
--- Input:
---    ideals I and element u of the same polynomial ring R OVER A PRIME FIELD.
---    a positive integer e
---    a prime p which is the characteristic of R
--- Output:
---    the smallest ideal J of R containing I with the property that u^(1+p+...+p^(e-1)) J is in J^{[p^e]}
Mstar = (e,A,U) ->(
	local answer;
	R:=ring(A); p:=char R;
	if (A==0) then
	{
		answer=A;
	}
	else
	{
		f:=true;
		Ne:=sum toList(apply(0..(e-1), i->p^i));
		lastA:= A;
		while (f) do
		{
			f=false;
			A1:=mEthRoot(e, mingens image ((U^Ne)*lastA));
			A1=A1 | lastA;
			t1:=compress ((A1))%((lastA));
			if (t1!=0) then 
			{
				f=true;
				lastA=mingens image A1;
			};
		};
		answer=mingens (image A1);
	};
	use(R);
	answer
)


--- end of MK ---------------------------------------------------------------------------------------------------

