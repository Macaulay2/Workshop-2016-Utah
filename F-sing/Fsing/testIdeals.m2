--**************************************************
--**************************************************
--This file houses functions related to test ideals.
--**************************************************
--**************************************************



--This internal function computes the test ideal of (R, f^(a/(p^e - 1))) 
--when R is a polynomial ring.

tauAOverPEMinus1Poly = (e,a,f) -> (
     R := ring f;
     p := char R;
     b := a % (p^e - 1);
     k := a // (pp^e - 1); 
     I := ethRootSafe(fm, ideal(f), a2, e);
     I = ascendIdealSafe(I, fm, a2, e1);
     I*ideal(f^k)
)

--------------------------------------------------------------------------------------------------------

--Computes the test ideal of (R, f^t), when R is a polynomial ring over a perfect field.

tauPoly = (t,f) -> (
     R := ring f; 
     p := char R;
     L := divideFraction(t,p); --L={a,b,c}, where t = a/(p^b*(p^c-1))
     local I;
     --first we compute tau(R, f^{a/(p^c-1)})
     if (L1#2 != 0) then 
     (
         I1 = tauAOverPEMinus1Poly(L#2,L#0,f);
         if (L#1 != 0) then I = ethRoot(I, L#1)
     )
     else 
     (
         if (L#1 != 0) then I = ethRootSafe(f, ideal(sub(1, R)), L#0, L#1)
     	    else I = ideal(f^(L#0))
     );
     I
)

--------------------------------------------------------------------------------------------------------

--Given a quotient R of a polynomial ring modulo an ideal, outputs an ideal containing a nonzero 
--test element.   It views it as an element of the ambient ring of R.
--(NOTE:  One could make this faster by not computing the entire Jacobian / singular locus
--instead, if we just find one element of the Jacobian not in I, then that would also work
--and perhaps be substantially faster.)

findTestElementAmbient = R -> (
     J := ideal singularLocus(R);
     if (isSubset(J, ideal R) == true) then 
          error "findTestElementAmbient: No test elements found, is the ring non-reduced?";
     J      
)

--------------------------------------------------------------------------------------------------------

--Outputs the test ideal of a (Q-)Gorenstein ring (with no pair or exponent), where 
--the index divides (p^e - 1).  It outputs the appropriate stable/fixed ideal inside 
--the ambient ring

tauQGorAmb = (e, R) -> (
     J := findTestElementAmbient(R);
     h := findQGorGen(e,R);	
     sub(ascendIdeal(e,h,J),R)  
)

--------------------------------------------------------------------------------------------------------

--Computes the test ideal of a Gorenstein ring

tauGorAmb = R -> (tauQGorAmb(1,R))

--------------------------------------------------------------------------------------------------------

--This internal function computes the test ideals of a pair (R, f^(a/p^e-1)) where
--R = S/I. Other inuputs are J, a nonzero ideal contained in the test ideal, h, the 
--multiple used to form phi of the ambient ring, and ek, the power associated with h.

tauAOverPEMinus1QGorAmb = (S, Jk, hk, ek, f, a, e1) -> (
     p := char S;
     et := lcm(ek, e1);
     
     ak1 := numerator ((p^et - 1)/(p^ek - 1)); --an exponent for hk
     a3 := numerator (a*(p^et - 1)/(p^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
     a2 := a3 % (p^et - 1);
     k2 := a3 // (p^et - 1);    --save value to use Skoda tau(f^(1+k)) = f*tau(f^k)     
     Jl := ascendIdealSafe(Jk, hk, 1, ek);
                      
        --          assert false;                             
     Iasc := ascendIdealSafeList(Jk*ideal(f)^(ceiling(a3/(p^et - 1))), (f, hk), (a2, numerator ((p^et - 1)/(p^ek - 1))), et);
--     assert false;
     Iasc*ideal(fastExp(k2,f))
)

--------------------------------------------------------------------------------------------------------

--Computes the test ideal of (Rk, fk^t1), where Rk is a (Q-)Gorenstein ring, where index of the
--canonical module divides (p^ek - 1)
tauQGor := (Rk, ek, fk, t1) -> (
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
     k2 := a3 // d1; --we will use Skoda tau(f^(k+t)) = f^k*tau(f^t)			  
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
