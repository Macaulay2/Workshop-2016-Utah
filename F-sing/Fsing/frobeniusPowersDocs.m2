--***********************************************
--***********************************************
--Documentation for frobeniusPowers.m2
--***********************************************
--***********************************************

doc ///
     Key
     	fastExp 
     Headline
        Computes powers of elements in rings of characteristic p>0 quickly.
     Usage
     	  fastExp(N,f) 
     Inputs
         N:ZZ
     	 f:RingElement
     Outputs
        :RingElement
     Description
	Text
	     In prime characteristic p > 0, raising a sum (a+b) to a power is more quickly done by simply computing a^p 
	     and b^p and adding them.  The basic strategy is to break up the exponent into its base p expansion, and then 
	     use the exponent rules.  For example, (x+y)^(3*p^2 + 5*p+2) = ((x+y)^3)^(p^2)*((x+y)^5)^p*(x+y)^2.
///

doc ///
     Key
     	 frobenius
     Headline
        The following raises an ideal or matrix (entry-wise) to the p^e-th power.
     Usage
     	  frobenius(e,I)
	  frobenius^e(I) 
	  frobenius(e,M)
	  frobenius^e(M)
     Inputs
         e:ZZ
     	 I:Ideal
	 M:Matrix
     Outputs
        :Ideal
	:Matrix
     Description
	Text
	     frobenius(e,I) or frobenius^e(I) outputs I^[p^e] and frobenius(e,M) or frobenius^e(M) outputs a matrix whose entries are p^e-th powers of
	     the entries of M.
///

doc ///
     Key
     	frobeniusPower 
     Headline
        Computes the generalized Frobenius power of an ideal
     Usage
     	  frobeniusPower(t,I)
     Inputs
     	     	t:QQ
         	I:Ideal
     Outputs
        :Ideal
     Description
     	Text
	   frobeniusPower(t,I) outputs the generalized Frobenius power I^[t].
 ///
 
 