"carryTest",  
    "basePExp",    
    "digit", 	   
    "denom",   
    "divideFraction",
    "firstCarry", 
    "floorlog",
    "fracPart", 
    "getCanVector",
    "getNumAndDenom", 
    "maxIdeal", 
    "multOrder",
    "num",
    "taxicabNorm",
    "truncatedBasePExp",
    
--***********************************************
--***********************************************
--Documentation for BasicFunctions.m2
--***********************************************
--***********************************************

doc ///
     Key
     	fastExp
     Headline
        Computes
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
     	 frobeniusPower
     Headline
        The following raises an ideal or matrix (entry-wise) to the p^e-th power.
     Usage
     	  frobeniusPower(e,I) 
	  frobeniusPower(e,M)
     Inputs
         e:ZZ
     	 I:Ideal
	 M:Matrix
     Outputs
        :Ideal
	:Matrix
     Description
	Text
	     frobeniusPower(e,I) outputs I^[p^e] and frobeniusPower(e,M) outputs a matrix whose entries are p^e-th powers of
	     the entries of M.
///

doc ///
     Key
     	genFrobeniusPower 
     Headline
        Computes the generalized Frobenius power of an ideal
     Usage
     	  genFrobeniusPower(t,I)
     Inputs
     	     	t:QQ
         	I:Ideal
     Outputs
        :Ideal
     Description
     	Text
	   genFrobeniusPower(t,I) outputs the generalized Frobenius power I^[t].
 ///
 
 
