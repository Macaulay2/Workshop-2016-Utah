"carryTest",  
    "basePExp",    
    "digit", 	   
    "denom",   
    "divideFraction",
    "firstCarry", 
    "floorLog",
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
     	floorLog
     Headline
        Computes the floor of the log base b of x
     Usage
     	 floorLog(b,x)
     Inputs 
     		b:ZZ
		x:ZZ		
     Outputs
         :ZZ
     Description
	Text
	    This differs from floor(log_b(x)) in that it corrects problems due to rounding.
/// 

doc ///
     Key
     	multOrder
     	(multOrder, ZZ, ZZ)
     Headline
        Computes the multiplicative order of a modulo b
     Usage
     	 multOrder(a,b)
     Inputs 
     		a:ZZ
		b:ZZ		
     Outputs
         :ZZ
     Description
	Text
	    This computes the multiplicative order of a modulo b.  If a and b are not relatively prime, it returns an error.
///

 
 
