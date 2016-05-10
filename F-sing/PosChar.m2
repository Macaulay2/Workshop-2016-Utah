newPackage( "PosChar",
Version => "0.2a", 
Date => "May 30th, 2015", 
Authors => {
     {Name => "Erin Bela",
     Email=> "ebela@nd.edu"
     },
     {Name => "David J. Bruce",
     Email => "djbruce@math.wisc.edu",
     HomePage => "http://www.math.wisc.edu/~djbruce/"
     },
     {Name => "Daniel Hernandez",
     Email => "dhernan@math.utah.edu",
     HomePage => "http://math.utah.edu/~dhernan/"
     },
     {Name => "Zhibek Kadyrsizova",
     Email => "zhikadyr@umich.edu"
     },
     {Name => "Mordechai Katzman",
     Email=> "m.katzman@sheffield.ac.uk",
     HomePage=> "http://www.katzman.staff.shef.ac.uk/"
     },
     {Name => "Sara Malec",
     Email=> "smalec@gsu.edu"
     },
     {Name => "Karl Schwede",
     Email => "schwede@math.psu.edu",
     HomePage => "http://math.utah.edu/~schwede/"
     },
     {Name => "Pedro Teixeira",
     Email => "pteixeir@knox.edu",
     HomePage => "http://www.knox.edu/academics/faculty/teixeira-pedro.html"
     },
     {Name=> "Emily Witt",
     Email=> "ewitt@umn.edu",
     HomePage => "http://math.umn.edu/~ewitt/"
     }
},
Headline => "A package for calculations of singularities in positive characteristic", 
DebuggingMode => true, 
Reload => true 
)
export{
    "aPower",
    "ascendIdeal", -- EthRoots.m2
    "ascendIdealSafe", -- EthRoots.m2
    "ascendIdealSafeList", -- EthRoots.m2
    "AscentCount", -- EthRoots.m2
    "basePExp",
    "basePExpMaxE",
    "BinomialCheck", -- FThresholds.m2
    "binomialFPT", -- FThresholds.m2
    "canonicalIdeal",
    "canVector", -- FThresholds.m2
    "carryTest",
    "digit", 	 
    "denom",
    "DiagonalCheck", -- FThresholds.m2
    "diagonalFPT", -- FThresholds.m2
    "divideFraction",
    "estFPT", -- FThresholds.m2
    "ethRoot", -- EthRoots.m2
    "ethRootSafe", -- EthRoots.m2
    "ethRootSafeList", -- EthRoots.m2   
    "factorList", -- FThresholds.m2
    "fancyEthRoot", -- EthRoots.m2 
    "fastExp",  --frobeniousPowers.m2
    "findCPBelow", -- FThresholds.m2
    "findGeneratingMorphisms",     --MK
    "findHSLloci",                 --MK
    "findTestElementAmbient", -- testIdeals
    "FinalCheck", -- FThresholds
    "findAllCompatibleIdeals", 	--- MK
    "findQGorGen", -- testIdeals
    "finduOfIdeal",
    "firstCarry", 
    "FPTApproxList", -- FThresholds.m2    
    "FPT2VarHomog", -- FThresholds.m2    
    "FPT2VarHomogInternal", -- FThresholds.m2
    "fracPart",
    "frobenius", --Other
    "frobeniusPower",  --frobeniousPowers.m2 
    "fSig", -- Other
    "FTApproxList",  -- FThresholds.m2
    "FTHatApproxList",  -- FThresholds.m2
    "FullMap",--specifies whether the full data should be returned
    "getNumAndDenom",
    "genFrobeniusPower",   --frobeniousPowers.m2 
    "guessFPT",  -- FThresholds.m2
    "HSL", --Other
    "imageOfRelativeCanonical", -- Other
    "imageOfTrace", --doesn't work! -- Other
    "isBinomial",  -- FThresholds.m2
    "isCP",  -- FThresholds.m2
    "isDiagonal",  -- FThresholds.m2
    "isFJumpingNumberPoly",  -- FThresholds.m2
    "isFPTPoly",  -- FThresholds.m2
    "isFPure", -- Other
    "isFRegularPoly", -- Other
    "isFRegularQGor", -- Other
    "isInLowerRegion",  -- FThresholds.m2
    "isInUpperRegion",  -- FThresholds.m2
--    "isJToAInIToPe",
    "isSharplyFPurePoly", -- Other
    "isMapSplit", -- Other
    "MaxExp",  -- FThresholds.m2
    "maxIdeal", -- BasicFunctions
    "minimalCompatible", -- EthRoots.m2 
---    "Mstar",			--- MK
    "multOrder",
    "MultiThread", --FThresholds
    "nonFInjectiveLocus",   --MK
    "Nontrivial", --SpecialFThresholds
    "nu", -- FThresholds.m2
    "nuAlt", -- FThresholds.m2
    "NuCheck", -- FThresholds.m2
    "nuHat", -- FThresholds.m2
    "nuHatList", -- FThresholds.m2
    "nuList", -- FThresholds.m2
    "nuListAlt", -- FThresholds.m2
    "nuListAlt1", -- FThresholds.m2
    "num",
    "Origin", --FThresholds
    "OutputRange", -- FThresholds
    "paraTestModule",
    "paraTestModuleAmbient",
    "PrintCP",  -- FThresholds.m2
    "setFTData",  -- FThresholds.m2
    "sigmaAOverPEMinus1Poly", -- Other
    "sigmaQGorAmb", --needs optimization -- Other
    "sigmaAOverPEMinus1QGor",      --needs optimization --Other
    "splittingField",  -- FThresholds.m2
    "tauPoly", -- testIdeals
    "tauNonPrincipalAOverPEPoly", -- testIdeals
    "tauAOverPEMinus1Poly", --testIdeals
    "tauGor",--needs optimization -- testIdeals
    "tauGorAmb",--needs optimization --testIdeals
    "tauQGor",--needs optimization -- testIdeals
    "tauQGorAmb",--needs optimization -- testIdeals
    "taxicabNorm",  -- FThresholds.m2
    "truncation",
    "truncationBaseP"
}
--This file has "finished" functions from the Macaulay2 workshop at Wake 
--Forest in August 2012.  Sara Malec, Karl Schwede and Emily Witt contributed
--to it.  Some functions, are based on code written by Eric Canton and Moty
-- Katzman
--
--UPDATE January 2014 at Macaulay2 workshop at MSRI:  Daniel Hernandez, Moty 
--Katzman, Karl Schwede, Pedro Teixeira, Emily Witt added more functionality.
--
--UPDATE May 2015 at Macaulay2 workshop atBoise State:  Erin Bela, DJ Bruce,
-- Daniel Hernandez, Zhibek Kadyrsizova, and Emily Witt improved/fixed/added 
--functionality.



--****************************************************--
--*****************Documentation**********************--
--****************************************************--

beginDocumentation()

doc ///
   Key
      PosChar 
   Headline
      A package for calculations in positive characteristic 
   Description
      Text    
         This will do a lot of cool stuff someday. 
///

doc ///
     Key
     	aPower
     Headline
        Finds the largest power of p dividing x.
     Usage
     	 aPower(x,p)
     Inputs 
		x:ZZ
		p:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the largest exponent e such that p^e divides x.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- START: Transfered to EthRootDoc
doc ///
     Key
     	ascendIdeal
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdeal(J, h, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.
///

doc ///
     Key
     	ascendIdealSafe
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdealSafe(J, h, a, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	a:ZZ
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h^a.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.  It differs from ascendIdeal in that it minimizes the exponents that h is raised to, this can make it faster or slower depending on the circumstances.
///
--- END: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	basePExp 
     Headline
        Base p Expansion of an integer N
     Usage
     	  basePExp(N,p) 
     Inputs
         N:ZZ
	 p:ZZ
     Outputs
        :List
     Description
     	Text
	     Given an integer N and a prime p, outputs the digits of the base p expansion of N in base p.
///

doc ///
     Key
     	basePExpMaxE
     Headline
        Computes the base-p expansion of N from digits zero to e-1.
     Usage
     	 basePExpMaxE(N,p,e)
     Inputs
     	 N:ZZ
	 	 p:ZZ
	 	 e:ZZ
     Outputs
         :List
     Description
	Text
	     This computes the base p expansion of N, from digits 0 to e-1.  The digits are given in a list, and come with leading zeros.  If fewer than e digits are required, the list is padded with zeros.  If more digits are required, the final digit lists them.  Little endian is first.  For example, if p=5 and N = 16, the basePExpMaxE(16,5,4) will return {1,3,0,0} (1 one, 3 fives, 0 twentyfives, 0 onehundred twentyfives).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	binomialFPT
     Headline
        Computes the F-pure threshold of a binomial polynomial.
     Usage
     	 binomialFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a binomial in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	carryTest
     Headline
        Finds the number of digits we must check to see whether x and y add without carrying.
     Usage
     	 carryTest(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  This function finds the number of digit places we must check to see if x and y add without carrying.
///

doc ///
     Key
     	denom
     	(denom,ZZ)
     	(denom,QQ)
     Headline
        Returns the denominator of a rational number.
     Usage
     	 denom(x)
     	 denom(y)
     Inputs 
		x:QQ
		y:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the denominator of a rational number or integer (in the latter case it returns 1).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	diagonalFPT
     Headline
        Computes the F-pure threshold of a diagonal polynomial.
     Usage
     	 diagonalFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a diagonal hypersurface in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	digit
	(digit,ZZ,QQ,ZZ)
	(digit,ZZ,List,ZZ)
     Headline
        Gives the e-th digit of the base p expansion 
     Usage
     	 d=digit(e,x,p), D=digit(e,X,p)
     Inputs
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   consisting of rational numbers
     Outputs
        d:ZZ
	    which is the e-th digit of the non-terminating base p expansion of x
	D:List
	    which contains the e-th digits of the entries of the list X
     Description
	Text
	     Gives the e-th digit, to the right of the decimal point, of the non-terminating base p expansion of x in [0,1]; threads over lists of rational numbers. 
///

doc ///
     Key
     	divideFraction
     Headline
        Converts a rational number into something of the form (a/(p^b p^(c-1)).
     Usage
     	 divideFraction(t, p)
     Inputs 
		t:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Given a rational number t and prime p, this function finds a list of integers {a,b,c} such that t= (a/(p^b p^(c-1)).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 estFPT
     Headline
         Atempts to compute the F-pure threshold, where e is the max depth to search in.  
     Usage
     	  estFPT(f,e,finalCheck=>V,Verbose=>W)
     Inputs
     	 f:RingElement
         e:ZZ
	 V:Boolean
	 W:Boolean
     Outputs
        L:List
	Q:QQ
     Description
     	  Text 
	      This tries to find an exact value for the fpt.  If it can, it returns that value.  Otherwise it should return a range of possible values (eventually).  It first checks to see if the ring is binonmial or diagonal.  In either case it uses methods of D. Hernandez.  Next it tries to estimate the range of the FPT using nu's.  Finally, it tries to use this to deduce the actual FPT via taking advantage of convexity of the F-signature function and a secant line argument.  finalCheck is a Boolean with default value True that determines whether the last isFRegularPoly is run (it is possibly very slow).  If FinalCheck is false, then a last time consuming check won't be tried.  If it is true, it will be.  Verbose set to true displays verbose output.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- START: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doc ///
     Key
     	 ethRoot
     Headline
        Computes $I^{[1/p^e]}$ in a polynomial ring over a perfect field
     Usage
     	  ethRoot(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     In a polynomial ring k[x1, ..., xn], I^{[1/p^e]} is the smallest ideal J such that J^{[p^e]} = FrobeniusPower(J,e) \supseteq I.  This function computes it.
///

doc ///
     Key
     	ethRootSafe
     Headline
        Computes (f^a*I)^{[1/p^e]} in such a way that we don not blow exponent buffers.
     Usage
     	 ethRootSafe(f, I, a, e)
     Inputs
     	 f:RingElement
	 I:Ideal
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the 1/p^e-th root of (f^a*I).  It does it while trying to minimize the power that f gets raised to (in case a is a large number).  This can either be faster or slower than ethRoot.
///
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transfered to EthRootDoc
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	fastExp 
     Headline
        Computes powers of elements in rings of characteristic p>0 quickly.
     Usage
     	  fastExp(f,N) 
     Inputs
     	 f:RingElement
         N:ZZ
     Outputs
        :RingElement
     Description
	Text
	     In prime characteristic p > 0, raising a sum (a+b) to a power is more quickly done by simply computing a^p and b^p and adding them.  The basic strategy is to break up the exponent into its base p expansion, and then use the exponent rules.  For example, (x+y)^(3*p^2 + 5*p+2) = ((x+y)^3)^(p^2)*((x+y)^5)^p*(x+y)^2.
///

doc ///
     Key
     	findQGorGen
     Headline
        If R = S/I where S is a polynomial ring, returns the ring element with I^{[p^e]} : I = (f) + I^{[p^e]}.
     Usage
     	 findQGorGen(R, e)
     Inputs
     	 R:Ring
     Outputs
         :RingElement
     Description
	Text
	     If R is Q-Gorenstein with index not divisible by p, then I^{[p^e]} : I = (f) + I^{[p^e]}.  For some e.  This function tries to find the f.  If the argument e is left out then e is assumed to be 1.
///

doc ///
     Key
     	firstCarry
     Headline
        Finds the first spot where (the eth digit of x) + (the eth digit of y) >= p.
     Usage
     	 firstCarry(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  Finds the first place where (the eth digit of x) + (the eth digit of y) >= p, in other words where the numbers add with carrying.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 FPTApproxList
	 (FPTApproxList,Ideal,ZZ)
	 (FPTApproxList,RingElement,ZZ)
     Headline
        Gives a list of nu_I(p^d)/p^d for d=1,...,e.
     Usage
     	  FPTApproxList(I,e)
	  FPTApproxList(f,e) 
     Inputs
     	 I:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I(p^d)/p^d for d = 1, ..., e.  The {nu_I(p^d)/p^d} converge to the F-pure threshold.	     
///

doc ///
     Key
     	FPT2VarHomog
	(FPT2VarHomog,RingElement)
	(FPT2VarHomog,List,List)
     Headline
        F-pure threshold of a form in two variables
     Usage
     	  fpt=FPT2VarHomog(G), fpt=FPT2VarHomog(factors,multiplicities)
     Inputs 
	factors:List
	    which contains the linear factors of a form G in two variables 
	multiplicities:List
	    which contains the multiplicities of those linear factors in G
	G:RingElement
	    a form in two variables
     Outputs
        fpt:QQ
     Description
	Text
	    FPT2VarHomog computes the F-pure threshold of a homogeneous polynomial G
	    	in two variables. 
	    The polynomial G can be entered directly, or if the user knows a factorization
	    	G=L1^(a1)...Ln^(an) into linear forms, that can be used for improved 
		performance: FPT2VarHomog({L1,...,Ln},{a1,...,an}).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 frobeniusPower
     Headline
        The following raises an ideal to the $p^e$th power.
     Usage
     	  frobeniusPower(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     If I = ideal(x1, ..., xn), then frobeniusPower(I,e) outputs ideal(x1^(p^e), ..., xn^(p^e)) where p is the characteristic of the ring.
///

doc ///

     Key
     	 fSig
     Headline
        Computes the F-signature for a specific value $a/p^e$.
     Usage
     	  fSig(f,a,e)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
     Outputs
        :QQ
     Description
	Text
	     This computes the F-signature $s(R, f^{a/p^e})$ if R is a polynomial ring over a perfect field.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 FTApproxList
	 (FTApproxList,Ideal,Ideal, ZZ)
	 (FTApproxList,RingElement,Ideal,ZZ)
     Headline
        Gives a list of nu_I^J(p^d)/p^d for d=1,...,e.
     Usage
     	  FTApproxList(I,J,e)
	  FTApproxList(f,J,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I^J(p^d)/p^d for d = 1, ..., e.  The {nu_I^J(p^d)/p^d} converge to the F-threshold.	     
///

doc ///
     Key
     	 FTHatApproxList
	 (FTHatApproxList,Ideal,Ideal, ZZ)
	 (FTHatApproxList,RingElement,Ideal,ZZ)
     Headline
        Gives a list of nuHat_I^J(p^d)/p^d for d=1,...,e.
     Usage
     	  FTHatApproxList(I,J,e)
	  FTHatApproxList(f,J,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nuHat_I^J(p^d)/p^d for d = 1, ..., e.  The {nuHat_I^J(p^d)/p^d} converge to the FHat-threshold.	     
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	genFrobeniusPower 
     Headline
        Computes the generalized Frobenius power of an ideal
     Usage
     	  frobeniusPower(I1,e1)  
     Inputs
         	I1:Ideal
	 	e1:ZZ
     Outputs
        :Ideal
     Description
     	Text
	     Computes I^[N] for an ideal I and an integer N, where I^[N] is defined as follows. If N's base P-expansion is N=n_0+n_1P+...+n_eP^e then I^[N]=I^(n_0)*(I^(n_1))^[P]*...*(I^(n_e))^[P^e]. When P is prime I^[P^e] is the usual Frobenius power.
 ///
 
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	guessFPT 
     Headline
        Tries to guess the FPT in a really naive way (this should be improved).
     Usage
     	  guessFPT(f,e,d) 
     Inputs
     	 f:RingElement
         e:ZZ
	 d:ZZ
     Outputs
        :List
     Description
	Text
	     This tries to guess the FPT.  In particular, it computes the number nu such that nu/(p^e - 1) <= FPT < (nu+1)/p^e.  It then outputs a list of all rational numbers with denominators less than or equal to d, which lie in that range.  WARNING:  There are several improvements which should be made to this function to rule out many of the possibilies.
///

doc ///
     Key
     	isBinomial 
     Headline
        Checks whether a polynomial is binomial.
     Usage
     	 isBinomial(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a binomial, otherwise returns false.
///

doc ///
     Key
     	isDiagonal 
     Headline
        Checks whether a polynomial is diagonal.
     Usage
     	 isDiagonal(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a diagonal, otherwise returns false.  Recall f is called diagonal if it is of the form x_1^(a_1)+...+x_n^(a_n) up to renumbering of the variables.
///

doc ///
     Key
     	isFJumpingNumberPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFJumpingNumberPoly(f,t,Verbose=>V)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is an F-jumping number, otherwise it returns false.
///

doc ///
     Key
     	isFPTPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFPTPoly(f,t,Verbose=>V,Origin=>W)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
		W:Origin
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is the FPT, otherwise it returns false.  If Origin is true, it only checks it at ideal(vars ring f).
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 isFPure 
     Headline
         Tests for a given ideal I, if R/I is F-pure. 
     Usage
     	 isFPure(I)
     Inputs
     	 I:Ideal
     Outputs
         :Boolean
     Description
	Text 
	    In the case where I is a complete intersection,this function applies Fedder's Criterion.
	    Otherwise, checks if I^[p]:I is contained in m^[p]. 
 	   
///

doc ///
     Key
     	 isFRegularPoly
     Headline
        Determines if a pair $(R, f^t)$ is F-regular when R is a polynomial ring. 
     Usage
     	  isFRegularPoly
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Boolean
     Description
	Text
	     This computes the test ideal.  The ring is F-regular if the test ideal is the whole ring, in which case this function returns true.  Otherwise, this function returns false.

///

doc ///
     Key
     	isFRegularQGor
	 (isFRegularQGor, ZZ, RingElement, QQ)
	 (isFRegularQGor, ZZ, RingElement, QQ, Ideal)
	 (isFRegularQGor, Ring, ZZ)
	 (isFRegularQGor, Ring, ZZ, Ideal)
     Headline
        Checks whether a ring or a pair is Q-Gorenstein.
     Usage
     	 isFRegularQGor(e,f,t)
     	 isFRegularQGor(e,f,t,Q)
     	 isFRegularQGor(R,e)
     	 isFRegularQGor(R,e,Q)
     Inputs 
		R:Ring
	 f:RingElement
	 e:ZZ
	 t:QQ
	 Q:Ideal
     Outputs
         :Boolean
     Description
	Text
	     Checks whether R, or the pair (R, f^t),  is strongly F-regular at Q (respectively the origin).  It assumes the Q-Gorenstein index divides (p^e - 1).
///

doc ///
     Key
     	 isSharplyFPurePoly
     Headline
        Checks whether (R, f^(a/(p^e - 1))) is F-pure at the prime ideal m.
     Usage
     	 isSharplyFPurePoly(f,a,e,m)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
	 m:Ideal
     Outputs
         :Boolean
     Description
	Text
	     This checks whether (R, f^(a/(p^e-1))) is F-pure at the prime ideal m at least in the case that R is a polynomial ring.
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- BEGIN: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	 nu
	 (nu,Ideal,Ideal,ZZ)
	 (nu,Ideal,ZZ)
	 (nu, RingElement,Ideal, ZZ)
	 (nu, RingElement, ZZ)
     Headline
        Gives $\(nu_I)^J(p^e)$ or $\(nu_f)^J(p^e)$
     Usage
     	  nu(I,J,e)
	  nu(I,e)
	  nu(f,J,e)
	  nu(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :ZZ
     Description
	Text
	    Given an ideal I in a polynomial ring k[x1, ..., xn], this function outputs the maximal integer nu such that I^nu is not in ideal J^[p^e].  If the input is (Ideal,ZZ) then the function computes the maximal integer nu such that I^nu in not in (x_1, ...,x_n)^[p^e]. If a RingElement is passed, it computes nu of the principal ideal generated by this element.This is used frequently to compute the F-pure threshold.
///

doc ///
     Key
     	 nuList
	 (nuList, Ideal,Ideal,ZZ)
	 (nuList, Ideal, ZZ)
	 (nuList, RingElement, Ideal, ZZ)
	 (nuList, RingElement, ZZ)
     Headline
        Lists $\(nu_I)^J(p^d)$ for d = 1,...,e.
     Usage
     	  nuList(I,J,e)
	  nuList(I,e)
	  nuList(f,J,e)
	  nuList(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :List
     Description
	Text
	     Given an ideal I in a polynomial ring k[x1,...,xn], this function computes nu(I,d) for d = 1,...,e. If a RingElement is passed, it computes nu of the principal ideal generated by this element for d=1,...,e
///

doc ///
     Key
     	 nuHat
	 (nuHat,Ideal,Ideal,ZZ)
	 (nuHat,Ideal,ZZ)
	 (nuHat, RingElement,Ideal, ZZ)
	 (nuHat, RingElement, ZZ)
     Headline
        Gives $\hat(nu_I)^J(p^e)$ or $\hat(nu_f)^J(p^e)$
     Usage
     	  nuHat(I,J,e)
	  nuHat(I,e)
	  nuHat(f,J,e)
	  nuHat(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :ZZ
     Description
	Text
	    Given an ideal I in a polynomial ring k[x1, ..., xn], this function outputs the maximal integer nu such that I^[nu] is not in ideal J^[p^e].  If the input is (Ideal,ZZ) then the function computes the maximal integer nu such that I^[nu] in not in (x_1, ...,x_n)^[p^e]. If a RingElement is passed, it computes nuHat of the principal ideal generated by this element.This is used frequently to compute the generalized F-pure threshold.
///

doc ///
     Key
     	 nuHatList
	 (nuHatList, Ideal,Ideal,ZZ)
	 (nuHatList, Ideal, ZZ)
	 (nuHatList, RingElement, Ideal, ZZ)
	 (nuHatList, RingElement, ZZ)
     Headline
        Lists $\hat(nu_I)^J(p^d)$ for d = 1,...,e.
     Usage
     	  nuHatList(I,J,e)
	  nuHatList(I,e)
	  nuHatList(f,J,e)
	  nuHatList(f,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :List
     Description
	Text
	     Given an ideal I in a polynomial ring k[x1,...,xn], this function computes nuHat(I,d) for d = 1,...,e. If a RingElement is passed, it computes nuHat of the principal ideal generated by this element for d=1,...,e
///

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- END: Transferred to FThresholds
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doc ///
     Key
     	sigmaAOverPEMinus1Poly
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}) when R is a polynomial ring.
     Usage
     	 sigmaAOverPEMinus1Poly (f, a, e, HSL=>W)
     Inputs 
		f:RingElement
	a:ZZ
	e:ZZ
	W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by f^a.  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	sigmaAOverPEMinus1QGor
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}).
     Usage
     	 sigmaAOverPEMinus1QGor(f, a, e, g,HSL=>W)
     Inputs 
		f:RingElement
		a:ZZ
		e:ZZ
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace of R by f^a (we assume that the Q-Gorenstein index of R divides p^g-1).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number of sigma(ring f) relative to f.
///

doc ///
     Key
     	sigmaQGor
     Headline
        Computes the non-sharply F-pure ideal of R, where R is Q-Gorenstein with index dividing (p^g-1).
     Usage
     	 sigmaQGor(R, g,HSL=>W)
     Inputs 
		R:Ring
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the  g-th Frobenius trace of R (we assume that g is the Q-Gorenstein index of R).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	tauAOverPEMinus1Poly
     Headline
        Computes the test ideal of f^(a/(p^e-1)) if f is in a polynomial ring.
     Usage
     	 tauAOverPEMinus1Poly(f, a, e)
     Inputs
     	 f:RingElement
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the test ideal tau(f^(a/(p^e-1)) ).  The basic idea first appeared in a paper of Mordechai Katzman.
///

doc ///
     Key
     	tauGor
     Headline
        Computes tau(R,f^t) for a Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauGor(R,f,t)
     Inputs
     	 R:Ring
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Gorenstein ring.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing the images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	tauGorAmb
     Headline
        Computes tau(R) for a Gorenstein ring.
     Usage
     	 tauGorAmb(R)
     Inputs
     	 R:Ring
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a quasi-Gorenstein ring R.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S. 
///

doc ///
     Key
     	 tauPoly
     Headline
        Computes the test ideal of $(R, f^t)$.
     Usage
     	  tauPoly(f,t) 
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Ideal
     Description
	Text
	     This computes the test ideal of (R, f^t) when R is a polynomial ring over a perfect field.  It is done as follows.  If t = a/(p^e - 1) then tau(R, f^t) is computed as a sum of (f^{lceil t rceil}*f^{lceil t(p^e-1) rceil})^{[1/p^e]} until the sum stabilizes.  For the more general case, we use the formula tau(R, f^t)^{[1/p^d]} = tau(R, f^{t/p^d}).
///

doc ///
     Key
     	tauQGorAmb
     Headline
        Computes tau(R) for a Q-Gorenstein ring with index not dividing p^e - 1.
     Usage
     	 tauQGorAmb(R,e)
     Inputs
     	 R:Ring
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a Q-Gorenstein ring R with index dividing p^e - 1.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S.  e is the index of the canonical divisor on R.
///

doc ///
     Key
     	tauQGor
     Headline
        Computes tau(R,f^t) for a Q-Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauQGor(R,e,f,t)
     Inputs
     	 R:Ring
	 e:ZZ
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Q-Gorenstein ring such that the index divides p^e -1.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	truncation
	(truncation,ZZ,QQ,ZZ)
	(truncation,ZZ,List,ZZ)
     Headline
        Truncations of base p expansions of rational numbers
     Usage
     	 t=truncation(e,x,p), T=truncation(e,X,p)
     Inputs 
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   which contains rational numbers
     Outputs
        t:QQ
	    which is the e-th truncation of the non-terminating base p expansion of x
	T:List
	    which contains the e-th truncations of the entries of the list X
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of a nonnegative rational number x, as a fraction; threads over lists of rational numbers.
///

doc ///
     Key
     	truncationBaseP
     Headline
        Gives the first e digits of the non-terminating base p expansion of x.
     Usage
     	 truncationBaseP(e,x,p)
     Inputs 
		e:ZZ
	x:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of x in [0,1], as a list.
///

end
