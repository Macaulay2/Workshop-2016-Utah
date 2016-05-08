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
