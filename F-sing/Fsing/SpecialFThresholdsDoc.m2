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

