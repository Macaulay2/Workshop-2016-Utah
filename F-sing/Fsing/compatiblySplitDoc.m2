--this doc file is done

doc ///
     Key
     	findAllCompatibleIdeals 
     Headline
        Finds all ideals compatible with a given Frobenius near-splitting
     Usage
     	 findAllCompatibleIdeals (u)
     Inputs
	u:RingElement
     Outputs
         :List
     Description
	Text
	     Given an element u of the polynomial ring R OVER A PRIME FIELD, this function returns a list of all prime ideals P such that: 
	     (a) u P \subseteq P^{[p]}, and 
	     (b) u is not in P^[p]. 
	     Condition (b) is equivalent to the non-vanishing of the corresponding Frobenius map on annihilator of P on the injective hull of the residue field of R. This is an implementation of the algorithm described in Moty Katzman and Karl Schwede's paper "An algorithm for computing compatibly Frobenius split subvarieties" J. Symbolic Comput. 47 (2012), no. 8, 996-1008. 


///


