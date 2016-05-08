doc ///
     Key
     	generatingMorphism
     Headline
        Compute a generating morphism for local cohomology modules of a polynomial ring with support at I
     Usage
     	 generatingMorphism (i,I) 
     Inputs
	i:ZZ
	I:Ideal
     Outputs
         :Sequence
     Description
	Text
	     Compute the map Ext^i(R/I,R) -> Ext^i(R/I^[p],R) induced by the surjection R/I^[p]->R/I. The output is a pair (A,U), where A is a matrix whose cokernel is Ext^i(R/I,R)  and U is a square matrix defining the desired map from coker A to coker A^[p].
///

doc ///
     Key
     	findGeneratingMorphisms
     Headline
        Compute a generating morphisms for all local cohomology modules of a polynomial ring with support at I
     Usage
     	 findGeneratingMorphisms (I) 
     Inputs
	I:Ideal
     Outputs
         :List
     Description
	Text
	     Compute the map Ext^i(R/I,R) -> Ext^i(R/I^[p],R) induced by the surjection R/I^[p]->R/I. The ith entry in the list output is a pair (A,U), where A is a matrix whose cokernel is Ext^i(R/I,R)  and U is a square matrix defining the desired map from coker A to coker A^[p].
///

