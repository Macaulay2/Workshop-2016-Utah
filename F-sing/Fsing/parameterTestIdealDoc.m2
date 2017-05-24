doc ///
    Key
        paraTestModuleAmbient 
    Headline
        (Ask Karl for various versions)
    Usage
        paraTestModuleAmbient(u)
    Inputs
        u:RingElement
    Outputs
        :List
    Description
        Text
            Ask Karl
///

doc ///
    Key
        canonicalIdeal 
    Headline
        Given a ring, produces an ideal isomorphic to the canonical module.
    Usage
        canonicalIdeal(R)
    Inputs
        R:Ring  
    Outputs
        :Ideal
    Description
        Text
            Given a ring $R$, typically a domain, this produces an ideal isomorphic to the canonical module of $R$.  This will not always produce the same ideal, especially in a non-domain.  
        Example
            S = QQ[x,y,u,v];
            T = QQ[a,b];
            f = map(T, S, {a^3, a^2*b, a*b^2, b^3});
            R = S/(ker f);
            canonicalIdeal(R)
        Text
            Here's an example in a non-domain.
        Example
            R = ZZ/13[x,y,z]/ideal(x*y, x*z, y*z);
            canonicalIdeal(R)
            canonicalIdeal(R)                    
///

doc ///
    Key
        findusOfIdeal
    Headline
        Finds the u, which in a polynomail ring, determines the Frobenius trace on canonical module of a quotient of that ring.
    Usage
        findusOfIdeal(canIdeal, defIdeal)
    Inputs
        canIdeal:Ideal
        defIdeal:Ideal
    Outputs
        :RingElement
    Description
        Text
            Given $R = S/I$, where $S$ is a polynomial ring, there is a map from the canonical module of $R$ back to itself, dual to the Frobenius on $R$.  This map comes from a $p$ inverse linear map on $S$, restricted appropriately.  But every $p$ inverse linear map on $S$ is a premultiple of the Grothendieck dual by some element $u$.  This function finds the $u$, or at least finds some elements, some linear combination of them is the actual $u$.
///

doc ///
    Key
        testModule
    Headline
        Finds the parameter test module of a reduced ring.
    Usage
        testModule(R)
    Inputs
        R:Ring
    Outputs
        :Sequence
    Description
        Text
            Computes the parameter test module (as a submodule of the canonical module).  The function returns three values, the parameter test submodule, the canonical module of which it is a subset, and the element u (or us) used to compute this ideal via the method findusOfIdeal.
///

