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
        (testModule, Ring)
        (testModule, Ring, Ideal)
        (testModule, QQ, RingElement)
        (testModule, ZZ, RingElement)
        (testModule, QQ, RingElement, Ideal, List)
        (testModule, List, List)
        (testModule, List, List, Ideal, List)
        [testModule, EthRootStrategy]
    Headline
        Finds the parameter test module of a reduced ring.
    Usage
        testModule(R)
        testModule(R, canIdeal)
        testModule(tt, ff)
        testModule(tt, ff, canIdeal, u1)
    Inputs
        R:Ring
        canIdeal:Ideal
        tt:QQ
        tt:ZZ
        u1:List
    Outputs
        :Sequence
    Description
        Text
            Computes the parameter test module (as a submodule of the canonical module).  The function returns three values, the parameter test submodule, the canonical module of which it is a subset, and the element u (or us) used to compute this ideal via the method findusOfIdeal.  
        Example
            R = ZZ/7[x,y,z]/ideal(x^3+y^3+z^3);
            testModule(R)
        Text
            Note this is a Gorenstein ring and so the ambient canonical module is the unit ideal.
        Example
            S = ZZ/3[x,y,u,v];
            T = ZZ/3[a,b];
            f = map(T, S, {a^3, a^2*b, a*b^2, b^3});
            R = S/(ker f);
            testModule(R)
        Text
            Note that the output in this case has the parameter test module equal to the canonical module, as it should be.  Let's a non-Gorenstein example which is not F-rational.
        Example
            R = ZZ/5[x,y,z]/ideal(y*z, x*z, x*y);
            paraTestMod = testModule(R)
            (paraTestMod#0) : (paraTestMod#1)
        Text
            This function can also be used to compute the parameter test module of a pair $(R, f^t)$.
        Example
            R = ZZ/7[x,y];
            f = y^2 - x^3;
            testModule(5/6, f)
            testModule(5/7, f)
        Text
            This can also be used to compute $(R, f^s g^t)$.
        Example
            R = ZZ/7[x,y];
            f = y^2 - x^3;
            g = x^2 - y^3;
            testModule({1/2, 1/2}, {f, g})        
        Text
            Finally, sometimes you would like to specify the ambient canonical module (and choice of u) across multiple calls of testModule.  Those are what the $canIdeal$ or $u1$ can be used to specify.
///
