doc ///
    Key
        findQGorGen
        (findQGorGen, ZZ, Ring)
        (findQGorGen, Ring)
    Headline
        finds an element representing the Frobenius trace map of a Q-Gorenstein ring
    Usage
        findQGorGen(e, R)
        findQGorGen(R)
    Inputs
        e: ZZ
        R: Ring
    Outputs
        :RingElement
    Description
        Text
            Suppose that $R$ is a ring such that $(p^e-1)K_R \sim 0$ (for example, if $R$ is $Q$-Gorenstein with index not divisible by $p$).  Then if we write $R = S/I$ where $S$ is a polynomial ring, we have that $I^{[p^e]} : I = u S + I^{[p^e]}$ for some $u \in S$.  By Fedder's criterion, this element $u$ represents the generator of the $R^{1/p^e}$-module $Hom(R^{1/p^e}, R)$.  For example if $I = (f)$ is principal, then $u = f^{p^e-1}$ works.
        Text
            This function produces the element $f$ described above.  If do not specify an integer e, it assumes e = 1.
        Example
            S = ZZ/3[x,y,z];
            f = x^2*y - z^2;
            I = ideal(f);
            R = S/I;
            u = findQGorGen(1, R)
            u%I^3 == f^2%I^3
        Text
            If Macaulay2 does not recognize that $I^{[p^e]} : I / I^{[p^e]}$ is principal, an error is thrown.  Note in the nongraded case, Macaulay2 is not guaranteed to do this simplification.
///

doc ///
    Key
        testElement
        (testElement, Ring)
    Headline
        finds a test element of a ring
    Usage
        testElement(R)
    Inputs
        R: Ring
    Outputs
        :RingElement
    Description
        Text
            Given $R = S/I$ where $S$ is a polynomial ring, this finds an element of $S$ that restricts to a test element of $R$.  This does this by finding a minor of the Jacobian whose determinant is not in any minimal prime of the defining ideal of $R$.  It looks at random minors until one is found instead of computing all of them.
        Example
            R = ZZ/5[x,y,z]/(x^3+y^3+z^3);
            testElement(R)
            testElement(R)
            testElement(R)
///

doc ///
    Key
        testIdeal
        (testIdeal, Ring)
    Headline
        computes the test element of a ring or pair
    Usage
        testIdeal(R)
    Inputs
        R: Ring
    Outputs
        :Ideal
    Description
        Text
            Given a normal (or G1 + S2) Q-Gorenstein ring R, this computes the test ideal.          
///
