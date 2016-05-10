TEST ///
R = ZZ/5[x,y,z,w]
I = ideal(x^27*y^10+3*z^28+4*x^2*y^15*z^35,x^17*w^30+2*x^10*y^10*z^35,x*z^50)
assert(ethRoot(1,I) == ideal(x^5*y^2+4*y^3*z^7,z^5,x^3*w^6,x^2*y^2*z^7,z^10))
assert(ethRoot(1,I,ethRootStrategy => MonomialBasis) == ideal(x^5*y^2+4*y^3*z^7,z^5,x^3*w^6,x^2*y^2*z^7,z^10))
assert(ethRoot(2,I) == ideal(x,z,w))
assert(ethRoot(2,I,ethRootStrategy => MonomialBasis) == ideal(x,z,w))
assert(ethRoot(3,I) == ideal(1_R))
assert(ethRoot(3,I,ethRootStrategy => MonomialBasis) == ideal(1_R))
///

TEST ///
R = GF(27)[x,y,z]
--The ambient ring of GF(27) is ZZ[a]/(a^3-a+1).
I = ideal(a^2*x^18+(a-1)*x^14*y^7*z^4 +x^2*y^10*z^10,(a^2-a)*x^5*y^9*z^8-y^21)
--a^(1/3) = a + 1
--a^(1/9) = a - 1
assert(ethRoot(1,I) == ideal(x^6,a*x^4*y^2*z+y^3*z^3,x*y^3*z^2,y^7))
assert(ethRoot(1,I,ethRootStrategy => MonomialBasis) == ideal(x^6,a*x^4*y^2*z+y^3*z^3,x*y^3*z^2,y^7))
assert(ethRoot(2,I) == ideal(x,y))    
assert(ethRoot(2,I,ethRootStrategy => MonomialBasis) == ideal(x,y))
assert(ethRoot(3,I) == ideal(1_R))
assert(ethRoot(3,I,ethRootStrategy => MonomialBasis) == ideal(1_R))    
///