R = ZZ/7[x,y,z]
phi = rationalMapping map(R, R, {y*z, x*z, x*y})
ident = rationalMapping map(R, R)
phi^2 == ident