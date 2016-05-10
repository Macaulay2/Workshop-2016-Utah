restart
load "Implicitization.m2"
software = "parametricM2"
--setDefault(Software=>PHCPACK)
NAGtrace 1
setRandomSeed 0
R = CC[s,t]
d = 66
F = flatten entries basis(d, R)
I = ideal(0_R)
numericalImageDegree(F, I, software) -- should be degree 6

F = delete(s^2*t^2, flatten entries basis(4, R))
numericalImageDegree(F, I, software) -- should be degree 4
S= QQ[s,t]
degree(ker map(S, QQ[x_1..x_366], basis(365,S)))
