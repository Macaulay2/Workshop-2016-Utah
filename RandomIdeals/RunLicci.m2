needsPackage "RandomIdeal"
needsPackage "ResidualIntersections"
L = apply(100,i->randomShelling(5,1,5))
LI = apply(100,i->idealFromShelling(L_i))
Licci = apply(100,i->isLicci(LI_i))
