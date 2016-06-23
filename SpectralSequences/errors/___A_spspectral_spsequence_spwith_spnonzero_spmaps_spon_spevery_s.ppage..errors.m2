-- -*- M2-comint -*- {* hash: -1889126683 *}

i1 : S = ZZ/101[x,y];

i2 : I = ideal(x^2,x*y,y^2);

o2 : Ideal of S

i3 : R = S/I;

i4 : kR = coker vars R;

i5 : kS = coker vars S;

i6 : CS = res kS;

i7 : CR = res(kR,LengthLimit=>6);

i8 : CS' = CS**R;

i9 : E = prune spectralSequence (CS' ** filteredComplex CR);

i10 : use ZZ[t]

o10 = ZZ[t]

o10 : PolynomialRing

i11 : easyPresentation = (P,n,m) -> (
         transpose matrix apply(n,
             i-> apply(m,
                 j-> (rank (P_{i,j}))*t^(
                     if (L = unique flatten degrees P_{i,j})!= {} then first L else 0)
                 )
             ));

i12 : easyPresentation(E_infinity,6,3)
stdio:14:32:(3):[4]: error: can't promote number to ring
