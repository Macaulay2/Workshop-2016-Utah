R = ZZ/101[a..d]
I = monomialCurveIdeal(R,{1,3,4})
C = res I
C.dd
betti C

R = ZZ/101[vars(0..17)]
m1 = genericMatrix(R,a,3,3)
m2 = genericMatrix(R,j,3,3)
J = ideal(m1*m2-m2*m1)
betti J
C = res J
betti C


restart
R = ZZ/101[vars(0..7)]
I = ideal fromDual matrix{{random(3,R)}}

elapsedTime C1 = res I
betti C1

J = ideal I_*
elapsedTime C2 = res(J, FastNonminimal=>true)
betti C2
betti(C2, Minimize=>true)


load "prym12.m2" --sets ideal I
elapsedTime C = res(I, FastNonminimal=>true)
betti C
betti(C, Minimize=>true)
-- the following takes perhaps 2-3 hours or so.
-- elapsedTime res(ideal I_*) 

load "prym14.m2"
elapsedTime C = res(I, FastNonminimal=>true)
betti C
betti(C, Minimize=>true)

restart
load "prym16.m2"
gbTrace=1
elapsedTime C = res(I, FastNonminimal=>true)
betti C
betti(C, Minimize=>true)
exit

elapsedTime C = res(ideal I_*, FastNonminimal=>true, LengthLimit=>7)
betti C
betti(C, Minimize=>true)

elapsedTime C = res(ideal I_*, FastNonminimal=>true, LengthLimit=>6)
betti C
betti(C, Minimize=>true)


run ///open "http://web.macaulay2.com"///
