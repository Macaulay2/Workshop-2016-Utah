needsPackage"SchurRings"
n = 3
R = schurRing(r,n)
d = dim r_(a,b)

kk = ZZ/101
S = kk[x_(1,1)..x_(n,n)]
M = genericMatrix(S,n,n)

lis = for i from 0 to d^2-1 list
(
    A = random(kk^n,kk^n);
    B = random(kk^n,kk^n);
    N = A * M * B;
    (det(N_{0,1}^{0,1}))^b * (N_0_0)^(a-b)
    );
J = ideal lis;
I = ideal mingens J;
if (numgens I != d^2) then error"wrong number of generators"

--F = resolution(I,DegreeLimit => 5)
--F = resolution I
F = res(I, FastNonminimal => true)
betti(F, Minimize => true)
end

restart
a = 4
b = 2
time load"syz12.m2"
betti F
