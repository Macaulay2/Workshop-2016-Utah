This was for running certain tests from the paper.

restart
loadPackage "RationalMaps"; loadPackage "Cremona"
n = 5; d = 6;
R = ZZ/101[x_0..x_n];
L = {x_0^d, x_1*x_0^(d-1)} | toList(apply(2..n, i -> (x_i*x_0^(d-1) + x_(i-1)^d)));
psi = map(R, R, L);
time inv = inverseOfMap(psi, AssumeDominant=>true, CheckBirational=>false, Verbosity=>0);
