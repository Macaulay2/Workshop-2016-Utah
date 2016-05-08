
binomPoly = method()
binomPoly (Ring, ZZ, ZZ) := RingElement => (R, a, b) -> (
	product(toList(1..b)/(i -> (first gens R + a + 1 - i)/i))
)

gotzmannBound = method()
gotzmannBound (RingElement, ZZ) := List => (P, s) -> (
	d := (degree P)#0;
	a := leadCoefficient P;
	if d <= 0 then return lift(a, ZZ);
	t := lift(d!*a, ZZ);
	return t + gotzmannBound(P - sum(toList(0..<t)/(i -> binomPoly(ring P, d - s - i, d))), t+s);
)

lexSegment = method()
lexSegment Ideal := Ideal => I -> (
	s := gotzmannBound(hilbertPolynomial(I, Projective => false), 0);
	R := ring I;
	L := ideal(0_R);
	hF := (flatten entries last coefficients hilbertSeries(I, Order => s+1))/(c -> lift(c, ZZ));
	for i from 1 to s do (
		A := basis(i, R/L);
		L = trim(L + ideal(submatrix'(lift(A, R), {(numcols A - hF#i)..numcols A})));
	);
	L
)

TEST ///
binomPoly(QQ[i], 3, 3)
///

TEST ///
R = QQ[x_0..x_3]
I = ideal(random(2, R), random(3, R)) -- canonical curve of genus 4
P = hilbertPolynomial(I, Projective => false)
gotzmannBound(P,0)
J = lexSegment I;
hilbertPolynomial(J, Projective => false) == hilbertPolynomial(I, Projective => false)
mingens J
///

TEST ///
R = QQ[x_0..x_4]
I = ideal(x_0^2, x_1^3)
P = hilbertPolynomial(I, Projective => false)
s = gotzmannBound(P,0)
time J = lexSegment I;
hilbertSeries(I, Order => s) == hilbertSeries(J, Order => s)
///