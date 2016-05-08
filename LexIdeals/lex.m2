
binomPoly = method()
binomPoly (Ring, ZZ, ZZ) := RingElement => (R, a, b) -> (
	product(toList(1..b)/(i -> (first gens R + a + 1 - i)/i))
)

gotzmannBound = method()
gotzmannBound (RingElement, ZZ) := List => (P, s) -> (
	d := (degree P)#0;
	a := leadCoefficient P;
	if d <= 0 then return lift(a, ZZ);
	if a < 0 then error "Polynomial has no Gotzmann representation.";
	t := lift(d!*a, ZZ);
	return t + gotzmannBound(P - sum(toList(0..<t)/(i -> binomPoly(ring P, d - s - i, d))), t+s);
)
gotzmannBound RingElement := List => P -> gotzmannBound(P, 0)
gotzmannBound Ideal := List => I -> gotzmannBound hilbertPolynomial(I, Projective => false)

lexSegment = method()
lexSegment Ideal := Ideal => I -> (
	s := max(append(flatten degrees I, gotzmannBound I)); -- needed if I is not saturated
	print("Gotzmann bound: " | s);
	R := ring I;
	L := ideal(0_R);
	hF := (flatten entries last coefficients hilbertSeries(I, Order => s+1))/(c -> lift(c, ZZ));
	print "Finished finding Hilbert function values up to Gotzmann bound. Computing lex ideal ...";
	for i from 1 to s do (
		A := basis(i, R/L);
		L = trim(L + ideal(submatrix'(lift(A, R), {(numcols A - hF#i)..numcols A})));
	);
	L
)

TEST ///
R = QQ[i]
assert(binomPoly(R, 3, 3) == 1/6*i^3 + i^2 + 11/6*i + 1)
///

TEST ///
R = QQ[x_0..x_2]
I = intersect(ideal(x_2), (ideal gens R)^5) -- non-saturated ideal, coordinate hyperplane in P^2
time L = lexSegment I;
assert(hilbertSeries I === hilbertSeries L)
///

TEST ///
R = QQ[x_0..x_3]
I = ideal(random(2, R), random(3, R)) -- canonical curve of genus 4
time L = lexSegment I;
assert(hilbertSeries I === hilbertSeries L)
///