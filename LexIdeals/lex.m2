restart
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

TEST ///
binomPoly(QQ[i], 3, 3)
///

TEST ///
R = QQ[x_0..x_3]
I = ideal(random(2, R), random(3, R)) -- canonical curve of genus 4
P = hilbertPolynomial(I, Projective => false)
gotzmannBound(P,0)
///

TEST ///
R = QQ[x_0..x_4]
I = ideal(x_0^2, x_1^3)
P = hilbertPolynomial(I, Projective => false)
gotzmannBound(P,0)
///