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
	s := max(append(flatten degrees I, gotzmannBound I)); -- necessary when I is not saturated
	print("Gotzmann bound: " | s);
	R := ring I;
	L := ideal(0_R);
	hF := (flatten entries last coefficients hilbertSeries(I, Order => s+1))/(c -> lift(c, ZZ));
	print "Computing lex ideal ...";
	for i from 1 to s do (
		A := basis(i, R/L);
		if not hF#?i then ( L = trim(L + (ideal gens R)^i); break; );
		n := numcols A;
		A = lift(A, R);
		K := ideal 0_R;
		for j from 0 to floor(n/80000) do ( -- to prevent long sequence error
			K = K + ideal(A_{j*80000..min(n - 1 - hF#i, (j+1)*80000 - 1)});
		);
		L = L + K;
		--L = trim(L + ideal(submatrix'(lift(A, R), {(numcols A - hF#i)..numcols A})));
	);
	L
)

-------------------------------------------------------
--DOCUMENTATION
-------------------------------------------------------

doc ///
     Key
     	  gotzmannBound
	  (gotzmannBound, RingElement)
	  (gotzmannBound, RingElement, ZZ)
	  (gotzmannBound, Ideal)
     Headline
     	  the Gotzmann bound for a saturated ideal
     Usage
     	  gotzmannBound P
	  gotzmannBound(P, d)
	  gotzmannBound I
     Inputs
     	  P:RingElement
		a Hilbert polynomial
	  d:ZZ
		an integer
	  I:Ideal
		an Ideal
     Outputs
     	  s:ZZ
		the Gotzmann bound of P (or of hilbertPolynomial(I))
     Description
     	  Text
		Returns the Gotzmann bound for regularity of a saturated ideal with Hilbert polynomial P,
		which is the number of binomial coefficients in the Gotzmann representation of P. The 
		Gotzmann bound of hilbertPolynomial(I) is computed when an ideal I is given as input.
		The optional integer input (default: 0) is for internal use of the function in recursive calls.
	  Example
	       R = QQ[x,y,z]
	       I = ideal(random(2, R), random(3, R))
	       s = gotzmannBound I
	       P = hilbertPolynomial(I, Projective => false)
	       s == gotzmannBound(P, 0) 
	       s == gotzmannBound P
     SeeAlso
	  lexSegment
     Caveat
	  The leading coefficient of P must be positive.
	  
///

doc ///
     Key
     	  binomPoly
	  (binomPoly, Ring, ZZ, ZZ)
     Headline
     	  binomial coefficient as a polynomial
     Usage
     	  binomPoly(R, a, b)
     Inputs
	  R:Ring
		The univariate ring for numerical polynomials (usually QQ[i])
	  a:ZZ
		an integer
	  b:ZZ
		an integer
     Outputs
     	  P:RingElement
		the binomial coefficient {\tt (i + a) choose b}, as a polynomial in i
     Description
     	  Text
		Returns the binomial coefficient {\tt (i + a) choose b}, as an element in QQ[i].
	  Example
	       P = binomPoly(QQ[i], 2, 3)
     SeeAlso
     	  gotzmannBound
	  lexSegment
///

doc ///
     Key
     	  lexSegment
	  (lexSegment, Ideal)
     Headline
     	  the lex segment ideal with the same Hilbert function
     Usage
     	  lexSegment I
     Inputs
	  I:Ideal
		a homogeneous ideal
     Outputs
     	  L:Ideal
		the lex segment ideal with the same Hilbert function as I
     Description
     	  Text
		Returns the lex segment ideal with the same Hilbert function as I, for any homogeneous
		ideal I (Artinian or not). Does not rely on lexIdeal() method. Uses the larger of the max 
		generating degree of I and the Gotzmann bound of hilbertPolynomial(I) as a stopping 
		criterion for the largest degree in which to compute the lex segment ideal.
	  Example
	       R = QQ[x_0..x_4]
	       I = ideal(random(2, R), random(3, R))
	       time L = lexSegment I
     SeeAlso
     	  lexIdeal
///


TEST ///
R = QQ[i]
assert(binomPoly(R, 3, 3) == 1/6*i^3 + i^2 + 11/6*i + 1)
///

TEST ///
R = QQ[x_0..x_2]
I = intersect(ideal(x_2), (ideal gens R)^5) -- non-saturated ideal, coordinate hyperplane in P^2
assert(hilbertSeries I === hilbertSeries lexSegment I)
///

TEST ///
R = QQ[x_0..x_3]
I = ideal(random(2, R), random(3, R)) -- canonical curve of genus 4
assert(hilbertSeries I === hilbertSeries lexSegment I)
///