
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

-------------------------------------------------------
--DOCUMENTATION gotzmannBound
-------------------------------------------------------

doc///
     Key
     	  gotzmannBound
	  (gotzmannBound,RingElement)
	  (gotzmannBound,RingElement,ZZ)
	  (gotzmannBound,Ideal)
     Headline
     	  the Gotzmann bound of a Hilbert polynomial
     Usage
     	  a=gotzmannBound(P,d) or a=gotzmannBound(P) or a=gotzmannBound(I)
     Inputs
     	  P:RingElement
		a Hilbert polynomial
	  d:ZZ
		an integer
	  I:Ideal
		an Ideal
     Outputs
     	  a:ZZ
		the Gotzmann bound of P (or of hilbertPolynomial(I))
     Description
     	  Text
		Returns the Gotzmann bound for regularity of a saturated ideal with Hilbert polynomial P,
		which is the number of binomial coefficients in the Gotzmann representation of P. The 
		Gotzmann bound of hilbertPolynomial(I) is computed when an ideal I is input. The optional
		integer input is meant for internal use of the function in recursive calls.
	  Example
	       a=gotzmannBound(P)
	       a=gotzmannBound(P,0)
	       a=gotzmannBound(I)
     SeeAlso
	  lexSegment     	  
///

doc///
     Key
     	  binomPoly
	  (binomPoly,Ring,ZZ,ZZ)
     Headline
     	  the Hilbert polynomial of k[x_1, ..., x_b](a)
     Usage
     	  p=binomPoly(R,a,b)
     Inputs
	  R:Ring
		The univariate ring for numerical polynomials (usually QQ[i])
	  a:ZZ
		an integer
	  b:ZZ
		an integer
     Outputs
     	  p:RingElement
		the Hilbert polynomial of k[x_1, ..., x_b](a)
     Description
     	  Text
		Returns the Hilbert polynomial of projective space P^(b-1), twisted by a, 
		as an element of QQ[i].
	  Example
	       p=binomPoly(QQ[i],3,3)
     SeeAlso
     	  gotzmannBound
	  lexSegment
///

doc///
     Key
     	  lexSegment
	  (lexSegment,Ideal)
     Headline
     	  the lex segment ideal with the same Hilbert function
     Usage
     	  L=lexSegment(I)
     Inputs
	  I:Ideal
		an Ideal
     Outputs
     	  L:Ideal
		the lex segment ideal with the same Hilbert function as I
     Description
     	  Text
		Returns the lex segment ideal with the same Hilbert function as I.  Does not rely
		on lexIdeal() method in the case of a non-artinian ideal. Uses the larger of the
		max generating degree of I and the Gotzmann bound of hilbertPolynomial(I) to 
		determine the lergest degree in which to compute the lex segment ideal.
	  Example
	       L=lexSegment(I)
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
time L = lexSegment I;
assert(hilbertSeries I === hilbertSeries L)
///

TEST ///
R = QQ[x_0..x_3]
I = ideal(random(2, R), random(3, R)) -- canonical curve of genus 4
time L = lexSegment I;
assert(hilbertSeries I === hilbertSeries L)
///