--*************************************************
--*************************************************
--*************************************************
--*************************************************
--This file is used for doing basic computations 
--i.e. things using only lists, numbers, etc.
-- that support other functions in the Fsing
--package.  
--*************************************************
--*************************************************
--*************************************************
--*************************************************

--*************************************************
--Basic Manipulations with Numbers 
--*************************************************
--===================================================================================

denom = method(); 

--Finds the denominator of a rational number.
--denom is always positive.
denom( QQ ) := x -> denominator(x); 

--Finds the denominator of an integer.
denom( ZZ ) := x -> 1; 

--===================================================================================

num = method(); 

--Finds the numerator of a rational number.
--Will be negative if x is negative.
num( QQ ) := x -> numerator(x); 

--Finds the numerator of an integer.
num( ZZ ) := x -> x; 

--===================================================================================

--Finds the fractional part of a number.
fracPart = x -> x - floor(x)

--===================================================================================

--Computes floor(log_b x), correcting problems due to rounding.
floorLog = ( b, x ) -> 
(
    if ( x < b ) then ( return 0 );
    flog := floor( log_b x );
    while b^flog <= x do flog = flog + 1;
    flog - 1       
)

-- the following function is faster than just floor(log_b x) if x is not too large compared
-- to b. 

-- for instance, the following function takes 1.2e-4 seconds to find log_2(13131231), 
-- whereas the old method takes 9.5e-5 seconds. 

fasterFloorLog = ( b, x ) -> (
    if ( x < b) then ( return 0 );
    flog := 1; 
    powerofb := b;
    oldpowerofb := 0;
    oldflog := 0;
    while powerofb <= x do (
        oldflog = flog;
        flog = flog * 2; -- just so we don't waste time dividing flog by 2 later
                         -- getting rid of this seems to make this function take
                         -- 7% more time. 

        oldpowerofb = powerofb; -- just so we don't waste time taking square roots
                               -- during the binary search part
                               -- taking the squareroot seems to slow down
                               -- this function by about 70%
        powerofb = powerofb^2;
    );

    -- binary search
    lowerbound := 0;
    upperbound := oldflog; 
    while (lowerbound + 1 < upperbound ) do ( --maybe the answer is between these two
        mid := ceiling ((lowerbound + upperbound)/2);
        if (oldpowerofb * b^mid > x)
        then ( upperbound = mid; )
        else ( lowerbound = mid; ); --if equality, now lowerbound is the answer
    );

    -- it's possible that x is == to b^upperbound, so we check. 
    if (oldpowerofb * b^upperbound == x) then (return oldflog + upperbound;);
    return oldflog + lowerbound;
)

--===================================================================================

multOrder = method()

--Finds the multiplicative order of a modulo b.
multOrder( ZZ, ZZ ) := ( a, b ) ->
(
    if gcd( a, b ) != 1 then error "multOrder: Expected numbers to be relatively prime.";
    n := 1;
    x := 1;
    while (x = (x*a) % b) != 1 do n = n+1;
    n	      
)     

--===================================================================================

divideFraction = method(Options => {NoZeroC=>false});

-- This function takes in a fraction t and a prime p and spits out a list
-- {a,b,c}, where t = a/(p^b*(p^c-1))
-- if c = 0, then this means that t = a/p^b
--alternately, if NoZeroC => true, then we will always write t = a/p^b(p^c - 1)
--even if it means increasing a. 
divideFraction( ZZ, QQ ) := o -> ( p, t ) -> 
(
    a := num t; -- finding a is easy, for now
    den := denom(t);
    b := 1;
    while den % p^b == 0 do b = b+1;
    b = b-1; 
    temp := denom( t*p^b );
    local c;
    if (temp == 1) then c = 0 else 
    (
        c = multOrder( p, temp );  
        a = lift( a*(p^c-1)/temp, ZZ ); -- fix a
    );
    if ((o.NoZeroC == true) and (c == 0)) then (
        a = a*(p-1);
        c = 1;
    );
    {a,b,c}
)
divideFraction( ZZ, ZZ ) := (p, t) -> divideFraction(p, t/1)

--===================================================================================
     
--Finds the a/p^e nearest t from above.
findNearPthPowerAbove = ( p, e, t ) -> ceiling( t*p^e )/p^e

--===================================================================================

--Finds the a/p^e nearest t from below.
findNearPthPowerBelow = ( p, e, t ) -> floor( t*p^e )/p^e

--===================================================================================

--*************************************************
--Information Regarding Factors and Factorization
--*************************************************

--===================================================================================

--Returns the power set of a  list removing the emptyset.  
nontrivialPowerSet = L -> delete( {}, subsets L )

--===================================================================================

--Returns a list of factors of a number with repeats.
numberToPrimeFactorList = n ->
(
     prod := factor n;
     flatten apply( toList prod, x -> toList( x#1:x#0 ) )
)

--===================================================================================

--Returns a list of all proper -- not one -- factors of number.
--Has funny order...
getFactorList = n ->
(
     if (n < 1) then error "getFactorList: expected an integer greater than 1.";
     powSet := nontrivialPowerSet( numberToPrimeFactorList( n ) ); 
     toList set apply( powSet, x -> product( x ) )
)

--===================================================================================

--*************************************************
--Finding Numbers in Given Range
--*************************************************

--===================================================================================

findNumberBetweenWithDenom = method(); 

--This function finds rational numbers in the range of the interval
--with the given denominator
findNumberBetweenWithDenom( ZZ, ZZ, ZZ) := ( myDenom, firstN, secondN) ->
(
     upperBound := floor((secondN)*myDenom)/myDenom; 
          --finds the number with denominator myDenom less than the second number
     lowerBound := ceiling((firstN)*myDenom)/myDenom; 
          --finds the number with denominator myDenom greater than the first number

     if (upperBound >= lowerBound) then (
	  --first we check whether there is anything to search for

	  apply( 1+numerator((upperBound-lowerBound)*myDenom), i-> lowerBound+(i/myDenom) )
     )
     else {}
)

--for backwards compatibility
--findNumberBetweenWithDenom( ZZ, List ) := (a, L) -> findNumberBetweenWithDenom(a, L#0, L#1);

--===================================================================================

findNumberBetween = method(); 

--This function finds rational numbers in the range of 
--the interval; the max denominator allowed is listed. 
findNumberBetween( ZZ, ZZ, ZZ) := ( maxDenom, firstN, secondN)->
(
     divisionChecks :=  new MutableList from maxDenom:true; 
         -- creates a list with maxDenom elements all set to true.
     outList := {};
     i := maxDenom;
     while (i > 0) do (
	  if ((divisionChecks#(i-1)) == true) then --if we need to do a computation..
	      outList = join(outList,findNumberBetweenWithDenom(i, firstN, secondN ));
	  factorList := getFactorList(i);
     	  apply(factorList, j-> (divisionChecks#(j-1) = false) );
	  i = i - 1;
     );
     sort(toList set outList)
)

--for backwards compatibility
--findNumberBetween( ZZ, List ) := ( maxDenom, myInterv )-> findNumberBetween( maxDenom, myInterv#0, myInterv#1);

--===================================================================================

--*************************************************
--Base p-expansions
--*************************************************

--===================================================================================

digit = method()

--Gives the e-th digit of the non-terminating base p expansion of x in [0,1].
digit ( ZZ, ZZ, QQ ) := ( p, e, x ) -> 
(
    if x < 0 or x > 1 then error "digit: Expected x in [0,1]";     	
    local y;
    if fracPart( p^e*x ) != 0 then y = floor( p^e*x ) - p*floor( p^(e-1)*x );
    if fracPart( p^e*x ) == 0 then y = floor( p^e*x ) - p*floor( p^(e-1)*x ) - 1;
    if fracPart( p^(e-1)*x ) == 0 then y = p-1;
    y	  
)

--Creates list containing e-th digits of non-terminating base p expansion of list of numbers.
digit ( ZZ, ZZ, List ) := ( p, e, u ) -> apply( u, x -> digit( p, e, x ) )

--===================================================================================

basePExp = method(); 

--Computes the terminating base p expansion of a positive integer.
--Gives expansion in reverse... so from left to right it gives
--the coefficient of 1, then of p, then of p^2, and so on
basePExp( ZZ, ZZ ) := ( p, N ) ->
(
    if N < 0 then error "basePExp: Expected N to be positive";
    if N < p then {N} else prepend( N % p, basePExp(p, N // p)) 
    -- would this be faster if it were tail-recursive? we could do this w/ a helper function.
)

--Creates a list of the first e digits of the non-terminating base p expansion of x in [0,1].
basePExp( ZZ, ZZ, QQ ) := ( p, e, x ) -> 
(
    if x < 0 or x > 1 then error "basePExp: Expected x in [0,1]";
    apply( e, i -> digit( p, i+1, x ) )
)

--===================================================================================

truncatedBasePExp = method()

--Gives the e-th truncation of the non-terminating base p expansion of a rational number.
truncatedBasePExp ( ZZ, ZZ, QQ ) := ( p, e, x ) -> 
(
    if x < 0 then error "truncatedBasePExp: Expected x>0";
    ( ceiling( p^e*x ) - 1 )/p^e    	
)

--truncation threads over lists.
truncatedBasePExp ( ZZ, ZZ, List ) := ( p, e, u ) -> apply( u, x -> truncatedBasePExp( p, e, x ) )

--===================================================================================

--- write n=a*p^e+a_{e-1} p^{e-1} + \dots + a_0 where 0\leq a_j <p 
--- DS: so it's just like doing basePExp but giving up after p^e and just returning whatever number's left
--- DS: this could be merged with basePExp. Should it be? 
--- note: I changed the calling order here should change to be consistent with basePExp
--- The change I made was switching the order of the first two arguments
baseP1 = ( p, n, e ) ->
(
    a:=n//(p^e);
    answer:=1:a; -- this generates the list (a)
    m:=n-a*(p^e);
    f:=e-1; 
    while (f>=0) do
    (
        d:=m//(p^f);
        answer=append(answer,d);
        m=m-d*(p^f);
        f=f-1;
    );
    answer
)	

--===================================================================================

--*************************************************
--Manipulations with Vectors   
--*************************************************

--===================================================================================

--Given a vector w of rational integers in [0,1], returns a number of digits such that
--it suffices to check to see if the components of w add without carrying in base p
carryTest = ( p, w ) ->
(
    if any( w, x -> x < 0 or x > 1 ) then 
        error "carryTest: Expected the second argument to be a list of rational numbers in [0,1]";
     div := apply( w, x -> divideFraction(p, x) );
     c := max (transpose div)#1; --max of second components of div
     v := selectNonzero (transpose div)#2; -- nonzero third components of div
     d := if v === {} then 1 else lcm v;
     c+d+1
)

--===================================================================================

--Given a vector w of rational integers in [0,1], returns the first spot 
--e where the the sum of the entries in w carry in base p
firstCarry = ( p, w ) ->
(   
    if any( w, x -> x < 0 or x > 1 ) then 
        error "firstCarry: Expected the second argument to be a list of rational numbers in [0,1]";
    if product( w ) == 0 then -1 else
    (
	i := 0;	
	d := 0;
	while d < p and i < carryTest(p,w) do 
	(
	    i = i + 1;
	    d = sum digit( p, i, w )
	);
        if i == carryTest(p,w) then -1 else i
     )
)

--===================================================================================

--Returns a vector of the reciprocals of the entries of a vector.
-- Probably not needed anymore.
reciprocal = w ->
(
    if product(w) == 0 then error "reciprocal: entries of vector must be non-zero.";
    apply(w, i -> 1/i)
)

--===================================================================================

getCanVector = method()

--canVector(i,n) returns the i-th canonical basis vector in dimension n
--Warning: for convenience, this uses Macaulay2's convention of indexing lists starting 
--with 0; so, for example, {1,0,0,0} is canVector(0,4), not canVector(1,4).
getCanVector ( ZZ, ZZ ) := ( i, n ) -> 
(
    if ( (i<0) or (i>=n) ) then error "canVector(i,n) expects integers i and n with 0<=i<n.";   
    apply( n, j -> if i==j then 1 else 0 )
)
 
--===================================================================================

getNumAndDenom = method()

-- Takes a rational vector u and returns a pair (a,q), where a
--is an integer vector and q an integer such that u=a/q.
getNumAndDenom ( List ) := u -> 
(
    den := lcm apply( u, denom );
    a := apply( u, n -> lift( n*den, ZZ ) );
    ( a, den )        
)

--===================================================================================

taxicabNorm = method()

--Computes the taxicab norm of a vector.
taxicabNorm ( List ) := u -> sum( u, abs )

--===================================================================================

--Selects or finds positions of nonzero, zero, positive entries in a list
selectNonzero = L -> select( L, x -> x != 0 )
selectPositive = L -> select( L, x -> x > 0 )
nonzeroPositions = L -> positions( L, x -> x != 0 )
zeroPositions = L -> positions( L, x -> x == 0 )

--===================================================================================

--*************************************************
--Tests for various types of polynomials   
--*************************************************

--===================================================================================

--isPolynomial(F) checks if F is a polynomial
isPolynomial = method( TypicalValue => Boolean )

isPolynomial (RingElement) := Boolean => F -> isPolynomialRing( ring F ) 

--===================================================================================

--isPolynomialOverPosCharField(F) checks if F is a polynomial over a field
--of positive characteristic
isPolynomialOverPosCharField = method( TypicalValue => Boolean )

isPolynomialOverPosCharField (RingElement) := Boolean => F ->
    isPolynomial F and isField( kk := coefficientRing ring F ) and ( char kk ) > 0

--===================================================================================

--isPolynomialOverFiniteField(F) checks if F is a polynomial over a finite field.
isPolynomialOverFiniteField = method( TypicalValue => Boolean )

isPolynomialOverFiniteField (RingElement) := Boolean => F ->
    isPolynomialOverPosCharField( F ) and isFinitePrimeField(coefficientRing ring F)

--===================================================================================

--Determines whether a polynomial f is a diagonal polynomial (i.e., of the form 
--x_1^(a_1)+...+x_n^(a_n)) over a field of positive characteristic 
isDiagonal = method( TypicalValue => Boolean )

isDiagonal (RingElement) := Boolean => f ->
    isPolynomialOverPosCharField( f ) and 
    ( product( exponents( f ), v -> #(positions( v, x -> x != 0 )) ) == 1 )

--===================================================================================

--Returns true if the polynomial is a monomial
isMonomial = method( TypicalValue => Boolean )

isMonomial (RingElement) := Boolean => f -> 
    isPolynomial f and #( terms f ) == 1

--===================================================================================

--Returns true if the polynomial is a binomial over a field of positive characteristic
isBinomial = method( TypicalValue => Boolean )

isBinomial (RingElement) := Boolean => f -> 
    isPolynomialOverPosCharField f and #( terms f ) == 2

--===================================================================================
  
--isBinaryForm(F) checks if F is a homogeneous polynomial in two variables.
--WARNING: what we are really testing is if the *ring* of F is a polynomial ring in two 
--variables, and not whether F explicitly involves two variables. (For example, if F=x+y 
--is an element of QQ[x,y,z], this test will return "false"; if G=x is an element of 
--QQ[x,y], this test will return "true".)
isBinaryForm = method( TypicalValue => Boolean )

isBinaryForm (RingElement) := Boolean => F ->
    isPolynomial F and numgens ring F == 2 and isHomogeneous F 

--===================================================================================

--isNonconstantBinaryForm(F) checks if F is a nonconstant homogeneous polynomial in two 
--variables. See warning under "isBinaryForm".
isNonConstantBinaryForm = method( TypicalValue => Boolean )

isNonConstantBinaryForm (RingElement) := Boolean => F -> 
    isBinaryForm F  and ( degree F )#0 > 0

--===================================================================================

--isLinearBinaryForm(F) checks if F is a linearform in two variables. See warning 
--under "isBinaryForm".
isLinearBinaryForm = method( TypicalValue => Boolean )

isLinearBinaryForm (RingElement) := Boolean => F -> 
    isBinaryForm F and ( degree F )#0 == 1

--===================================================================================

--*************************************************
--Partitions
--*************************************************

---------------------------------------------------------------------------------------
--- The following code was written in order to more quickly compute eth roots of (f^n*I)
--- It is used in fancyEthRoot
----------------------------------------------------------------------------------------
--- Find all ORDERED partitions of n with k parts
allPartitions = ( n, k )->
(
	PP0:=matrix{ toList(1..k) };
	PP:=mutableMatrix PP0;
	allPartitionsInnards (n,k,PP,{})
)

allPartitionsInnards = ( n, k, PP, answer)->
(
	local i;
	if (k==1) then 
	(
		PP_(0,k-1)=n;
		answer=append(answer,first entries (PP));
	)
	else
	(
		for i from 1 to n-(k-1) do
		(
			PP_(0,k-1)=i;
			answer=allPartitionsInnards (n-i,k-1,PP,answer)	;	
		);
	);
	answer
)

--===================================================================================

--*************************************************
--Miscelaneous
--*************************************************

--===================================================================================

-- maxIdeal returns the ideal generated by the variables of a polynomial ring
maxIdeal = method( TypicalValue => Ideal )

maxIdeal ( PolynomialRing ) := Ideal => R -> monomialIdeal R_*

maxIdeal ( RingElement ) := Ideal => f -> maxIdeal ring f

maxIdeal ( Ideal ) := Ideal => I -> maxIdeal ring I

--===================================================================================

--Finds the x-intercept of a line passing through two points
xInt = ( x1, y1, x2, y2 ) ->
(
    if x1 == x2 then error "xInt: x1==x2 no intersection";
    x1-(y1/((y1-y2)/(x1-x2)))
)

--===================================================================================

