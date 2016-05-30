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
    lowerbound = 0;
    upperbound = oldflog; 
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

divideFraction = method();

-- This function takes in a fraction t and a prime p and spits out a list
-- {a,b,c}, where t = a/(p^b*(p^c-1))
-- if c = 0, then this means that t = a/p^b
divideFraction( ZZ, QQ ) := ( p, t ) -> 
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
     {a,b,c}
)

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
     flatten (apply(#prod, i -> toList(((prod#(i))#1):((prod#(i))#0)) ))
)

--===================================================================================

--Returns a list of all proper -- not one -- factors of number.
--Has funny order...
getFactorList = n ->
(
     if (n < 1) then error "getFactorList: expected an integer greater than 1.";
     powSet := nontrivialPowerSet(numberToPrimeFactorList(n)); 
     toList ( set apply(#powSet, i->product(powSet#i)) )
)


--===================================================================================

--*************************************************
--*************************************************
--Finding Numbers in Given Range
--*************************************************
--*************************************************

--===================================================================================

findNumberBetweenWithDenom = method(); 

--This function finds rational numbers in the range of the interval
--with the given denominator
findNumberBetweenWithDenom( ZZ, List ) := ( myDenom, myInterv ) ->
(
     upperBound := floor((myInterv#1)*myDenom)/myDenom; 
          --finds the number with denominator myDenom less than the upper 
	  --bound of myInterv
     lowerBound := ceiling((myInterv#0)*myDenom)/myDenom; 
          --finds the number with denominator myDenom greater than the lower
	  -- bound of myInterv
     if (upperBound >= lowerBound) then (
	  --first we check whether there is anything to search for
	  apply( 1+numerator((upperBound-lowerBound)*myDenom), i-> lowerBound+(i/myDenom) )
     )
     else {}
)

--===================================================================================

findNumberBetween = method(); 

--This function finds rational numbers in the range of 
--the interval; the max denominator allowed is listed. 
findNumberBetween( ZZ, List ) := ( maxDenom, myInterv )->
(
     divisionChecks :=  new MutableList from maxDenom:true; 
         -- creates a list with maxDenom elements all set to true.
     outList := {};
     i := maxDenom;
     while (i > 0) do (
	  if ((divisionChecks#(i-1)) == true) then --if we need to do a computation..
	      outList = join(outList,findNumberBetweenWithDenom(i,myInterv));
	  factorList := getFactorList(i);
     	  apply(factorList, j-> (divisionChecks#(j-1) = false) );
	  i = i - 1;
     );
     sort(toList set outList)
)

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
--Gives expansion in reverse...
basePExp( ZZ, ZZ ) := ( p, N ) ->
(
    if N < 0 then error "basePExp: Expected N to be positive";
    if N < p then {N} else prepend( N % p, basePExp(p, N // p))
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
    if x<0 then error "truncatedBasePExp: Expected x>0";
    ( ceiling( p^e*x ) - 1 )/p^e    	
)

--truncation threads over lists.
truncatedBasePExp ( ZZ, ZZ, List ) := ( p, e, u ) -> apply( u, x -> truncatedBasePExp( p, e, x ) )

--===================================================================================

--- write n=a*p^e+a_{e-1} p^{e-1} + \dots + a_0 where 0\leq e_j <p 
baseP1 = ( n, p, e ) ->
(
    a:=n//(p^e);
    answer:=1:a;
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


--Given a vector w={x,y}, x and y rational in [0,1], returns a number of digits 
--such that it suffices to check to see if x and y add without carrying in base p
carryTest = ( p, w ) ->
(
     if w#0 < 0 or w#0 > 1 or w#1 < 0 or w#1 > 1 then error "basePExp: Expected w in [0,1]^2";
     c := 0; for i from 0 to #w-1 do c = max(c, (divideFraction(p, w#i))#1);
     d := 1; for j from 0 to #w-1 do if ((divideFraction(p, w#j))#2)!=0 then d = lcm(d,(divideFraction(p,w#j))#2);
     c+d+1
)

--Given a vector w={x,y} of rational integers in [0,1], returns the first spot 
--e where the x and y carry in base p; i.e., 
--(the e-th digit of x)+(the e-th digit of y) >= p
-- counting is weird...
firstCarry = ( p, w ) ->
(   
    if w#0 < 0 or w#0 > 1 or w#1 < 0 or w#1 > 1 then error "firstCarry: Expected w in [0,1]^2";
    i:=0;
    d:=0;
    carry:=0;
    zeroTest := false;
    for j from 0 to #w-1 do if w#j == 0 then zeroTest=true;
    if zeroTest == true then carry = -1 else
    (
	i = 0; while d < p and i < carryTest(p,w) do 
	(
	    i = i + 1;
	    d = 0; for j from 0 to #w-1 do  d = d + digit(p, i,w#j);
	);
        if i == carryTest(p,w) then i = -1;
        carry = i;
     );
     carry
)

--===================================================================================

--Returns a vector of the reciprocals of the entries of a vector.
--Mutable list....
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
    apply( n, j-> if i==j then 1 else 0 )
)
 
--===================================================================================

getNumAndDenom = method()

-- Takes a rational vector u and returns a pair (a,q), where a
--is an integer vector and q an integer such that u=a/q.
getNumAndDenom ( List ) := u -> 
(
    den := lcm apply( u, denom );
    a := apply( u, n-> lift( n*den, ZZ ) );
    ( a, den )        
)

--===================================================================================

taxicabNorm = method()

--Computes the taxicab norm of a vector.
taxicabNorm ( List ) := u -> sum( u, abs )

--===================================================================================

--Finds the x-intercept of a line passing through two points
xInt = ( x1, y1, x2, y2 ) ->
(
    if x1 == x2 then error "xInt: x1==x2 no intersection";
    x1-(y1/((y1-y2)/(x1-x2)))
)
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

-- maxIdeal returns the ideal generated by the variables of a polynomial ring
maxIdeal = method()

maxIdeal ( PolynomialRing ) := R -> monomialIdeal R_*


