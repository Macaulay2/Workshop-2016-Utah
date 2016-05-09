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
denom( QQ ) := x -> denominator x; 

--Finds the denominator of an integer.
denom( ZZ ) := x -> 1; 

--===================================================================================

num = method(); 

--Finds the numerator of a rational number.
num( QQ ) := x -> numerator x; 

--Finds the numerator of an integer.
num( ZZ ) := x -> x; 

--===================================================================================

--Finds the fractional part of a number.
fracPart = ( x ) -> (x - floor(x))

--===================================================================================

--Computes floor(log_b x), correcting problems due to rounding.
floorLog = ( b, x ) -> 
(
    flog := floor( log_b x );
    while b^flog <= x do flog = flog + 1;
    flog - 1       
)

--===================================================================================

multOrder = method()

--Finds the multiplicative order of a modulo b.
multOrder( ZZ, ZZ ) := ( a, b ) ->
(
    if gcd( a, b ) != 1 then error "Expected numbers to be relatively prime.";
    n := 1;
    x := 1;
    while  (x = (x*a) % b) != 1  do n = n+1;
    n	      
)     

--===================================================================================

divideFraction = method();

-- This function takes in a fraction t and a prime p and spits out a list
-- {a,b,c}, where t = (a/p^b)(1/(p^c-1))
-- if c = 0, then this means that t = (a/p^b)
divideFraction( ZZ, QQ ) := ( p, t ) -> 
(
     a := num t; -- finding a is easy, for now
     den:=denom(t);
     b := 1;
     while den%p^b==0 do b=b+1;
     b = b-1; 
     temp := denom(t*p^b); --find the p^c-1 part of the denominator
     pow := 0; --we will look around looking for the power of pp that is 1 mod temp. 
     done := false; --when we found the power, this is set to true.
     if (temp == 1) then done = true; --if there is nothing to do, do nothing.
     while (done==false)  do (
          pow = pow + 1;
	  if (p^pow % temp == 1) then done = true
     );
     c := pow; --we found c, now we return the list
     if (c > 0) then a = lift(a*(p^c-1)/temp, ZZ); --after we fix a
     {a,b,c}
)

--===================================================================================
     
--Finds the a/pp^e1 nearest t1 from above.
findNearPthPowerAbove = ( p, e, t ) -> 
(
     ceiling(t*p^e)/p^e
)

--===================================================================================

--Finds the a/pp^e1 nearest t1 from below.
findNearPthPowerBelow = ( p, e, t ) -> 
(
     floor(t*p^e)/p^e
)

--===================================================================================

--*************************************************
--Information Regarding Factors and Factorization
--*************************************************

--===================================================================================

--Returns the digits in nn which are nonzero in binary 
--for example, 5 in binary is 101, so this would return {0,2}
--the second term tells me where to start the count, so passing
--5,0 gives {0,2} but 5,1 is sent to {1,3}.  i should be
--used only for recursive purposes
getNonzeroBinaryDigits = ( i, n ) -> 
(
    halfsies := n//2;
    val1 := n%2;
    val2 := false; 
    if (halfsies > 0) then val2 = (getNonzeroBinaryDigits(i+1,n//2));
    if ( (val1 != 0) and (not (val2 === false))) then (
	 flatten {i, val2}
    )
    else if (val1 != 0) then (
	 {i}
    )
    else if ( not (val2 === false)) then (
	 flatten {val2}
    )
    else(
	 false
    )
)

--===================================================================================


getSublistOfList = method();

--Returns the entries of myList specified by their position.
getSublistOfList( List, List ) := ( entryList, myList ) -> 
(
     if entryList == {} then(
	  error "getSublistOfList expected non-empty list" 
	  )
     else(
	 apply(#entryList, i->myList#(entryList#i))
    )
)

--===================================================================================

--Returns the power set of a  list removing the emptyset.  
nontrivialPowerSet = ( myList ) ->
(
     apply(2^(#myList)-1, i-> getSublistOfList(getNonzeroBinaryDigits(0,i+1),myList))
)

--===================================================================================

--Returns a list of factors of a number with repeats.
numberToPrimeFactorList = ( n )->
(
     prod := factor n;
     flatten (apply(#prod, i -> toList(((prod#(i))#1):((prod#(i))#0)) ))
)

--===================================================================================

--Returns a list of all proper factors of number.
getFactorList = ( n ) ->
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
findNumberBetweenWithDenom( ZZ, List ) := ( myDenom, myInterv )->
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
     else(
	  {}
     )
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
     	  apply(#factorList, j-> (divisionChecks#( (factorList#j)-1) = false) );
	  i = i - 1;
     );
     sort(toList set outList)
)

--===================================================================================

--*************************************************
--Base p-expansions
--*************************************************

--===================================================================================

basePExp = method(); 

--Computes the terminating base p expansion of an integer.
basePExp( ZZ, ZZ ) := ( p, N ) ->
(
    if N < p then return {N};
    prepend( N % p, basePExp(p, N // p))
)

--Computes terminating base p expansion of an integer 
--from digits zero to e-1 (little-endian first).
basePExp( ZZ, ZZ, ZZ ) := ( p, f, N ) ->
(
    e:=f-1;
    E:=new MutableList;
    scan(0..e,i-> 
    	(
     	    a := N//p^(e-i);
     	    E#(e-i) = a;
     	    N = N - (a)*p^(e-i);
    	)
    );
    new List from E
)

--Creates a list of the first e digits of the non-terminating base p expansion of x in [0,1].
basePExp( ZZ, ZZ, QQ ) := ( p, e, x ) -> 
(
    if x < 0 or x > 1 then(
	 error "basePExp: Expected x in [0,1]"
     )
     else(
     	 L := new MutableList;
     	 for i from 0 to e-1 do L#i = digit(p,i+1,x);
     	 L
     )
)

--===================================================================================

digit = method()

--Gives the e-th digit of the non-terminating base p expansion of x in [0,1].
digit ( ZZ, ZZ, QQ ) := ( p, e, x) -> 
(
    if x < 0 or x > 1 then(
	 error "digit: Expected x in [0,1]"
	 )
    else(
     	y := 0;
     	if fracPart(p^e*x) != 0 then y = floor(p^e*x) - p*floor(p^(e-1)*x);
     	if fracPart(p^e*x) == 0 then y = floor(p^e*x) - p*floor(p^(e-1)*x) - 1;
     	if fracPart(p^(e-1)*x) == 0 then y = p-1;
     	y	  
     )
)

--Creates list containing e-th digits of non-terminating base p expansion of list of numbers.
digit ( ZZ, ZZ, List ) := ( p, e, u ) -> 
(
    apply(u,x->digit(p,e,x))
)

--===================================================================================

truncatedBasePExp = method()

--Gives the e-th truncation of the non-terminating base p expansion of a rational number.
truncatedBasePExp ( ZZ, ZZ, QQ ) := ( p, e, x ) -> 
(
    if x<0 then (error "truncation: Expected x>0")
    else (ceiling(p^e*x)-1)/p^e    	
)

--truncation threads over lists.
truncatedBasePExp ( ZZ, ZZ, List) := (p,e,u) -> 
(
    apply(u,x->truncation(p,e,x))
)

--===================================================================================

--- write n=a*p^e+a_{e-1} p^{e-1} + \dots + a_0 where 0\leq e_j <p 
baseP1 = ( n, p, e )->
(
	a:=n//(p^e);
	answer:=1:a;
	m:=n-a*(p^e);
	f:=e-1; 
	while (f>=0) do
	{
		d:=m//(p^f);
		answer=append(answer,d);
		m=m-d*(p^f);
		f=f-1;
	};
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
     c := 0; for i from 0 to #w-1 do c = max(c, divideFraction(p,w#i)#1);
     d := 1; for j from 0 to #w-1 do if (divideFraction(p, w#j)#2)!=0 then d = lcm(d,dividFraction(p,w#j)#2);
     c+d+1
)

--Given a vector w={x,y} of rational integers in [0,1], returns the first spot 
--e where the x and y carry in base p; i.e., 
--(the e-th digit of x)+(the e-th digit of y) >= p
firstCarry = ( p, w ) ->
(     
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
reciprocal = ( w ) ->
(
     v := new MutableList from w;
     for c from 0 to #w-1 do v#c = 1/w#c;
     v
)

--===================================================================================

getCanVector = method()

--canVector(i,n) returns the i-th canonical basis vector in dimension n
--Warning: for convenience, this uses Macaulay2's convention of indexing lists starting 
--with 0; so, for example, {1,0,0,0} is canVector(0,4), not canVector(1,4).
getCanVector ( ZZ, ZZ ) := ( i, n ) -> 
(
    if ((i<0) or (i>=n)) 
        then (error "canVector(i,n) expects integers i and n with 0<=i<n.");   
    apply(n,j->if i==j then 1 else 0)
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
taxicabNorm ( List ) := u -> 
(
    sum( u, abs )
)

--===================================================================================

--Finds the x-intercept of a line passing through two points
xInt = ( x1, y1, x2, y2 ) ->  x1-(y1/((y1-y2)/(x1-x2)))
 
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
	{
		PP_(0,k-1)=n;
		answer=append(answer,first entries (PP));
	}
	else
	{
		for i from 1 to n-(k-1) do
		{
			PP_(0,k-1)=i;
			answer=allPartitionsInnards (n-i,k-1,PP,answer)	;	
		};
	};
	answer
)

--===================================================================================
