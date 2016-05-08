--*************************************************
--*************************************************
--This file is used for doing computations with 
--integers that support other functions in the Fsing
--package.  
--*************************************************
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
divideFraction( ZZ, QQ ) := ( pp, t1 ) -> 
(
     a := num t1; -- finding a is easy, for now
     den:=denom(t1);
     b := 1;
     while den%p^b==0 do b=b+1;
     b = b-1; 
     temp := denom(t1*pp^b); --find the p^c-1 part of the denominator
     pow := 0; --we will look around looking for the power of pp that is 1 mod temp. 
     done := false; --when we found the power, this is set to true.
     if (temp == 1) then done = true; --if there is nothing to do, do nothing.
     while (done==false)  do (
          pow = pow + 1;
	  if (pp^pow % temp == 1) then done = true
     );
     c := pow; --we found c, now we return the list
     if (c > 0) then a = lift(a*(pp^c-1)/temp, ZZ); --after we fix a
     {a,b,c}
)

--===================================================================================
     
--Finds the a/pp^e1 nearest t1 from above.
findNearPthPowerAbove = ( t1, pp, e1 ) -> 
(
     ceiling(t1*pp^e1)/pp^e1
)

--===================================================================================

--Finds the a/pp^e1 nearest t1 from below.
findNearPthPowerBelow = ( t1, pp, e1 ) -> 
(
     floor(t1*pp^e1)/pp^e1
)

--===================================================================================

--Returns the digits in nn which are nonzero in binary 
--for example, 5 in binary is 101, so this would return {0,2}
--the second term tells me where to start the count, so passing
--5,0 gives {0,2} but 5,1 is sent to {1,3}.  i should be
--used only for recursive purposes
getNonzeroBinaryDigits = ( nn, i ) -> 
(
    halfsies := nn//2;
    val1 := nn%2;
    val2 := false; 
    if (halfsies > 0) then val2 = (getNonzeroBinaryDigits(nn//2,i+1));
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
getSublistOfList( List, List ) = ( entryList, myList ) -> 
(
     apply( #entryList, i->myList#(entryList#i) )
)

--===================================================================================

--Returns the power set of a  list removing the emptyset.  
nontrivialPowerSet = ( myList ) ->
(
     apply(2^(#myList)-1, i-> getSublistOfList(myList, getNonzeroBinaryDigits(i+1,0) ) )
)

--===================================================================================

--Returns a list of factors of a number with repeats.
numberToPrimeFactorList = ( nn )->
(
     prod := factor nn;
     flatten (apply(#prod, i -> toList(((prod#(i))#1):((prod#(i))#0)) ))
)

--===================================================================================

--Returns a list of all proper factors of number.
getFactorList = ( nn ) ->
(
     if (nn < 1) then error "getFactorList: expected an integer greater than 1.";
     powSet := nontrivialPowerSet(numberToPrimeFactorList(nn)); 
     toList ( set apply(#powSet, i->product(powSet#i)) )
)


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
	      outList = join(outList,findNumberBetweenWithDenom(myInterv, i));
	  factorList := getFactorList(i);
     	  apply(#factorList, j-> (divisionChecks#( (factorList#j)-1) = false) );
	  i = i - 1;
     );
     sort(toList set outList)
)

--===================================================================================

basePExp = method(); 

--Computes the terminating base p expansion of an integer.
basePExp( ZZ, ZZ ) := ( p, N ) ->
(
    if N < p then return {N};
    prepend( N % p, basePExp( N // p, p ))
)

--Computes terminating base p expansion of an integer 
--from digits zero to e-1 (little-endian first).
basePExp( ZZ, ZZ, ZZ ) := ( p, e1, N ) ->
(
    e:=e1-1;
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
