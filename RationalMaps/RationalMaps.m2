newPackage( "RationalMaps",
Version => "0.1", Date => "May 7th, 2016", Authors => {
     {Name => "Karl Schwede",
     Email=> "kschwede@gmail.com",
     HomePage=> "http://www.math.utah.edu/~schwede"
     },
}, --this file is in the public domain
Headline => "A package for working with Weil divisors.", DebuggingMode => true, Reload=>true)
export{
	"isBirationalMap",
	"imageOfMap",
	"baseLocusOfMap",
	"dimImage",
	"isRegular",
	"invertMap",
	"isEmbedding"
}

----------------------------------------------------------------
--************************************************************--
--Structure of our divisor objects and their display------------
--************************************************************--
----------------------------------------------------------------

dimImage = method();

dimImage(Ideal) := (I1) -> (
);


imageOfMap = method();
imageOfMap(Matrix,Ideal,Ideal) := (f,a,b) -> (
	h = map((ring a)/a,(ring b)/b,f);
	ker h
	);

--****************************************************--
--*****************Documentation**********************--
--****************************************************--

beginDocumentation();

doc /// 
	 Key
		RationalMaps
     Headline
     	A package for computations with rational maps.
     Description
    	Text   
    	 A package for computations with rational maps.
///

TEST ///
assert true
///



----FUTURE PLANS------

