newPackage( "RationalMaps",
Version => "0.1", Date => "May 7th, 2016", Authors => {
     {Name => "Karl Schwede",
     Email=> "kschwede@gmail.com",
     HomePage=> "http://www.math.utah.edu/~schwede"
     },
     {Name => "Daniel Smolkin",
     Email=> "smolkin@math.utah.edu",
     HomePage=> "http://www.math.utah.edu/~smolkin"
     },
     {Name => "S. Hamid Hassanzadeh",
     Email => "hassanzadeh.ufrj@gmail.com"},
     {Name => "C.J. Bott",
     Email => "cjamesbott@gmail.com"}
}, --this file is in the public domain
Headline => "A package for working with Weil divisors.", DebuggingMode => true, Reload=>true)
export{
	"isBirationalMap",
	"imageOfMap",
	"baseLocusOfMap",
	"dimImage",
	"isRegular",
	"invertMap",
	"isEmbedding",
	"blowUpIdeals",
	"relationType",
	"dgi",
	"isSameDegree"
}

----------------------------------------------------------------
--************************************************************--
--Structure of our divisor objects and their display------------
--************************************************************--
----------------------------------------------------------------

imageOfMap = method();
imageOfMap(Matrix,Ideal,Ideal) := (f,a,b) -> (
	h := map((ring a)/a,(ring b)/b,f);
	-- the image of f is the same as the kernel of its pullback on the 
	-- coordinate rings. h is this pullback
	im := ker h;
	im
	);

dimImage = method();
dimImage(Matrix,Ideal,Ideal) := (f,a,b) ->(
	I := imageOfMap(f,a,b);
	dim I - 1
	-- substract 1 from the dimension of the image since in projective space
	);

baseLocusOfMap = method();

baseLocusOfMap(Matrix) := (L1) -> ( --L1 is a row matrix
    M:= gens ker transpose presentation image L1;
    -- this matrix gives all the "equivalent"
    -- ways to write the map in question (e.g. (xy : xz) is 
    -- equivalent to (y : z) ). So we do this to get the 
    -- representation of our map that's defined on the biggest
    -- set of points (e.g. (y : z) extends (xy : xz) to the locus where
    -- x is zero). C.f. proposition 1.1 of the paper
    -- "Cremona Transformations and some Related Algebras" by Aron Simis, 
    -- J. Algebra 280 (2004)
    
    
    L:= apply(entries M, ll->ideal(ll));
    saturate fold(L, plus)
    -- these two commands create an ideal for the base 
    -- locus from the information
    -- given in the matrix above. We take the saturation to get
    -- the biggest ideal that gives the same variety. 
    
);



 blowUpIdeals=method();
  
  --this is to compute the ideal definition of the blowup of a subvariety Z 
  --in the projective variety X
  --the defining ideal of the variety X=a
  --the subvariety Z = L list of element in the ring of X
  
  blowUpIdeals(Ideal, BasicList):=(a,L)->(
    r:=length  L;
    S:=ring L_0;
    n:=numgens ambient  S;
    K:=coefficientRing S;
    yyy:=local yyy;
    ttt:=local ttt;
    mymon:=monoid[({ttt}|gens ambient S|toList(yyy_0..yyy_(r-1))), MonomialOrder=>Eliminate 1];
    tR:=K(mymon);
   -- tR:=K[t,gens ambient S,vars(0..r-1),   MonomialOrder=>Eliminate 1];
    f:=map(tR,S,submatrix(vars tR,{1..n}));
    F:=f(matrix{L});
    myt:=(gens tR)#0;
    J:=ideal(f(gens a))+ideal apply(1..r,j->(gens tR)_(n+j)-myt*F_(0,(j-1)));
    L2:=ideal selectInSubring(1,gens gb J);
    W:=local W;
    nextmon:=monoid[(gens ambient  S|toList(W_0..W_(r-1))), Degrees=>{n:{1,0},r:{0,1}}];
    R:=K(nextmon);
    g:=map(R,tR,0|vars R);
    trim g(L2)); 

blowUpIdeals(Ideal, Ideal):=(a,b)->(
    blowUpIdeals(a, first entries gens b)
    );

--Matrix M consists of elements in the ring of a 
blowUpIdeals(Ideal, Matrix):=(a,M)->(
    blowUpIdeals(a, first entries gens M)
    );


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 relationType=method();
 --this function computes the "relation tye" of an ideal in a ring R.
 --Let R be the ring given bythe  ideal a and L be a list of elements in R.
 --the relation type is the biggest degree in terms of new variables in the
 --defining ideal of the rees algebra of I over R. 
 --  
 
 relationType(Ideal,BasicList):=(a,L)->(
     S:=ring L_0;
     J:=blowUpIdeals(a,L);
     n:=numgens J;
     L2:={};
     for i from 0 to n-1 do L2=append(L2,(degree J_i)_1);
     max L2);
 
 relationType(Ideal,Ideal):=(a,b)->(
     relationType(Ideal,first entries gens b)
     );
 
 relationType(Ring,Ideal):=(R1,b)->(
     relationType(Ideal R1,first entries gens b)
     );
     
     
     
 --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dgi=method();
 --this function compute the degeneration index of an ideal a which is the 
 --number of t linear generators among the generators of a.
 --dgi measures the number of hyperPlanes which cut the variety
 -- defined by a.
  

 dgi(Ideal):=(a)->(
     S := ring a; 
     n:=numgens a;
     d:=0;
     for i from 0 to n-1 do (
	 if (a_i != sub(0, S)) then (
	     if (degree a_i)=={1} then d=d+1
	 );
     );
     d);
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isSameDegree=method();
isSameDegree(BasicList):=(L)->(
   n:=#L;
   flag := true;
   if n!=0 then (
       d:=degree L#0;
       i := 0;
       while ((i < n) and (flag == true)) do(
	   if (isHomogeneous(L#i) == false) then flag = false;
	   if (degree(L#i) != d) then flag = false;
       	   i = i+1;
       );
    );  
    flag
);
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isBirationalMap = method();

--this checks whether a map X -> Y is birational.

--X = Proj R
--Y = Proj S
--This madfp is given by a list of elements in R, all homogeneous
--of the same degree.  
--Below we have defining ideal of X = di
--defining ideal of Y = im
--list of elements = bm
isBirationalMap(Ideal,Ideal,BasicList) :=(di,im,bm)->(
    if isSameDegree(bm)==false then error "Expected a list of homogenous elements of the same degree";
    R:=ring di;
    r:=numgens ambient R;
    K:=coefficientRing R;
    Jr:= blowUpIdeals(di,bm);
    n:=numgens Jr;
    L:={};
    for i from 0 to (n-1) do if  (degree Jr_i)_0==1 then
      L=append(L, Jr_i);
   JD:=diff(transpose ((vars ambient ring Jr)_{0..(r-1)}) ,gens ideal L);
   S:=ring im;
   vS:=gens ambient S;
   g:=map(S,ring Jr, toList(apply(0..r-1,z->0))|vS);
   barJD:=g(JD);
   jdd:=(numgens ambient R)-1+dgi(di);
   not(isSubset(minors(jdd,barJD),im))
   );    

isBirationalMap(Ring,Ring,BasicList) := (R1, S1, bm)->(
    isBirationalMap(ideal R1, ideal S1, bm)
    );

--isBirationalMap(Ideal,Ring,BasicList) := (di, S1, bm)->(
--    isBirationalMap(di, ideal S1, bm)
--    );

--isBirationalMap(Ring,Ideal,BasicList) := (R1, im, bm)->(
--    isBirationalMap(ideal R1, im, bm)
--    );

isBirationalMap(RingMap) :=(f)->(
    isBirationalMap(target f, source f, first entries matrix f)
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

doc ///
	Key 
		imageOfMap
	Headline
		Finds defining equations for the image of a rational map
	Usage
		image = imageOfMap(f,a,b)
	Inputs
		f:Matrix
			projective rational map given by polynomial represenative
		a:Ideal
			defining equations for X
		b:Ideal
			defining equations for Y
	Outputs
		im:Ideal
			defining equations for the image of f
	Description
		Text
			Defines the pullback map on the coordinate rings of X and Y. The kernel of this pullback map gives the image of the original map f
		Example
			S = QQ[x,y,z]
			a = ideal(x^2+y^2+z^2)
			T = QQ[u,v]
			b = ideal(u^2+v^2)
			f = matrix{{x*y,y*z}}
			imageOfMap(f,a,b)
///
			
doc ///
        Key
                dimImage
        Headline
                Computes dimension of image of rational map of projective varieties
        Usage
                dim = dimImage(f,a,b)
        Inputs
                f: Matrix
                        a rational map between 2 projective varieties, f: X -> Y
                        assumed that f is given by a polynomial representation
                a: Ideal
                        defining equations for X
                b: Ideal
                        defining equations for Y
        Outputs
                dim:ZZ
			dimension of image
///

doc ///
    Key
        baseLocusOfMap
    Headline
        Computes base locus of a map from a projective variety to projective space
    Usage
        I = baseLocusOfMap(L)
    Inputs
        L: Matrix
            Row matrix whose entries correspond to the coordinates of your map to projective space.
    Outputs
        I: Ideal
            The saturated defining ideal of the baselocus
///

TEST ///
	------------------------------------
	------- Tests for imageOfMap -------
	------------------------------------   

	S = QQ[x,y,z]
        a = ideal(x^2+y^2+z^2)
        T = QQ[u,v]
        b = ideal(u^2+v^2)
        f = matrix{{x*y,y*z}}
        image = imageOfMap(f,a,b)  
	assert(image == ideal(v^4,u*v^3))

	-------------------------------------
	-- Tests for baseLocusOfMap ---------
	-------------------------------------

    R = QQ[x,y,z]	
	M = matrix{{x^2*y, x^2*z, x*y*z}}
	I = ideal(x*y, y*z, x*z)
	assert(I == baseLocusOfMap(M))
	
/// 
----FUTURE PLANS------

