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
Headline => "A package for working with rational maps.", DebuggingMode => true, Reload=>true)
export{
	"isBirationalMap",
	"idealOfImageOfMap",
	"baseLocusOfMap",
	"dimImage",
	"isRegularMap",
	"isEmbedding",
	"blowUpIdeals",
	"relationType",
	"dgi",
	"isSameDegree",
	"isBirationalOntoImage",
	"nonZeroMinor",
	"inverseOfMap",
	"mapOntoImage"
}

----------------------------------------------------------------
--************************************************************--
-------------------- Function Defitions ------------------------
--************************************************************--
----------------------------------------------------------------


--***Karl:  I dislike how this returns stuff.   The image of a map is not an ideal.  Hence I renamed it
idealOfImageOfMap = method();

idealOfImageOfMap(Ideal,Ideal,Matrix) := (a,b,f) -> (
	h := map((ring a)/a, ring b ,f);
	-- the image of f is the same as the kernel of its pullback on the 
	-- coordinate rings. h is this pullback
        idealOfImageOfMap(h)
	);

idealOfImageOfMap(Ideal,Ideal,BasicList) := (a,b,g) -> (
	h:=map((ring a)/a, ring b ,g);
	idealOfImageOfMap(h)
	);

idealOfImageOfMap(Ring,Ring,Matrix) := (R,S,f) -> (
	h := map(R,S,f);
	idealOfImageOfMap(h)	
	);

idealOfImageOfMap(Ring,Ring,BasicList) := (R,S,g) -> (
        h := map(R,S,g);
        idealOfImageOfMap(h)
	);

idealOfImageOfMap(RingMap) := (p) -> (
        --idealOfImageOfMap(target p, source p, first entries matrix p)
        h := map(target p, ambient source p,p);
        im := ker h;
	im
	);

dimImage = method();

dimImage(Ideal,Ideal,Matrix) := (a,b,f) -> (
	I := idealOfImageOfMap(a,b,f);
	dim I - 1
	-- substract 1 from the dimension of the image since in projective space
	);

dimImage(Ideal,Ideal,BasicList) := (a,b,g) -> (
	dimImage(a,b,g)
	);

dimImage(Ring,Ring,Matrix) := (R,S,f) -> (
	dimImage(ideal R, ideal S, f)
	);

dimImage(Ring,Ring,BasicList) := (R,S,g) -> (
	dimImage(ideal R, ideal S, g)
	);

dimImage(RingMap) := (p) -> (
	dimImage(target p, source p, first entries matrix p)
	);

baseLocusOfMap = method();

baseLocusOfMap(Matrix) := (L1) -> ( --L1 is a row matrix
    --maybe check all the maps in L1 are of the same degree?

    --just need to convert L1 to a basic list, I guess
    --if isSameDegree(L1)==false then error "Expected a matrix of homogenous elements of the same degree";
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

baseLocusOfMap(List) := (L) ->(
    baseLocusOfMap(matrix{L})
);

baseLocusOfMap(RingMap) := (ff) ->(
    mm := sub(matrix ff, target ff);  
    baseLocusOfMap(mm)
);

isRegularMap = method();

isRegularMap(Matrix) := (L1) -> ( --L1 is a row matrix
    I:= baseLocusOfMap(L1);
    I == ideal 1_(ring I)
);

isRegularMap(RingMap) := (ff) ->(
        I:=baseLocusOfMap(ff);
        I == ideal 1_(ring I)
);

 blowUpIdeals=method();
  
  --this is to compute the ideal definition of the blowup of a subvariety Z 
  --in the projective variety X
  --the defining ideal of the variety X=a
  --the subvariety Z = L list of element in the ring of X
  
  blowUpIdeals(Ideal, BasicList):=(a,L)->(
    r:=length  L;
    SS:=ring a;
    LL:=apply(L,uu->sub(uu, SS));
    n:=numgens ambient  SS;
    K:=coefficientRing SS;
    yyy:=local yyy;
    ttt:=local ttt;
    mymon:=monoid[({ttt}|gens ambient SS|toList(yyy_0..yyy_(r-1))), MonomialOrder=>Eliminate 1];
    tR:=K(mymon);
   -- tR:=K[t,gens ambient SS,vars(0..r-1),   MonomialOrder=>Eliminate 1];
    f:=map(tR,SS,submatrix(vars tR,{1..n}));
    F:=f(matrix{LL});
    myt:=(gens tR)#0;
    J:=sub(a,tR)+ideal apply(1..r,j->(gens tR)_(n+j)-myt*F_(0,(j-1)));
    L2:=ideal selectInSubring(1,gens gb J);
    W:=local W;
    nextmon:=monoid[(gens ambient  SS|toList(W_0..W_(r-1))), Degrees=>{n:{1,0},r:{0,1}}];
    RR:=K(nextmon);
    g:=map(RR,tR,0|vars RR);
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
 
  --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isBirationalOntoImage = method();

isBirationalOntoImage(Ideal,Ideal, BasicList) :=(di,im,bm)->(
      tar:=idealOfImageOfMap(di,im, bm);
      isBirationalMap(di,tar,bm)
      );
  
 isBirationalOntoImage(Ring,Ring,BasicList) := (R1, S1, bm)->(
    isBirationalMap(ideal R1, ideal S1, bm)
    ); 

isBirationalOntoImage(RingMap) :=(f)->(
    isBirationalOntoImage(target f, source f, first entries matrix f)
    );
    
    
    --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nonZeroMinor=method();
nonZeroMinor(Matrix,ZZ):=(M,ra)->(
    cc:=numColumns(M);
    ro:=numRows(M);
    col:=apply(0..cc-1,i->i);
    row:=apply(0..ro-1,i->i);
    Collist:=subsets(col,ra);
    Rowlist:=subsets(row,ra);
    nzlist:={};
    for i from 0 to (#Collist)-1 do  (
       if nzlist!={} then break; 
       for j from 0 to (#Rowlist)-1 do (
	   if det(submatrix(M,Rowlist#j,Collist#i))!=0 then (
	       nzlist=append(nzlist,{Collist#i,Rowlist#j});
	       break;
	   );
       );
   );
flatten nzlist);  
   
 --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
inverseOfMap = method();

--this checks whether a map X -> Y is birational.

--X = Proj R
--Y = Proj S
--This madfp is given by a list of elements in R, all homogeneous
--of the same degree.  
--Below we have defining ideal of X = di
--defining ideal of Y = im
--list of elements = bm
inverseOfMap(Ideal,Ideal,BasicList) :=(di,im,bm)->(
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
   if not(isSubset(minors(jdd,barJD),im))==false then error "The map is not Birational";
   Col:=(nonZeroMinor(barJD,jdd))#0;
   SbarJD:=submatrix(barJD,,Col);
   Inv:={};
   for i from 0 to jdd do Inv=append(Inv,(-1)^i*det(submatrix'(SbarJD,{i},)));
   map(S/im, R/di, Inv)
);    

inverseOfMap(Ring,Ring,BasicList) := (R1, S1, bm)->(
    inverseOfMap(ideal R1, ideal S1, bm)
    );

--isBirationalMap(Ideal,Ring,BasicList) := (di, S1, bm)->(
--    isBirationalMap(di, ideal S1, bm)
--    );

--isBirationalMap(Ring,Ideal,BasicList) := (R1, im, bm)->(
--    isBirationalMap(ideal R1, im, bm)
--    );

inverseOfMap(RingMap) :=(f)->(
    inverseOfMap(target f, source f, first entries matrix f)
    );
    
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mapOntoImage = method(); --given a map f : X -> Y, this creates the map f : X -> f(X).

mapOntoImage(RingMap) := (f)->(
        newMap := map(target f, ambient source f, matrix f);
        kk := ker(newMap);
        map(target newMap, (source newMap)/kk, matrix newMap)
        
);

mapOntoImage(Ring, Ring, BasicList) := (R,S,l)->(
        newMap := map(R, ambient S, l);
        mapOntoImage(newMap)
);
    
mapOntoImage(Ideal, Ideal, BasicList) := (a,b,l)->(
        newMap := map((ring a)/a, (ring b)/b, l);
        mapOntoImage(newMap)
);
    
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
isEmbedding = method(); --checks whether a map is a closed embedding.

isEmbedding(Ideal, Ideal, BasicList) := (a1, b1, f1)->(
        newMap := map((ring a1)/a1, (ring b1)/b1, f1);
        isEmbedding(newMap)
);

isEmbedding(Ring, Ring, BasicList) := (R1, S1, f1)->(
        newMap:=map(R1,S1,f1);
        isEmbedding(newMap)
);

isEmbedding(RingMap) := (f1)->(
        f2 := mapOntoImage(f1);
        flag := true;
        try ( h := inverseOfMap(f2) ) then (
                flag = isRegularMap(f2);
                if (flag == true) then(
                        flag = isRegularMap(h);
                );                
        )
        else (
                flag = false
        );
        flag
        
);
    
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 

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
		isBirationalMap
		(isBirationalMap, Ideal, Ideal, BasicList)
		(isBirationalMap, Ring, Ring, BasicList)
		(isBirationalMap, RingMap)
    Headline
        Checks if a map X -> Y between projective varieties is birational.
    Usage
		val = isBirationalMap(a,b,f)
		val = isBirationalMap(R,S,f)
		val = isBirationalMap(Pi)
	Inputs
		a:Ideal
			defining equations for X			
		b:Ideal
			defining equations for Y
		f:BasicList
            A list of where to send the variables in the ring of b, to in the ring of a.
        R:Ring
            the homogeneous coordinate ring of X
        S:Ring
            the homogeneous coordinate ring of Y
        Pi:RingMap
            A ring map S to R corresponding to X mapping to Y
    Outputs
        val:Boolean
            true if the map is birational, false if otherwise
    Description
        Text   
            This checks if a map between projective varieties is birational.  There are a number of ways to call this.  A simple one is to have a map between two graded rings.  In this case, the variables should be sent to elements of a single fixed degree.  Let's check that the plane quadratic cremona transformation is birational.
        Example
            R=QQ[x,y,z];
            S=QQ[a,b,c];
            Pi = map(R, S, {x*y, x*z, y*z});
            isBirationalMap(Pi)
        Text   
            We can also verify that a cover of $P^1$ by an elliptic curve is not birational.
        Example
            R=QQ[x,y,z]/(x^3+y^3-z^3);
            S=QQ[s,t];
            Pi = map(R, S, {x, y-z});
            isBirationalMap(Pi)
///                     

doc ///
	Key 
		idealOfImageOfMap
		(idealOfImageOfMap,Ideal,Ideal,Matrix)
		(idealOfImageOfMap,Ideal,Ideal,BasicList)
		(idealOfImageOfMap,Ring,Ring,Matrix)
		(idealOfImageOfMap,Ring,Ring,BasicList)
		(idealOfImageOfMap,RingMap)
	Headline
		Finds defining equations for the image of a rational map between varieties or schemes
	Usage
		im = idealOfImageOfMap(a,b,f)
		im = idealOfImageOfMap(a,b,g)
		im = idealOfImageOfMap(R,S,f)
		im = idealOfImageOfMap(R,S,g)
		im = idealOfImageOfMap(p)
	Inputs
		a:Ideal
			defining equations for X
		b:Ideal
			defining equations for Y
		f:Matrix
                        projective rational map given by polynomial representatives
		g:BasicList
			projective rational map given by polynomial representatives
		R:Ring
			coordinate ring of X
		S:Ring
			coordinate ring of Y
		p:RingMap
			projective rational map given by polynomial representatives
	Outputs
		im:Ideal
			defining equations for the image of f
	Description
		Text
			Given $f : X->Y \subset P^N$, this returns the defining ideal of $f(x) \subseteq P^N$. It should be noted for inputs that all rings are quotients of polynomial rings, and all ideals and ring maps are of these.  In particular, this function returns an ideal defining a subset of the  the ambient projective space of the image.  In the following example we consider the image of $P^1$ inside $P^1 \times P^1$.
		Example
			S = QQ[x,y,z,w];
			b = ideal(x*y-z*w);
			R = QQ[u,v];
			a = ideal(sub(0,R));
			f = matrix {{u,0,v,0}};
			idealOfImageOfMap(a,b,f)
///
			
doc ///
        Key
                dimImage
		(dimImage,Ideal,Ideal,Matrix)
		(dimImage,Ideal,Ideal,BasicList)
		(dimImage,Ring,Ring,Matrix)
		(dimImage,Ring,Ring,BasicList)
		(dimImage,RingMap)
        Headline
                Computes dimension of image of rational map of projective varieties
        Usage
                d = dimImage(a,b,f)
		d = dimImage(a,b,g)
		d = dimImage(R,S,f)
		d = dimImage(R,S,g)
		d = dimImage(p)
        Inputs 
                a: Ideal
                        defining equations for X
                b: Ideal
                        defining equations for Y
                f:Matrix
                        projective rational map given by polynomial represenative
        	g:BasicList
                        projective rational map given by polynomial representatives
                R:Ring
                        coordinate ring of X
                S:Ring
                        coordinate ring of Y
                p:RingMap
                        projective rational map given by polynomial representatives
	Outputs
                d:ZZ
			dimension of image
	Description
                Text
                        Gives the dimension of the image of a rational map. It should be noted for inputs that all rings are quotients of polynomial rings, and all ideals and ring maps are of these
                Example
                        S = QQ[x,y,z]
                        a = ideal(x^2+y^2+z^2)
                        T = QQ[u,v]
                        b = ideal(u^2+v^2)
                        f = matrix{{x*y,y*z}}
                        dimImage(a,b,f)
///


doc ///
        Key
                mapOntoImage
                (mapOntoImage, RingMap)
                (mapOntoImage, Ideal, Ideal, BasicList)
                (mapOntoImage, Ring, Ring, BasicList)
        Headline
                Given a map of rings, correspoing to $f : X -> Y$, this returns the map of rings corresponding to $f : X -> f(X)$.
        Usage
                h = mapOntoImage(f)
                h = mapOntoImage(a,b,l)
                h = mapOntoImage(R,S,l)                
        Inputs
                a:Ideal
                        defining equations for X
                b:Ideal
                        defining equations for Y
		l:BasicList
                        projective rational map given by polynomial represenatives of the same degree
                f:RingMap
                        the ring map corresponding to $f : X -> Y$
                R:Ring
                        coordinate ring for X
                S:Ring
                        coordinate ring for Y
                
        Outputs
                h:RingMap
			the map of rings corresponding to $f : X -> f(X)$.  
	Description
	        Text
	                This function is really simple, given $S -> R$, this just returns $S/kernel -> R$.  
	        Example 
	                R = QQ[x,y];
	                S = QQ[a,b,c];
	                f = map(R, S, {x^2, x*y, y^2});
	                mapOntoImage(f)
	                mapOntoImage(R,S,{x^2,x*y,y^2})
///

doc ///
        Key
                isEmbedding
                (isEmbedding, RingMap)
                (isEmbedding, Ideal, Ideal, BasicList)
                (isEmbedding, Ring, Ring, BasicList)
        Headline
                Given a map of rings, correspoing to $f : X -> Y$, this determines if this map embeds $X$ as a closed subscheme into $Y$.
        Usage
                val = isEmbedding(f)
                val = isEmbedding(a,b,l)
                val = isEmbedding(R,S,l)                
        Inputs
                a:Ideal
                        defining equations for X
                b:Ideal
                        defining equations for Y
		l:BasicList
                        projective rational map given by polynomial represenatives of the same degree
                f:RingMap
                        the ring map corresponding to $f : X -> Y$
                R:Ring
                        coordinate ring for X
                S:Ring
                        coordinate ring for Y
                
        Outputs
                val:Boolean
			true if the map is an embedding, otherwise false.
	Description
	        Text
	                Consider the Veronese embedding.
	        Example 
	                R = QQ[x,y];
	                S = QQ[a,b,c];
	                f = map(R, S, {x^2, x*y, y^2});
	                isEmbedding(f)
	        Text
	                Now consider the projection from a point on the plane to the line at infinity.
	        Example
	                R=QQ[x,y,z];
	                S=QQ[a,b];
	                f=map(R, S, {y,z});
	                isEmbedding(f)
	        Text 
	                That is obviously not an embedding.  It is even not an embedding when we restrict to a quadratic curve, even though it is a regular map.
	        Example
	                R=QQ[x,y,z]/(x^2+y^2-z^2);
	                S=QQ[a,b];
	                f=map(R,S, {y,z});
	                isRegularMap(f)
	                isEmbedding(f)
///



doc ///
    Key
        baseLocusOfMap
        (baseLocusOfMap, Matrix)
        (baseLocusOfMap, List)
        (baseLocusOfMap, RingMap)
    Headline
        Computes base locus of a map from a projective variety to projective space
    Usage
        I = baseLocusOfMap(M)
        I = baseLocusOfMap(L)
        I = baseLocusOfMap(h)
    Inputs
        M: Matrix
            Row matrix whose entries correspond to the coordinates of your map to projective space.
        L: List
            A list whose entries correspond to the coordinates of your map to projective space.
        h: RingMap
            A ring map corresponding to a map of projective varieties.
    Outputs
        I: Ideal
            The saturated defining ideal of the baselocus of the corresponding maps.
    Description
        Text
            This defines the locus where a given map of projective varieties is not defined.  For instance, consider the following rational map from $P^2$ to $P^1$
        Example
            R = QQ[x,y,z];
            S = QQ[a,b];
            f = map(R, S, {x,y});
            baseLocusOfMap(f)
        Text
            Observe it is not defined at the point [0:0:1], which is exactly what one expects.  However, we can restrict the map to a curve on $P^2$ and then it will be defined everywhere.
        Example
            R=QQ[x,y,z]/(y^2*z-x*(x-z)*(x+z));
            S=QQ[a,b];
            f=map(R,S,{x,y});
            baseLocusOfMap(f)
        Text
            Let us next consider the quadratic Cremona transformation.
        Example
            R=QQ[x,y,z];
            S=QQ[a,b,c];
            f=map(R,S,{y*z,x*z,x*y});
            J=baseLocusOfMap(f)
            minimalPrimes J
        Text
            The base locus is exactly the three points one expects.
///

doc ///
    Key
        isRegularMap
    Headline
        Checks whether a map to projective space is regular
    Usage
        b = isRegularMap(M)
    Inputs
        M: Matrix
            Row matrix whose entries correspond to the coordinates of your map to projective space
    Outputs
        b: Boolean
    Description
        Text
            This function just runs baseLocusOfMap(M) and checks if the ideal defining the base locus is the whole ring.
///  

doc ///
    Key
        inverseOfMap
		(inverseOfMap, Ideal, Ideal, BasicList)
		(inverseOfMap, Ring, Ring, BasicList)
		(inverseOfMap, RingMap)
    Headline
        Computes the inverse map of a given birational map between projective varieties. Returns an error if the map is not birational
    Usage
        f = inverseOfMap(I, J, L)
        f = inverseOfMap(R, S, L)
        f = inverseOfMap(g)
    Inputs
        I: Ideal
            Defining ideal of source
        J: Ideal
            Defining ideal of target
        L: List
            List of polynomials that define the coordinates of your birational map
        g: RingMap
            Your birational map
    Outputs
        f: RingMap
            Inverse function of your birational map
    Description
        Text
            Finds the inverse function of your birational map
    Caveat
        Only works for irreducible varieties right now
        
///

--******************************************
--******************************************
--******TESTS TESTS TESTS TESTS TESTS*******
--******************************************
--******************************************

TEST ///
	------------------------------------
	------- Tests for idealOfImageOfMap -------
	------------------------------------   
--Karl: I have no idea why this kernel should be (u^4, u*v^3)
	S = QQ[x,y,z]
        a = ideal(x^2+y^2+z^2)
        T = QQ[u,v]
        b = ideal(u^2+v^2)
        f = matrix{{x*y,y*z}}
        im = idealOfImageOfMap(a,b,f)  
	assert(im == ideal(v^4,u*v^3))
///

TEST ///
	S = QQ[x0,x1]
	a = ideal(0*x0)
	T = QQ[y0,y1,y2]
	b = ideal(0*y1)
	f = matrix{{x0^4,x0^2*x1^2,x1^4}}
	im = idealOfImageOfMap(a,b,f)
	assert(im == ideal(y2^2-y1*y^3))
///

TEST ///
	-- Since in Projective Space, check to make sure different representations give the same result
	S = QQ[x,y]
	a = ideal(0*x)
	T = QQ[u,v]
	b = ideal(0*v)
	f1 = matrix{{x,y}}
	f2 = matrix{{x^3*y^2,x^2*y^3}}
	assert(idealOfImageOfMap(a,b,f1)==idealOfImageOfMap(a,b,f2))
///

TEST ///
	-------------------------------------
	------ Tests for dimImage -----------
	-------------------------------------

	S = QQ[x,y,z]
        a = ideal(x^2+y^2+z^2)
        T = QQ[u,v]
        b = ideal(u^2+v^2)
        f = matrix{{x*y,y*z}}
        d = dimImage(a,b,f)
        assert(d == -1)
///

TEST ///        
        S = QQ[x0,x1]
        a = ideal(0*x0)
        T = QQ[y0,y1,y2]
        b = ideal(0*y1)
        f = matrix{{x0^4,x0^2*x1^2,x1^4}}
        d = dimImage(a,b,f)
        assert(d == 1)	
///

TEST ///
    -- Since in Projective Space, check to make sure different representations give the same result
    S = QQ[x,y]
    a = ideal(0*x)
    T = QQ[u,v]
    b = ideal(0*v)
    f1 = matrix{{x,y}}
    f2 = matrix{{x^3*y^2,x^2*y^3}}
    assert(dimImage(a,b,f1)==dimImage(a,b,f2))
///


	-------------------------------------
	-- Tests for baseLocusOfMap ---------
	-------------------------------------
TEST ///	
    R = QQ[x,y,z]	
	M = matrix{{x^2*y, x^2*z, x*y*z}}
	I = ideal(x*y, y*z, x*z)
	assert(I == baseLocusOfMap(M))
///

TEST ///	
    R = QQ[x,y,z]	
	L = {x^2*y, x^2*z, x*y*z}
	I = ideal(x*y, y*z, x*z)
	assert(I == baseLocusOfMap(L))
///

TEST ///
	-- reducible source
	R = QQ[x,y,z]/(x*y)
	M = matrix{{x^2, x*y, y^2}}
	I = ideal(x,y)
	assert(I == baseLocusOfMap(M))
///


	-------------------------------------
	----- isRegularMap -----------------
	-------------------------------------
TEST ///
	R = QQ[x,y,z,w]/(x*y - z*w)
	M = matrix{{1, 0, 0}}
	assert(isRegularMap(M))
///

TEST ///
    R = QQ[x,y]/(x*y)
    M = matrix{{x,y}}
    assert(isRegularMap(M))
///

TEST ///
    R = QQ[x,y,z]/(x^3 + y^3 - z^3)
    M = matrix{{y-z, x}}
    assert(isRegularMap(M) == false)
///


	-------------------------------------
	----- inverseOfMap  -----------------
	-------------------------------------



    
 
----FUTURE PLANS------

