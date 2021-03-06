newPackage(
	"PolymakeInterface",
    	Version => "0.40", 
    	Date => "May 10, 2016",
    	Authors => {{Name => "Josephine Yu", 
		     Email => "josephine.yu@math.gatech.edu", 
		     HomePage => "http://people.math.gatech.edu/~jyu67"},
	            {Name => "Nathan Ilten",
	             HomePage => "http://math.berkeley.edu/~nilten",
	             Email => "nilten@math.berkeley.edu"},
	            {Name => "Qingchun Ren", 
		    Email => "qingchun.ren@gmail.com", 
		    HomePage => "http://math.berkeley.edu/~qingchun/"},
		    {Name => "David Swinarski",
		     Email => "dswinarski@fordham.edu",
		     HomePage => "http://faculty.fordham.edu/dswinarski/"},
		    {Name => "Madeline Brandt",
		     Email => "brandtm@berkeley.edu",
		     HomePage => "http://math.berkeley.edu/~brandtm/"}
		},
    	Headline => "a package for interfacing with polymake",
	PackageExports => {"PolyhedralObjects"},
    	DebuggingMode => true
    	)

export {
     "runPolymake",
     "ParseAllProperties",
     "BinaryOperationOutputType",
     "hasProperty",
     "getProperty",
     "getPropertyNames",
     "parseAllAvailableProperties"
     }


{*
Major changes from version 0.3 (Wake Forest workshop 2012) 
to 0.4 (Utah workshop 2016):

- We use PackageExports to load PolyhedralObjects.  This removes 
the need for several needsPackage("PolyhedralObjects") commands

- We changed the exports from symbols to strings as required in 
more recent versions of Macaulay2

- Polymake seems to have changed how it handles lattice points 
in order to accommodate unbounded polyhedra.  See the new polymake 
command LATTICE_POINTS_GENERATORS and M2 command LatticePointsGenerators 

- We added a function runPolymake(PolyhedralObject,PolyhedralObject,String) 
that allows a user to run a polymake binary operation, e.g. 
minkowski_sum($P,$Q).  We also added an option "BinaryOperationOutputType" 
so that the user can tell the runPolymake function whether to expect a 
Boolean (e.g. equal_polyhedra($P,$Q)) or a PolyhedralObject (e.g. 
minkowski_sum($P,$Q))

- We added tests.  

*}

---------------------------------------------------------------------------
-- Code
---------------------------------------------------------------------------

runPolymakePrefix = "polymake"

-- May 6, 2016: polymake 3.0 on some Macs returns extra 
-- characters with the output, ending in a bell (ascii 7).  
runPolymake = method(Options => {"ParseAllProperties" => false,"BinaryOperationOutputType"=>PolyhedralObject})
runPolymake(String) := o -> (script) -> (
     filename := temporaryFileName()|currentTime()|".poly";
     filename << script << endl << close;
     s:=get("!"|runPolymakePrefix|" "|filename);
     n:=regex(ascii(7),s);
     if n === null then return s;
     if #n > 1 then error "Parsing error; more than one bell";
     return substring(s,n_0_0+1,#s-1)
)

------------------------------------------------------------------------------
--Types of properties in polymake (hard coded)
------------------------------------------------------------------------------

emptyMatrixWithSource = (sourceDimensionPropertyName) -> (
     (P) -> (
          sourceDimension := runPolymake(P,sourceDimensionPropertyName);
          map(QQ^0,QQ^sourceDimension,0)
	  )
     )
propertyTypes = {
     {    
	  "M2PropertyName" => "AffineHull",
	  "PolymakePropertyName" => "AFFINE_HULL",
	  "ValueType" => "Matrix",
	  "EmptyMatrixFallback" => emptyMatrixWithSource("ConeAmbientDim")
	  },   
     {    
	  "M2PropertyName" => "AmbientDim",
	  "PolymakePropertyName" => "AMBIENT_DIM",
	  "ValueType" => "Integer"
	  }, 
    {    
	  "M2PropertyName" => "BoundaryLatticePoints",
	  "PolymakePropertyName" => "BOUNDARY_LATTICE_POINTS",
	  "ValueType" => "Matrix"
	  },
     {    
	  "M2PropertyName" => "Bounded",
	  "PolymakePropertyName" => "BOUNDED",
	  "ValueType" => "Boolean"
	  },
     {    
	  "M2PropertyName" => "ConeAmbientDim",
	  "PolymakePropertyName" => "CONE_AMBIENT_DIM",
	  "ValueType" => "Integer"
	  },
     {    
	  "M2PropertyName" => "ConeDim",
	  "PolymakePropertyName" => "CONE_DIM",
	  "ValueType" => "Integer"
	  },    
     {    
	  "M2PropertyName" => "EhrhartPolynomialCoeff",
	  "PolymakePropertyName" => "EHRHART_POLYNOMIAL_COEFF",
	  "ValueType" => "Array"
	  },
     {   
	  "M2PropertyName" => "Equations",
	  "PolymakePropertyName" => "EQUATIONS",
	  "ValueType" => "Matrix"
	  }, 
     {   
	  "M2PropertyName" => "Facets",
	  "PolymakePropertyName" => "FACETS",
	  "ValueType" => "Matrix"
	  },
     {    
	  "M2PropertyName" => "Feasible",
	  "PolymakePropertyName" => "FEASIBLE",
	  "ValueType" => "Boolean"
	  },
     {    
	  "M2PropertyName" => "FVector",
	  "PolymakePropertyName" => "F_VECTOR",
	  "ValueType" => "Vector"
	  },
     {    
	  "M2PropertyName" => "HilbertBasis",
	  "PolymakePropertyName" => "HILBERT_BASIS",
	  "ValueType" => "Matrix"
	  },  
     {   
	  "M2PropertyName" => "Inequalities",
	  "PolymakePropertyName" => "INEQUALITIES",
	  "ValueType" => "Matrix"
	  },
     {   
	  "M2PropertyName" => "InputLineality",
	  "PolymakePropertyName" => "INPUT_LINEALITY",
	  "ValueType" => "Matrix",
	  "EmptyMatrixFallback" => emptyMatrixWithSource("ConeAmbientDim")
	  },
     {    
	  "M2PropertyName" => "InputRays",
	  "PolymakePropertyName" => "INPUT_RAYS",
	  "ValueType" => "Matrix"
	  },
     {    
	  "M2PropertyName" => "InteriorLatticePoints",
	  "PolymakePropertyName" => "INTERIOR_LATTICE_POINTS",
	  "ValueType" => "Matrix"
	  },
     {    
	  "M2PropertyName" => "LatticePoints",
	  "PolymakePropertyName" => "LATTICE_POINTS",
	  "ValueType" => "Matrix"
	  },      
     {    
	  "M2PropertyName" => "LatticePointsGenerators",
	  "PolymakePropertyName" => "LATTICE_POINTS_GENERATORS",
	  "ValueType" => "LatticePointsGenerators"
	  },
     {    
	  "M2PropertyName" => "LatticeVolume",
	  "PolymakePropertyName" => "LATTICE_VOLUME",
	  "ValueType" => "Integer"
	  },     
     {   
	  "M2PropertyName" => "LinealitySpace",
	  "PolymakePropertyName" => "LINEALITY_SPACE",
	  "ValueType" => "Matrix",
	  "EmptyMatrixFallback" => emptyMatrixWithSource("ConeAmbientDim")	  },     
     {    
	  "M2PropertyName" => "LinearSpan",
	  "PolymakePropertyName" => "LINEAR_SPAN",
	  "ValueType" => "Matrix",
	  "EmptyMatrixFallback" => emptyMatrixWithSource("ConeAmbientDim")
	  },
     {   
	  "M2PropertyName" => "Points",
	  "PolymakePropertyName" => "POINTS",
	  "ValueType" => "Matrix"
	  },
     {    
	  "M2PropertyName" => "Rays",
	  "PolymakePropertyName" => "RAYS",
	  "ValueType" => "Matrix"
	  },                  
     {    
	  "M2PropertyName" => "Vertices",
	  "PolymakePropertyName" => "VERTICES",
	  "ValueType" => "Matrix"
	  }
}
            

propertyTypes = apply(propertyTypes, x -> new HashTable from x);
polymakePropertyNameToValueType = new HashTable from apply(propertyTypes, x -> (x#"PolymakePropertyName",x#"ValueType"));
M2PropertyNameToPolymakePropertyName = new HashTable from apply(propertyTypes, x -> (x#"M2PropertyName",x#"PolymakePropertyName"));
polymakePropertyNameToM2PropertyName = new HashTable from apply(propertyTypes, x -> (x#"PolymakePropertyName",x#"M2PropertyName"));
polymakePropertyNameToEmptyMatrixFallback = new HashTable from apply(select(propertyTypes,x -> (x#?"EmptyMatrixFallback")), x -> (x#"PolymakePropertyName",x#"EmptyMatrixFallback"));

---------------------------------------------------------------------------
----- Utilities

makeList = (str) -> (
     t := separateRegexp(" +",str);
     t = apply(t,value);
     select(t, x-> x =!= null)
     )

makeMatrix = (str) -> (
     L := lines str;
     L = select (L, x -> x!="");
     promote(matrix(L/makeList),QQ)
     )

------------------------------------------------------------------------------
--Parse properties into M2 format
------------------------------------------------------------------------------

parseUnknownProperty = method(TypicalValue => String)
parseUnknownProperty(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 print $object->"|propertyName|";";
     runPolymake(script)
     )

parseBooleanProperty = method(TypicalValue => Boolean)
parseBooleanProperty(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 my $v = $object->"|propertyName|";
	 if($v){print(\"true\")}else{print(\"false\");}";
     (runPolymake(script)=="true")
     )

parseScalarProperty = method()
parseScalarProperty(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 print $object->"|propertyName|";";
     value(runPolymake(script))
     )

parseListProperty = method(TypicalValue => List)
parseListProperty(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 print $object->"|propertyName|";";
     makeList(runPolymake(script))
     )

parseMatrixProperty = method(TypicalValue => Matrix)
parseMatrixProperty(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 print $object->"|propertyName|";";
     result := runPolymake(script);
     isEmpty := #(select (lines result, x -> x!=""))==0;
     if (isEmpty and polymakePropertyNameToEmptyMatrixFallback#?propertyName) then (
	  (polymakePropertyNameToEmptyMatrixFallback#propertyName)(P)
	  )
     else (
          makeMatrix(result)
	  )
     )

-- The polymake function LATTICE_POINTS_GENERATORS returns a list 
-- of three matrices (possibly empty)
parseLatticePointsGenerators = method(TypicalValue => List)
parseLatticePointsGenerators(PolyhedralObject,String) := (P, propertyName) -> (
     script := "use application \"polytope\";
         my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	 print $object->"|propertyName|";";
     results := runPolymake(script);
     results = separate(">\n",results);
     results=delete("",results);
     results = apply(results, i -> substring(i,1,#i-1));  
     for result in results list (  
         if result=="" then emptyMatrixWithSource("ConeAmbientDim") else makeMatrix(result)
     )
)

parseProperty = method()
parseProperty(PolyhedralObject,String) := (P, propertyName) -> (
     valueType := polymakePropertyNameToValueType#propertyName;
     if (valueType=="Boolean") then (
	  parseBooleanProperty(P,propertyName)
	  )
     else if (valueType=="Integer" or valueType=="Float" or valueType=="Scalar") then (
	  parseScalarProperty(P,propertyName)
	  )
     else if (valueType=="Vector" or valueType=="Array") then (
	  parseListProperty(P,propertyName)
	  )
     else if (valueType=="Matrix") then (
	  parseMatrixProperty(P,propertyName)
	  )
     else if (valueType=="LatticePointsGenerators") then (
	  parseLatticePointsGenerators(P,propertyName)
     )
     else (
	  parseUnknownProperty(P,propertyName)
	  )
     )

------------------------------------------------------------------------------
--Getting properties from cache
------------------------------------------------------------------------------

hasParsedProperty = method()
hasParsedProperty(PolyhedralObject,String) := (P, propertyName) -> (
     P#?propertyName
     )

hasCachedProperty = method()
hasCachedProperty(PolyhedralObject,String) := (P, propertyName) -> (
     P#?cache and P#cache#?"AvailableProperties" and M2PropertyNameToPolymakePropertyName#?propertyName and P#cache#"AvailableProperties"#?(M2PropertyNameToPolymakePropertyName#propertyName)
     )

hasProperty = method()
hasProperty(PolyhedralObject,String) := (P, propertyName) -> (
     hasParsedProperty(P,propertyName) or hasCachedProperty(P,propertyName)
     )

getProperty = method()
getProperty(PolyhedralObject,String) := (P, propertyName) -> (
     if (hasParsedProperty(P,propertyName)) then (
	  P#propertyName
	  )
     else if (hasCachedProperty(P,propertyName)) then (
	  propertyValue := parseProperty(P,M2PropertyNameToPolymakePropertyName#propertyName);
	  P#propertyName = propertyValue;
	  propertyValue
	  )
     else (
	  error ("The polyhedral object does not have property "|propertyName);
	  )
     )

getParsedPropertyNames = method()
getParsedPropertyNames(PolyhedralObject) := (P) -> (
     new Set from select(keys(P),x->(class(x)===String))
     )

getCachedPropertyNames = method()
getCachedPropertyNames(PolyhedralObject) := (P) -> (
     if (P#?cache and P#cache#?"AvailableProperties") then (
	  polymakePropertyNames := select(keys(P#cache#"AvailableProperties"),x->(polymakePropertyNameToM2PropertyName#?x));
	  new Set from apply(polymakePropertyNames,x->polymakePropertyNameToM2PropertyName#x)
	  )
     else (
	  new Set from {}
	  )
     )

getPropertyNames = method()
getPropertyNames(PolyhedralObject) := (P) -> (
     getParsedPropertyNames(P)+getCachedPropertyNames(P)
     )

parseAllAvailableProperties = method()
parseAllAvailableProperties(PolyhedralObject) := (P) -> (
     for propertyName in toList(getCachedPropertyNames(P)-getParsedPropertyNames(P)) do (
	  polymakePropertyName := M2PropertyNameToPolymakePropertyName#propertyName;
	  P#propertyName = parseProperty(P,polymakePropertyName);
	  )
     )

---------------------------------------------------------------------------
-- Create a polymake input file

toPolymakeFormat = method(TypicalValue => String)
toPolymakeFormat(String, Boolean) := (propertyName, x) -> (
     if x === null then ""
     else (
	  S := propertyName|"\n";
	  if x then (
	       S = S|"1";
	       );
	  S
	  )
     )
toPolymakeFormat(String, Matrix) := (propertyName, M) -> (
     if M === null then ""
     else(
     	  S := propertyName|"\n";
     	  if numRows M > 0 then
	     S = S|replace("\\|", "", toString net M);
     	  S
     	  )
     )
toPolymakeFormat(String,Vector) := (propertyName,V) -> (
     if V === null then ""
     else(
     	  S := propertyName|"\n";
     	  if length V > 0 then
              S = S|replace("\\|", "", toString net matrix{V});     
     	  S
     	  )
     )
toPolymakeFormat(String,List) := (propertyName,L) -> (
     if L === null then ""
     else(
     	  S := propertyName|"\n";
     	  S = S|concatenate(apply(L,x->(toString(x)|" ")));
     	  S
     	  )
     )
toPolymakeFormat(String,ZZ) := (propertyName,x) -> (
     if x === null then ""
     else(
     	  S := propertyName|"\n"|x|"\n";
     	  S
     	  )
     )
toPolymakeFormat(PolyhedralObject) := (P) -> (
     parsedPropertyNames := select(toList(getParsedPropertyNames(P)), x->(M2PropertyNameToPolymakePropertyName#?x));
     parsedProperties := apply(parsedPropertyNames,x->(M2PropertyNameToPolymakePropertyName#x,P#x));
     concatenate apply(parsedProperties, a -> toPolymakeFormat(a)|"\n\n")
     )

writePolymakeFile = method(TypicalValue => Nothing)
writePolymakeFile(PolyhedralObject, String) := (P, header) ->(
    if P#?cache then (
        if P.cache#?"PolymakeFile" then return 	P.cache#"PolymakeFile"
    );
    fileName := temporaryFileName()|currentTime()|".poly";
    << "using temporary file " << fileName << endl;
    fileName << header << endl << endl << toPolymakeFormat(P) << endl << close;
    if (not(P#?cache)) then (
        P#cache = new MutableHashTable;
    );
    P#cache#"PolymakeFile" = fileName;
    fileName	  
)
writePolymakeFile(PolyhedralObject) := (P) -> (
     writePolymakeFile(P, "")
     )
writePolymakeFile(Cone) := (C) -> (
     header := "_type Cone\n";
     writePolymakeFile(C, header)
     )

---------------------------------------------------------------------------
----- Run Polymake with PolyhedralObject




runPolymake(PolyhedralObject,String) := o -> (P,propertyName) -> (
     if (not(M2PropertyNameToPolymakePropertyName#?propertyName)) then (
	  error ("Property does not exist:"|propertyName)
	  )
     else (
	  polymakePropertyName := M2PropertyNameToPolymakePropertyName#propertyName;
          if (not(hasProperty(P,propertyName))) then (
	       if (not(P#?cache)) then (
		    P#cache = new MutableHashTable;
		    );
               if (not(P#cache#?"PolymakeFile")) then (
                    P#cache#"PolymakeFile" = writePolymakeFile(P);
	            )
	       else if (not(isSubset(getParsedPropertyNames(P),getCachedPropertyNames(P)))) then (
		    parseAllAvailableProperties(P);
		    P#cache#"PolymakeFile" = writePolymakeFile(P);
		    );		
               script := "use application \"polytope\";
                    my $object = load(\""|(P#cache#"PolymakeFile")|"\");
	            $object->"|polymakePropertyName|";
	            save($object,\""|(P#cache#"PolymakeFile")|"\");
	            my @properties = $object->list_properties;
	            my $numberOfProperties = scalar @properties;
	            for(my $i=0;$i<$numberOfProperties;$i++)
	            {print \"$properties[$i]\\n\";}";
               propertiesString := runPolymake(script);
               propertiesList := lines propertiesString;
               P#cache#"AvailableProperties" = new Set from propertiesList;
	       if member(propertyName,{"AmbientDim","LinearSpan","LatticePoints"}) then ( 
		   if o#"ParseAllProperties" then parseAllAvailableProperties(P);
		   return runPolymakeExceptions(P,polymakePropertyName)
	       );
	  );
	  if (o#"ParseAllProperties") then (
	       parseAllAvailableProperties(P);
	       );
	  getProperty(P,propertyName)
	  )
     )
 
-- For some properties, it seems that polymake does not save them in the file 
-- If we get here: then the propertyName is in {"AmbientDim","LinearSpan","LatticePoints"}
-- P has a cache and a PolymakeFile already
-- Write a script to print the property.  Then parse the polymake output
runPolymakeExceptions = (P,polymakePropertyName) -> (
    if polymakePropertyName == "AMBIENT_DIM" then (
	return parseScalarProperty(P,"AMBIENT_DIM")
    );
    if polymakePropertyName == "LINEAR_SPAN" then (
	return parseMatrixProperty(P,"LINEAR_SPAN")
    );
    if polymakePropertyName == "LATTICE_POINTS" then (
	return parseMatrixProperty(P,"LATTICE_POINTS")
    );	           
); 
------------------------------------------------------------------------------
-- Run any binary operation in polymake on two PolyhedralObjects
------------------------------------------------------------------------------
 
runPolymake(PolyhedralObject,PolyhedralObject,String):= o ->(P,Q,binaryOperationName) -> (
    writePolymakeFile(P);
    writePolymakeFile(Q);   
    fn:=temporaryFileName()|currentTime()|".poly";
    script:="";
    if o#"BinaryOperationOutputType"===Boolean then (      
        script="use application \"polytope\";
                 my $P = load(\""|P.cache#"PolymakeFile"|"\");
                 my $Q = load(\""|Q.cache#"PolymakeFile"|"\");
	         my $v = "|binaryOperationName|"($P,$Q);
                 if($v){print(\"true\")}else{print(\"false\");}";
        return(runPolymake(script)=="true")
    );
    script="use application \"polytope\";
             my $P = load(\""|P.cache#"PolymakeFile"|"\");
             my $Q = load(\""|Q.cache#"PolymakeFile"|"\");		 
             my $R = "|binaryOperationName|"($P,$Q);
             $R->CONE_AMBIENT_DIM;
             my $fn=\""|fn|"\";
	     save($R,$fn);
             my @properties = $R->list_properties;
             my $numberOfProperties = scalar @properties;
             for(my $i=0;$i<$numberOfProperties;$i++)
             {print \"$properties[$i]\\n\";}";
    propertiesString := runPolymake(script);
    propertiesList := lines propertiesString;
    R:=new Polyhedron from {cache=>new MutableHashTable from {"PolymakeFile"=>fn,"AvailableProperties"=>new Set from propertiesList}};
    if (o#"ParseAllProperties") then (
        parseAllAvailableProperties(R);
    );
    R
); 


---------------------------------------------------------------------------
-----------------------DOCUMENTATION----------------------
---------------------------------------------------------------------------

beginDocumentation()

doc ///
  Key
    PolymakeInterface
  Headline
    Interface for Polymake
  Description
   Text 
     {\tt polymake} is a software for convex polyhedra, simplicial complexes, and other discrete geometric objects, written by Ewgenij Gawrilow and Michael Joswig.  It is available at @HREF "http://www.math.tu-berlin.de/polymake/"@. The user should have {\tt polymake} installed on their machine.

   Text 
     Warning: this package is still under development. It may not function properly on some machines. Use it at your own risk.

   Text
     The current polymake interface assumes that you can run the command {\tt polymake [script_file]} in the terminal to run a polymake script. To check if this works, just type {\tt polymake} in a terminal window. However, this does not work on some Mac machines. If you have this problem, see @HREF "https://github.com/Macaulay2/Workshop-2016-Utah/wiki/Running-polymake-on-a-Mac"@.

   Text
     We can use the interface to get properties of a @TO "PolyhedralObjects::PolyhedralObject"@.  Here is a list of @TO "Available properties"@.

   Example
     --needsPackage "PolyhedralObjects";
     P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
     runPolymake(P, "FVector")

   Text
     Instead of performing the computation every time, we store the result in a cache, and reuse the result every time the interface is called with the same object.

   Example
     --needsPackage "PolyhedralObjects";
     P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
     runPolymake(P, "FVector")
     runPolymake(P, "FVector")
   Text
     Look in the cache for the properties already computed.
   Example
     --needsPackage "PolyhedralObjects";
     P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
     runPolymake(P, "FVector")
     getPropertyNames(P)
     hasProperty(P, "FVector")
     getProperty(P, "FVector")
   Text
     We can save the {\tt Polymake} output as a boolean, an integer, a list, or a matrix, depending on the type of the output.
   Example
     --needsPackage "PolyhedralObjects";
     P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
     runPolymake(P, "Bounded")
     runPolymake(P, "ConeDim")
     runPolymake(P, "FVector")
     runPolymake(P, "Facets")
  Caveat
     You should not manually remove the temporary files created in {\tt /tmp} before all computation is complete.
     You should not change the polytope object essentially after calling the Polymake interface. You may get wrong results if the polytope object is not consistent with the Polymake cache.
  SeeAlso
     PolyhedralObjects
///


doc ///
  Key
    "Available properties"
  Headline
    Available properties
  Description
   Text 
     The following PolyhedralObject output properties are available and have been tested: 
   
     "AffineHull"
     
     "AmbientDim"
     
     "BoundaryLatticePoints"
     
     "Bounded"
     
     "ConeAmbientDim"
     
     "ConeDim"
     
     "EhrhartPolynomialCoeff"
     
     "Facets"
     
     "Feasible"
     
     "FVector"
     
     "HilbertBasis"
     
     "InteriorLatticePoints"

     "LatticePoints"
          
     "LatticePointsGenerators"
     
     "LatticeVolume"
     
     "LinealitySpace"
     
     "LinearSpan"
    
     "Rays"
     
     "Vertices"
     
     The following binary operations in polymake have been tested:

     "conv" (compute the convex hull of two polyhedra)
          
     "equal_polyhedra"
     
     "included_polyhedra"
     
     "intersection"
     
     "minkowski_sum"
     
     "product"
     
     For more binary operations, see the polymake documentation at @HREF "https://polymake.org/release_docs/3.0/polytope.html"@  However these are untested; use at your own risk.  
///

doc ///
    Key
        runPolymake
	(runPolymake,String)
	(runPolymake,PolyhedralObject,String)
	(runPolymake,PolyhedralObject,PolyhedralObject,String)
    Headline
        perform computation using Polymake
    Usage
        result = runPolymake(polymakeScript)
	result = runPolymake(P,propertyName)
	result = runPolymake(P,Q,binaryOperationName)
	result = runPolymake(P,propertyName,"ParseAllProperties"=>true)
    Inputs
        polymakeScript:String
	P:PolyhedralObject
	propertyName:String
    Outputs
        result:Thing
	    the result returned from Polymake
    Description
        Text
	    runPolymake(polymakeScript) runs a Polymake script and returns whatever Polymake prints in its standard output as a String.
	Example
	    script = "use application \"polytope\"; my $a = cube(2,2); print $a->F_VECTOR;";
	    runPolymake(script)
        Text
	    runPolymake(P,propertyName) runs Polymake on a PolyhedralObject (such as Polyhedron, Cone) P to get a property.
	Example
            --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector")
	Text
	    runPolymake(P,Q,binaryOperationName) runs a binary operation in polymake on the two input PolyhedralObjects
	Example
	    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}}
            Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,0,1}}}
	    R = runPolymake(P,Q,"minkowski_sum")
	Text
	    When Polymake computes a property of a polyhedral object, it may compute other properties in the process. If the "ParseAllProperties" option is set {\tt true}, then all properties that are computed in the process are parsed into M2 format. Otherwise, only the needed property is parsed, and the other properties are stored in a temporary file in Polymake format.
	Example
	    --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector","ParseAllProperties"=>true)
	    Q = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(Q, "FVector","ParseAllProperties"=>false)
	Text
	    The type of the output varies according to the type of the property.
	Example
            --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
	    runPolymake(P, "Feasible")
	    runPolymake(P, "ConeDim")
            runPolymake(P, "FVector")
	    runPolymake(P, "Facets")
///

doc ///
    Key
        hasProperty
    Headline
        check if a polyhedral object has a property in its Polymake cache
    Usage
        result = hasProperty(P,propertyName)
    Inputs
	P:PolyhedralObject
	propertyName:String
    Outputs
        result:Boolean
	    whether the polyhedral object has the property in its Polymake cache
    Description
        Text
	    Checks if a polyhedral object has a property in its Polymake cache.
	Example
            --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector")
	    hasProperty(P, "Facets")
///
-- May 6, 2016: the getProperty command below has not been working for a long time
doc ///
    Key
        getProperty
    Headline
        get the value of a property in the Polymake cache
    Usage
        result = getProperty(P,propertyName)
    Inputs
	P:PolyhedralObject
	propertyName:String
    Outputs
        result:Thing
	    the value of a property
    Description
        Text
	    Gets the value of a property in the Polymake cache.
	Example
            --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector")
--	    getProperty(P, "Facets")
///

doc ///
    Key
        getPropertyNames
    Headline
        get all property names in the Polymake cache
    Usage
        result = getProperty(P)
    Inputs
	P:PolyhedralObject
	propertyName:String
    Outputs
        result:Set
	    all property names in the Polymake cache
    Description
        Text
	    Gets all property names in the Polymake cache.
	Example
            --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector")
	    getPropertyNames(P)
///

doc ///
    Key
        parseAllAvailableProperties
    Headline
        parse values of all property in the Polymake cache to M2 format
    Usage
        parseAllAvailableProperties(P)
    Inputs
	P:PolyhedralObject
    Description
        Text
	    Parses values of all property in the Polymake cache to M2 format.
	Example
	    --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector")
	    parseAllAvailableProperties(P)
///

doc ///
    Key
        "ParseAllProperties"
    Headline
        optional parameter for runPolymake
    Description
        Text
	    When Polymake computes a property of a polyhedral object, it may compute other properties in the process. If the "ParseAllProperties" option is set {\tt true}, then all properties that are computed in the process are parsed into M2 format. Otherwise, only the needed property is parsed, and the other properties are stored in a temporary file in Polymake format.
	Example
	    --needsPackage "PolyhedralObjects";
            P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
            runPolymake(P, "FVector","ParseAllProperties"=>true)
///

---------------------------------------------------------------------------
------------------------- TESTS ---------------------------
---------------------------------------------------------------------------

-- Test runPolymake(script)

TEST ///
    result = runPolymake("print 123");
    assert(result=="123");
///

TEST ///
    str = runPolymake("use application \"polytope\"; my $a = cube(2,2); print $a->F_VECTOR;");
    t = separateRegexp(" +",str);
    t = apply(t,value);
    assert(t=={4,4});
///

-----------------------------------
-- Test runPolymake(P,propertyName)
-----------------------------------

-- The tests are in alphabetic order by propertyName

-- Test "AffineHull"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1}}};
    result = runPolymake(P,"AffineHull");
    assert(result==map((QQ)^1,(QQ)^4,{{-1, 1, 0, 0}}) );
    assert(class(result)===Matrix);
///

-- Test "AmbientDim"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"AmbientDim");
    assert(class(result)===ZZ);
    assert(result== 2);
///


-- Test "BoundaryLatticePoints"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,3},{1,3,0},{1,3,3}}};
    result = runPolymake(P,"BoundaryLatticePoints");
    assert(class(result)===Matrix);
    L1 = set entries substitute(result,ZZ);
    L2 = set {{1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 0, 3}, {1, 1, 0}, {1, 1, 3}, {1, 2, 0}, {1, 2, 3}, {1, 3, 0}, {1, 3, 1}, {1, 3, 2}, {1, 3, 3}};
    assert(L1 === L2);
///

-- Test "Bounded"
TEST ///
    P = new Polyhedron from {"Inequalities" => matrix{{0,1,0},{0,0,1}}};
    result = runPolymake(P,"Bounded");
    assert(class(result)===Boolean);
    assert(result === false);
    P = new Polyhedron from {"Inequalities" => matrix{{0,1,0},{0,0,1},{1,-1,0},{1,0,-1}}};
    result = runPolymake(P,"Bounded");
    assert(class(result)===Boolean);
    assert(result === true);
///

-- Test "ConeAmbientDim"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"ConeAmbientDim");
    assert(class(result)===ZZ);
    assert(result === 3);
///

-- Test "ConeDim"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"ConeDim");
    assert(class(result)===ZZ);
    assert(result === 3);
///

-- Test "EhrhartPolynomialCoeff"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"EhrhartPolynomialCoeff");
    assert(class(result)===List);
    assert(result === {1, 2, 1});
///

-- Test "Equations"
-- Note: EQUATIONS is an input property only
TEST ///
    P = new Polyhedron from {"Inequalities" => matrix{{0,1,0,0},{0,0,1,0},{1,-1,0,0},{1,0,-1,0}}, "Equations"=> matrix {{5,0,0,-1}}};
    result = runPolymake(P,"Vertices");
    assert(class(result)===Matrix);
    L1=set entries substitute(result,ZZ);
    L2=set {{1, 1, 0, 5}, {1, 1, 1, 5}, {1, 0, 1, 5}, {1, 0, 0, 5}};
    assert(L1 === L2)
///

-- Test "Facets"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"Facets");
    assert(class(result)===Matrix);
    L1 = set entries result;
    L2 = set { (1/1)*{0,1,0}, (1/1)*{0,0,1}, (1/1)*{1,-1,0}, (1/1)*{1,0,-1}};
    assert(L1 === L2);
///

-- Test "Feasible"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"Feasible");
    assert(class(result)===Boolean);
    assert(result === true);
    P = new Polyhedron from {"Inequalities" => matrix{{0,1,0},{0,0,1},{-1,-1,0}}};
    result = runPolymake(P,"Feasible");
    assert(result === false);
///

-- Test "FVector"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"FVector");
    assert(result=={4,4});
    assert(class(result)===List);
///

-- Test "HilbertBasis"
--TEST ///
--
--///


-- Test "Inequalities"
-- Note: INEQUALITIES is an input property only
TEST ///
    P = new Polyhedron from {"Inequalities" => matrix{{0,1,0},{0,0,1},{1,-1,0},{1,0,-1}}};
    result = runPolymake(P,"Vertices");
    assert(class(result)===Matrix);
    L1 = set entries substitute(result,ZZ);
    L2 = set {{1, 1, 0}, {1, 1, 1}, {1, 0, 1}, {1, 0, 0}};
    assert(L1 === L2)
///

-- Test "InputLineality"
--TEST ///
--
--///

-- Test "InputRays"
--TEST ///
--
--///

-- Test "InteriorLatticePoints"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,3},{1,3,0},{1,3,3}}};
    result = runPolymake(P,"InteriorLatticePoints");
    assert(class(result)===Matrix);
    L1 = set entries substitute(result,ZZ);
    L2 = set {{1, 1, 1}, {1, 1, 2}, {1, 2, 1}, {1, 2, 2}};
    assert(L1 === L2);
///

-- Test "LatticePoints"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,2},{1,2,0},{1,2,2}}};
    result = runPolymake(P,"LatticePoints");
    assert(class(result)===Matrix);
    L1 = set entries substitute(result,ZZ);
    L2 = set {{1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 2, 0}, {1, 2, 1}, {1, 2, 2}};
    assert(L1 === L2);
///

-- Test "LatticePointsGenerators"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,2},{1,2,0},{1,2,2}}};
    result=runPolymake(P,"LatticePointsGenerators");
    assert( class(result) === List);
    assert( #result === 3);
    assert( result_0 ==
matrix {{1/1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 2, 0}, {1, 2, 1}, {1, 2, 2}})
    assert( (result_1)(P) == map(QQ^0,QQ^3,0))    
    assert( (result_2)(P) == map(QQ^0,QQ^3,0))    
///

-- Test "LatticeVolume"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,3},{1,3,0},{1,3,3}}};
    result=runPolymake(P,"LatticeVolume");
    assert( class(result) === ZZ);
    assert( result === 18);
///

-- Test "LinealitySpace"
--TEST ///
--
--///


-- Test "LinearSpan"
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0,0},{1,1,0,0},{1,0,1,0},{1,1,1,0}}};
    result=runPolymake(P,"LinearSpan");
    assert(class(result)===Matrix);
    assert(result === map((QQ)^1,(QQ)^4,{{0, 0, 0, 1}}));
///


-- Test "Points" and Vertices
-- Note: Points is an input property only
TEST ///
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,2},{1,2,0},{1,2,2},{1,1,1}}};
    result=runPolymake(P,"Vertices");
    assert( class(result) === Matrix);
    L1=set entries substitute(result,ZZ);
    L2=set {{1,0,0},{1,0,2},{1,2,0},{1,2,2}};
    assert( L1===L2);
///

-- Test "Rays"
--TEST ///
--
--///

-- Test "Vertices"
-- See test of Points above

-------------------------
-- Test binary operations
-------------------------

-- Test "conv"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,0,1}}};
    result = runPolymake(P,Q,"conv","ParseAllProperties"=>true);
    assert(class(result)===Polyhedron);
    L1 = set entries substitute(result#"Points",ZZ);
    L2 = set {{1,0,0},{1,0,1},{1,1,0}};    
    assert(L1===L2);
///


-- Test "equal_polyhedra"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0},{1,2,0}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,2,0}}};
    result = runPolymake(P,Q,"equal_polyhedra","BinaryOperationOutputType"=>Boolean);
    assert(class(result)===Boolean);
    assert(result===true);
    R = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    result = runPolymake(P,R,"equal_polyhedra","BinaryOperationOutputType"=>Boolean);
    assert(class(result)===Boolean);
    assert(result===false);
///


-- Test "included_polyhedra"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,2,0}}};
    result = runPolymake(P,Q,"included_polyhedra","BinaryOperationOutputType"=>Boolean);
    assert(class(result)===Boolean);
    assert(result===true);
    R = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,-1,0}}};
    result = runPolymake(P,R,"equal_polyhedra","BinaryOperationOutputType"=>Boolean);
    assert(class(result)===Boolean);
    assert(result===false);
///

-- Test "intersection"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0},{1,0,1},{1,1,1}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,1,0},{1,1,1},{1,2,0},{1,2,1}}};
    result = runPolymake(P,Q,"intersection","ParseAllProperties"=>true);
    assert(class(result)===Polyhedron);
    result = runPolymake(result,"Vertices");
    L1 = set entries substitute(result,ZZ);
    L2 = set {{1,1,0},{1,1,1}};    
    assert(L1===L2);
///

-- Test "minkowski_sum"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,0,1}}};
    result = runPolymake(P,Q,"minkowski_sum","ParseAllProperties"=>true);
    L1 = set entries substitute(result#"Points",ZZ);
    L2 = set {{1,0,0},{1,0,1},{1,1,0},{1,1,1}};    
    assert(L1===L2);
///

-- Test "product"
TEST ///
    P = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    Q = new Polyhedron from {"Points"=> matrix {{1,0,0},{1,1,0}}};
    result = runPolymake(P,Q,"product","ParseAllProperties"=>true);
    L1 = set entries substitute(result#"Vertices",ZZ);
    L2 = set {{1, 0, 0, 0, 0}, {1, 0, 0, 1, 0}, {1, 1, 0, 0, 0}, {1, 1, 0, 1, 0}};    
    assert(L1===L2);
///


----------------------------
-- More interesting tests
----------------------------

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"Facets");
    result = runPolymake(P,"Feasible");
    assert(result);
    assert(class(result)===class(true));
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"Feasible");
    runPolymake(P,"Feasible");
    runPolymake(P,"Feasible");
    runPolymake(P,"FVector");
    runPolymake(P,"Facets");
    runPolymake(P,"FVector");
    result = runPolymake(P,"Feasible");
    assert(result);
    assert(class(result)===class(true));
///

-- Test hasProperty

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    assert(hasProperty(P,"Points"));
    assert(not(hasProperty(P,"POINTS")));
    assert(not(hasProperty(P,"FVector")));
    assert(not(hasProperty(P,"PolymakeFile")));
    assert(not(hasProperty(P,"AvailableProperties")));
    assert(not(hasProperty(P,"blahblahblahblahblahblahblahblahblahblahblahblahblah")));
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"FVector");
    assert(hasProperty(P,"Points"));
    assert(hasProperty(P,"FVector"));
    assert(not(hasProperty(P,"POINTS")));
    assert(not(hasProperty(P,"F_VECTOR")));
    assert(not(hasProperty(P,"PolymakeFile")));
    assert(not(hasProperty(P,"AvailableProperties")));
    assert(not(hasProperty(P,"blahblahblahblahblahblahblahblahblahblahblahblahblah")));
///


-- Test getProperty
TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = getProperty(P,"Points");
    assert(result==matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}});
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"FVector");
    resultPoints = getProperty(P,"Points");
    assert(resultPoints==matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}});
    resultFVector = getProperty(P,"FVector");
    assert(resultFVector=={4,4});
///

-- Test getPropertyNames

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = getPropertyNames(P);
    assert(result#?"Points");
    assert(not(result#?cache));
    assert(not(result#?"POINTS"));
    assert(not(result#?"FVector"));
    assert(not(result#?"PolymakeFile"));
    assert(not(result#?"blahblahblahblahblahblahblahblahblahblahblahblah"));
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"FVector");
    result = getPropertyNames(P);
    assert(result#?"Points");
    assert(result#?"FVector");
    assert(not(result#?cache));
    assert(not(result#?"POINTS"));
    assert(not(result#?"F_VECTOR"));
    assert(not(result#?"PolymakeFile"));
    assert(not(result#?"blahblahblahblahblahblahblahblahblahblahblahblah"));
///


-- Test "ParseAllProperties"
TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"FVector","ParseAllProperties"=>true);
    assert(result=={4,4});
    assert(class(result)===class({4,4}));
    propertyNames = getPropertyNames(P);
    for key in keys(propertyNames) do (
	assert(P#?key);
	);
///

-- Test "parseAllAvailableProperties"

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    result = runPolymake(P,"FVector");
    parseAllAvailableProperties(P);
    assert(result=={4,4});
    assert(class(result)===class({4,4}));
    propertyNames = getPropertyNames(P);
    for key in keys(propertyNames) do (
	assert(P#?key);
	);
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Feasible"=>true};
    result = runPolymake(P,"Feasible");
    assert(result);
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Feasible"=>false};
    result = runPolymake(P,"Feasible");
    assert(not result);
///

TEST ///
    --needsPackage "PolyhedralObjects";
    C = new Cone from {"InputRays" => matrix{{0,0,1},{0,1,1},{1,0,1},{1,1,1}}};
    result = runPolymake(C,"LinearSpan");
    assert(class(result)===Matrix);
    assert(numRows(result)==0);
    assert(numColumns(result)==3);
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"Feasible");
    P#"Facets" = matrix{{0,0,1},{0,1,0},{1,0,-1},{1,-1,0}};
    result = runPolymake(P,"FVector");
    assert (result=={4,4});
///



TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}};
    runPolymake(P,"FVector");
    result = runPolymake(P,"ConeDim");
    assert(class(result)===class(1));
///

TEST ///
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"Points" => matrix{{1,0,0},{1,0,1},{1,1,0},{1,1,1}},"FVector" => {4,4}};
    result = runPolymake(P,"ConeDim");
    assert(class(result)===class(1));
///

TEST ///
    inputRays = matrix {{1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 1,
       0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
       0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, {1, 0, 0, 0,
       1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
       1, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {1, 0, 1,
       0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 1, 0, 0}, {1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {1, 0,
       1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 0, 0, 1, 0}, {1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1,
       0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 0, 0, 0, 0, 1}};
    --needsPackage "PolyhedralObjects";
    P = new Polyhedron from {"InputRays"=>inputRays};
    runPolymake(P,"FVector");
    result = runPolymake(P,"ConeDim");
    assert(class(result)===class(1));
///


end

---------------------------------------------------------------------------
------------------------- END ---------------------------
---------------------------------------------------------------------------

--Failed examples

///
    --This is a bug in Polymake. Need to contact its developers.
    --needsPackage "PolyhedralObjects"
    P = new Polyhedron from {"ConeDim"=>1};
    --ConeDim is insufficient to determine Points
    runPolymake(P,"Points");
    --So should not have property Points after calling runPolymake
    assert(not(hasProperty(P,"Points")));
///

--For Emacs with F11 hotkey

restart
--needsPackage "PolymakeInterface"
installPackage "PolymakeInterface"
check PolymakeInterface
help PolymakeInterface
viewHelp PolymakeInterface

---------------------------------------------------------------------------
------------------------- TO DO ---------------------------
---------------------------------------------------------------------------

{* Version 0.3 To do list
-- There is 1 known bug
-- Need to handle error from Polymake
-- Add new properties to the table
-- Parse IncidenceMatrix,SimplicialComplex objects and other types
-- Parse all properties in one round to save the overhead
-- Documentation
*}

{* Version 0.4 To do list:
Add function to get faces of a specified dimension
Add function to get a triangulation of a polytope 

Many properties above are only tested on a bounded polytope.  
Where appropriate, add an example of an unbounded polytope to the test.  

The following properties for unbounded polyhedra are not tested at all yet:
  HilbertBasis
  InputLineality
  InputRays
  LinealitySpace

Add functions and tests for other kinds of PolyhedralObjects (Cone, Fan, etc.)

*}