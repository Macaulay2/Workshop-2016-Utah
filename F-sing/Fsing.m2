newPackage( "Fsing",
Version => "0.1", 
Date => "May 8th, 2016", 
Authors => {
     {Name => "Erin Bela",
     Email=> "ebela@nd.edu"
     },
     {Name => "DJ Bruce",
     Email => "djbruce@math.wisc.edu",
     HomePage => "http://www.math.wisc.edu/~djbruce/"
     },
     {Name => "Daniel Hernandez",
     Email => "dhernan@math.utah.edu",
     HomePage => "http://math.utah.edu/~dhernan/"
     },
     {Name => "Zhibek Kadyrsizova",
     Email => "zhikadyr@umich.edu"
     },
     {Name => "Mordechai Katzman",
     Email=> "m.katzman@sheffield.ac.uk",
     HomePage=> "http://www.katzman.staff.shef.ac.uk/"
     },
     {Name => "Sara Malec",
     Email=> "smalec@gsu.edu"
     },
     {Name => "Karl Schwede",
     Email => "schwede@math.psu.edu",
     HomePage => "http://math.utah.edu/~schwede/"
     },
     {Name => "Pedro Teixeira",
     Email => "pteixeir@knox.edu",
     HomePage => "http://www.knox.edu/academics/faculty/teixeira-pedro.html"
     },
     {Name=> "Emily Witt",
     Email=> "ewitt@umn.edu",
     HomePage => "http://math.umn.edu/~ewitt/"
     }
},
Headline => "A package for calculations of singularities in positive characteristic", 
DebuggingMode => true, 
Reload => true 
)


--- *** I SORTED THESE ALPHABETICALLY TO FIND MY WAY AROUND *** MK
export{
--IntegerComputations (IntegerComps.m2)
	"floorlog",
	"multOrder",
--ethRootFunctions (EthRoots.m2)
    "ascendIdeal", 
    "ascendIdealSafe",
    "ascendIdealSafeList",
    "AscentCount",
	"basePExpMaxE",
    "BinomialCheck",
    "binomialFPT",
    "canVector",
    "DiagonalCheck", 
    "diagonalFPT",
	"dividFraction",
    "estFPT",
    "ethRoot",
    "ethRootSafe", 	       
    "ethRootSafeList",    
    "EthRootStrategy",    
    "fancyEthRoot",
    "MonomialBasis",	
    "Substitution",
    --F-thresholds computations (FThresholds.m2)
    "factorList",    
    "FinalCheck",    
    "findAllCompatibleIdeals", ---MK	
    "findCPBelow",
	"findTestElementAmbient", ---Karl
	"FFiniteSupport", ---MK
	"findGeneratingMorphisms", ---MK
	"frobeniusPower",
    "FPTApproxList",     
    "FPT2VarHomog",     
    "FPT2VarHomogInternal",
	"FPureIdeals",
   "FTApproxList",
    "FTHatApproxList", 
	"FullMap", ---Karl
	"generatingMorphism", ---MK
	"generatingRoot", ---MK
    "guessFPT",
    "isBinomial",
    "isCP",
    "isDiagonal",
    "isFJumpingNumberPoly",
    "isFPTPoly",
    "isInLowerRegion",
    "isInUpperRegion",
    "MaxExp",
    "minimalCompatible",
    "MultiThread",    
    "Nontrivial",    
    "nu",
    "nuAlt",
    "NuCheck",
    "nuHat",
    "nuHatList",
    "nuList",
    "nuListAlt",
    "nuListAlt1",
    "Origin",
    "OutputRange",
    "PrintCP",
    "setFTData",
    "splittingField",
    "taxicabNorm",
    "paraTestModule", ---MK
    "paraTestModuleAmbient" ---MK
}

--*************************************************
--*************************************************
--This is the revised (and cleaned up) version
--of the PosChar package, which has been under 
--continuous development since the Wake Forest 
--Macaulay2 workshop of August 2012.
--Only well documented and working functions are 
--migrated to this package.
--*************************************************
--*************************************************

load "./Fsing/BasicFunctions.m2"

load "./Fsing/EthRoots.m2"

load "./Fsing/generatingMorphism.m2"


load "./Fsing/compatiblySplit.m2"


load "./Fsing/FPure.m2"

load "./Fsing/FFiniteSupport.m2"



load "./Fsing/parameterTestIdeal.m2"


load "./Fsing/FThresholds.m2"

load "./Fsing/SpecialFThresholds.m2"



load "./Fsing/FsingDocs.m2"

beginDocumentation()

load "./Fsing/EthRootsDoc.m2"

load "./Fsing/FThresholdsDoc.m2"

load "./Fsing/SpecialFThresholdsDoc.m2"

load "./Fsing/compatiblySplitDoc.m2"

load "./Fsing/FFiniteSupportDoc.m2"

load "./Fsing/generatingMorphismDoc.m2"

load "./Fsing/parameterTestIdealDoc.m2"

load "./Fsing/Fpure.m2"

