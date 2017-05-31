newPackage( "Fsing",
Version => "0.1", 
Date => "May 30th, 2017", 
Authors => {
     {Name => "Erin Bela",
     Email=> "ebela@nd.edu"
     },
     {Name => "Alberto F. Boix",
     Email=> "alberto.fernandezb@upf.edu"
     },
     {Name => "David J. Bruce",
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
Reload => true,
AuxiliaryFiles=>true 
)


--- *** I SORTED THESE ALPHABETICALLY TO FIND MY WAY AROUND *** MK
--- *** Reorganized by subpackages, so we know where to find stuff *** PT

export{
--BasicFunctions (BasicFunctions.m2) 
    "carryTest",  
    "basePExp",    
    "digit", 	   
    "denom",   
    "divideFraction",
    "firstCarry", 
    "floorlog",
    "fracPart", 
    "getCanVector",
    "getNumAndDenom", 
    "maxIdeal", 
    "multOrder",
    "num",
    "taxicabNorm",
    "truncatedBasePExp",
    "NoZeroC", --option to force certain behavior from a function
    
--ethRootFunctions (EthRoots.m2)
    "ascendIdeal", 
    "AscentCount",
    "boundLargestCompatible", ---MK
    "ethRoot",
    "ethRootRingElements",   
    "EthRootStrategy",    
    "minimalCompatible",
    "MonomialBasis",	
    "Substitution",
    "getFieldGenRoot",
    
--Frobenius Powers (frobeniusPowers.m2)
    "fastExp",
    "frobeniusPower",
    "genFrobeniusPower",    
    "gfpStrategy",
    "Naive", 
    "Safe", 
    
--F-thresholds computations (FThresholds.m2)
    "BinomialCheck",
    "DiagonalCheck", 
    "estFPT",
    "FinalCheck",    
    "FPTApproxList",     
    "FTApproxList",
    "FTHatApproxList", 
    "guessFPT",
    "isFJumpingNumberPoly",
    "isFPTPoly",
    "MultiThread",    
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

--F-thresholds of special families of polynomials (SpecialFThresholds.m2)
    "binomialFPT",
    "diagonalFPT",
    "factorList",    
    "findCPBelow",
    "FPT2VarHomog",     
    "FPT2VarHomogInternal",
    "isBinomial",
    "isCP",
    "isDiagonal",
    "isInLowerRegion",
    "isInUpperRegion",
    "MaxExp",
    "Nontrivial",    
    "PrintCP",
    "setFTData",
    "splittingField",

-- parameterTestIdeal.m2
    "canonicalIdeal",
    "findusOfIdeal",
    "testModule",
    "randomSubset",
    "isCohenMacaulay",
    "isFrational",
    "IsLocal", --an option for isCohenMacaulay, isFrational, etc.
    "AssumeCM", --an option for isFrational, if true then the function won't check if the ring is CM.

-- testIdeals.m2
    "findQGorGen",
    "findTestElementAmbient",
    "tauAOverPEMinus1Poly",
    "tauGor", --needs optimization
    "tauGorAmb",--needs optimization
    "tauNonPrincipalAOverPEPoly",    
    "tauPoly",
--    "tauQGor",    --Karl removed  since it is subsumed by the new testIdeal
--    "tauQGorAmb", --Karl removed  since it is subsumed by the new testIdeal
    "MaxCartierIndex", --the cartier index limit in the test ideal method
    "testIdeal",

-- Other.m2
    "frobenius", 
    "fSig",
    "HSL", 
    "imageOfRelativeCanonical",
    "imageOfTrace", --doesn't work! 
    "isFPure",  
    "isFRegularPoly",  
    "isFRegularQGor",  
    "isMapSplit",
    "isSharplyFPurePoly",
    "sigmaAOverPEMinus1Poly", 
    "sigmaQGorAmb", --needs optimization  
    "sigmaAOverPEMinus1QGor",  --needs optimization 
 
-- Other
       "findAllCompatibleIdeals", ---MK	
    
	"FFiniteSupport", ---MK
	"findGeneratingMorphisms", ---MK
	"FPureIdeals",
 	"FullMap", ---Karl
	"generatingMorphism", ---MK
	"generatingRoot", ---MK
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

load "./Fsing/frobeniusPowers.m2"

load "./Fsing/compatiblySplit.m2"


load "./Fsing/FPure.m2"

load "./Fsing/FFiniteSupport.m2"



load "./Fsing/parameterTestIdeal.m2"


load "./Fsing/FThresholds.m2"

load "./Fsing/SpecialFThresholds.m2"


load "./Fsing/testIdeals.m2"

load "./Fsing/Other.m2"

load "./Fsing/FsingDocs.m2"

beginDocumentation()

load "./Fsing/EthRootsDoc.m2"

load "./Fsing/FThresholdsDoc.m2"

load "./Fsing/SpecialFThresholdsDoc.m2"

load "./Fsing/compatiblySplitDoc.m2"

load "./Fsing/FFiniteSupportDoc.m2"

load "./Fsing/generatingMorphismDoc.m2"

load "./Fsing/parameterTestIdealDoc.m2"

load "./Fsing/FPureDoc.m2"

load "./Fsing/EthRootsTest.m2"
