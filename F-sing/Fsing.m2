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
     {Name => "Drew Ellingson",
     Email => "drewtell@umich.edu"
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
    "basePExp",    
    "carryTest",  
    "denom",   
    "digit", 	   
    "divideFraction",
    "fasterFloorLog",
    "firstCarry", 
    "floorLog",
    "fracPart", 
    "getCanVector",
    "getNumAndDenom", 
    "maxIdeal", 
    "multOrder",
    "NoZeroC", --option to force certain behavior from a function
    "num",
    "taxicabNorm",
    "truncatedBasePExp",
    
--ethRootFunctions (EthRoots.m2)
    "ascendIdeal", --Karl (still needs more tests / documentation)
    "AscentCount",
    "boundLargestCompatible", ---MK
    "ethRoot",
    "ethRootRingElements",   
    "EthRootStrategy",  
    "frobeniusRoot",  
    "FrobeniusRootStrategy",  
    "getFieldGenRoot",
    "minimalCompatible",
    "MonomialBasis",	
    "Substitution",
    
--Frobenius Powers (frobeniusPowers.m2)
    "fastExp",
    "frobenius",
    "frobeniusPower",
    "gfpStrategy",
    "Naive", 
    "Safe", 
    
--F-thresholds computations (FThresholds.m2)
    "BinaryFormCheck",
    "binarySearch1",
    "binarySearchRecursive",
    "BinomialCheck",
    "ComputePreviousNus",
    "DiagonalCheck", 
    "estFPT", --Karl (and others, Pedro?, maybe should just be called fpt?)
    "FinalCheck", 
    "fpt",   
    "fpt1",   
    "FPTApproxList",     
    "FTApproxList",
    "FTHatApproxList", 
    "guessFPT", --Karl (probably should be incorporated into estFPT
    "isFJumpingNumberPoly", --Karl (should be redone, so as not to assume a polynomial ring)
    "isFPTPoly", --Karl (should be redone, so as not to assume a polynomial ring)
    "linearSearch",
    "MultiThread",
    "newNu",
    "newNuHat", 
    "newNuHatList",
    "newNuList",   
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
    "SearchFunction",
    "TestFunction",
    "testGenFrobeniusPower",
    "testPower",
    "testRoot",
    "UseColonIdeals",

--F-thresholds of special families of polynomials (SpecialFThresholds.m2)
    -- Eventually, only binomialFPT, diagonalFPT, and binaryFormFPT should  
    -- be exported from this section **PT
    "binaryFormFPT",     
    "binaryFormFPTInternal",
    "binomialFPT",
    "diagonalFPT",
    "factorList",    
    "findCPBelow",
    "isCP",
    "isInLowerRegion",
    "isInUpperRegion",
    "MaxExp",
    "Nontrivial",    
    "PrintCP",
    "setFTData",
    "splittingField",

-- parameterTestIdeal.m2
    "AssumeCM", --an option for function, if true, then the function will do less work.
    "AssumeReduced", --an option telling functions to assume a ring is reduced.
    "AssumeNormal", --an option telling functions to assume a ring is normal.
    "canonicalIdeal", --Karl (still needs more tests / documentation), this is based on Moty's old code.
    "findusOfIdeal", --Karl (this is Moty's find u function, but it returns a list if Macaulay2 doesn't identify 1 element).
    "isCohenMacaulay", --Karl (added recently, if anyone has ideas to improve this...)
    "isFrational", --Karl (added recently).
    "IsLocal", --an option for isCohenMacaulay, isFrational, etc.
    "randomSubset",
    "testModule", --Karl (this subsumes a bunch of older functions)
    
-- Finjective.m2
    "HSLGModule", --produces the non-F-injective module, ie the submodule of the canonical module
    "isFinjective",
    "CanonicalStrategy", --how to check F-injectivity on the canonical module (Ext or Katzman)
    "Katzman", --an option for CanonicalStrategy

-- testIdeals.m2
    "findQGorGen", --Karl (this finds y such that I^{[p^e]} : I = (y) + I^{[p^e]}, if it exists)
    "findTestElementAmbient", --Karl (my students Marcus and Dan did some improvements on this recently, it doesn't compute the whole Jacobian, it just looks at random minors until it finds a good one, it can be much much faster)
    "MaxCartierIndex", --the cartier index limit in the test ideal method
    "tauAOverPEMinus1Poly", --Karl (this should be removed)
    "tauGor", --Karl (this should be removed)
    "tauGorAmb",--Karl (this should be removed)
    "tauNonPrincipalAOverPEPoly", --Karl (this should be removed)
    "tauPoly", --Karl (this should be removed)
--    "tauQGor",    --Karl removed  since it is subsumed by the new testIdeal
--    "tauQGorAmb", --Karl removed  since it is subsumed by the new testIdeal
    "testIdeal", --Karl (the new version)

-- Other.m2
    "fSig",
    "HSL", 
    "imageOfRelativeCanonical",
    "imageOfTrace", --doesn't work! 
    "isFPure",  
    "isFRegularPoly",  --Karl : this should be removed / replaced with isFRegular
    "isFRegularQGor",  --Karl : this should be removed / replaced with isFRegular
    "isMapSplit",
    "isSharplyFPurePoly", --Karl needs to be redone
    "sigmaAOverPEMinus1Poly",  --Karl needs to be redone
    "sigmaAOverPEMinus1QGor",  --Karl needs to be redone 
    "sigmaQGorAmb", --Karl needs to be redone
 
-- Other
    "FFiniteSupport", ---MK
    "findAllCompatibleIdeals", ---MK	   
    "findGeneratingMorphisms", ---MK
    "FPureIdeals",
    "FullMap", ---Karl
    "generatingMorphism", ---MK
    "generatingRoot" ---MK
--    "paraTestModule", ---MK
--    "paraTestModuleAmbient" ---MK  
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

load "./Fsing/Finjective.m2"

load "./Fsing/FThresholds.m2"

load "./Fsing/SpecialFThresholds.m2"

load "./Fsing/testIdeals.m2"

load "./Fsing/Other.m2"

beginDocumentation()

load "./Fsing/BasicFunctionsDoc.m2"

load "./Fsing/frobeniusPowersDocs.m2"

load "./Fsing/FsingDocs.m2"

load "./Fsing/EthRootsDoc.m2"

load "./Fsing/FThresholdsDoc.m2"

load "./Fsing/SpecialFThresholdsDoc.m2"

load "./Fsing/compatiblySplitDoc.m2"

load "./Fsing/FFiniteSupportDoc.m2"

load "./Fsing/generatingMorphismDoc.m2"

load "./Fsing/parameterTestIdealDoc.m2"

load "./Fsing/FPureDoc.m2"

-- TESTS

load "./Fsing/BasicFunctionsTests.m2"

load "./Fsing/EthRootsTest.m2"

load "./Fsing/SpecialFThresholdsTest.m2"

load "./Fsing/frobeniusPowersTest.m2"
