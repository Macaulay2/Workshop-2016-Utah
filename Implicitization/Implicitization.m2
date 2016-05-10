newPackage("Implicitization",
    Headline => "Implicitization",
    Version => "0.1",
    Date => "August 10, 2015",
    Authors => {
        {Name => "Justin Chen"},
        {Name => "Joe Kileel"}
        },
    PackageImports => {"Bertini"},
    PackageExports => {"NumericalAlgebraicGeometry"},
    DebuggingMode => true,
    HomePage => "https://math.berkeley.edu/~jchen/"
    )
    export {
        "rationalForm",
        "numericalEvaluate",
            "isPolyMap",
	"numericalImageSample",
	"numericalImageDim",
            "Tol",
	"numericalImageDegree",
	"numericalImageMembership",
	-- "numericalImageDecompose",
	"numericalImageHilbert",
	-- "numericalImageEquations"
        "getRandomImagePoint",
	"monodromyLoop"
    }

debug NumericalAlgebraicGeometry

rationalForm = method() -- Writes a map as a list of 2-element sequences, given by (numerator, denominator)
rationalForm List := List => F -> (
    F/(f -> if instance(f, Sequence) then f else (f, 1_(ring f)))
)

numericalEvaluate = method(Options => {symbol isPolyMap => false})
numericalEvaluate (List, Thing) := Thing => opts -> (F, p) -> (
    if opts.isPolyMap == true then return point evaluate(polySystem F, p);
    F = rationalForm F;
    numerators := flatten entries evaluate(polySystem(F/(f -> f#0)), p);
    denominators := flatten entries evaluate(polySystem(F/(f -> f#1)), p);
    point{toList(0..<#F)/(i -> (numerators#i)/(denominators#i))}
)

numericalImageSample = method()
numericalImageSample (List, NumericalVariety, ZZ) := List => (F, NV, n) -> (
    (components NV)/(W -> bertiniSample(n, W, Verbose => false))/(c -> c/(p -> (p, numericalEvaluate(F, p))))
)
numericalImageSample (List, Ideal, ZZ) := List => (F, I, n) -> (
    {*if dim I == 0 then (
        print "Running BertiniZeroDimSolve ...";
        NV := bertiniZeroDimSolve(I_*, Verbose => false);
        return NV/(isolatedPt -> {(isolatedPt, numericalEvaluate(F, isolatedPt))});
    ) else *}if not I == 0 then (
        print "Running BertiniPosDimSolve ...";
        return numericalImageSample(F, bertiniPosDimSolve(I_*, Verbose => false), n)
    ) else (
        p := point{getRandomPoint(gens ring I, coefficientRing ring I)};
        return {{(p, numericalEvaluate(F, p))}}
    );
)
numericalImageSample (List, NumericalVariety) := List => (F, NV) -> numericalImageSample(F, NV, 1)
numericalImageSample (List, Ideal) := List => (F, I) -> numericalImageSample(F, I, 1)

numericalImageSample2 = method()
numericalImageSample2 (List, Ideal, ZZ) := List => (F, I, N) -> (
    NV := bertiniPosDimSolve(I_*);
    numericalImageSample(F, NV, ceiling(5*N/#components(NV)))     
)

numericalImageDim = method(Options => {symbol Tol => 1e2, symbol isPolyMap => false})
numericalImageDim (List, Ideal, List) := ZZ => opts -> (F, I, samplePoints) -> ( --I = radical I;
    GraphIdeal := makeGraphIdeal(F, I);
    jacobianGraph := jacobian polySystem GraphIdeal_*;
    print "Evaluating Jacobians ...";
    time evaluatedJacobians := samplePoints/(c -> evaluate(jacobianGraph, point{join(((c#0)#0)#Coordinates, ((c#0)#1)#Coordinates)}));
    print "Computing ranks ...";
    time max(evaluatedJacobians/(j -> (#gens ring GraphIdeal - (numericalRank(j, Threshold=>opts.Tol)) - (#gens ring I - numericalRank(j_{0..<#gens ring I},Threshold=>opts.Tol)))))
)
numericalImageDim (List, Ideal) := ZZ => opts -> (F, I) -> (
    numericalImageDim(F, I, numericalImageSample(F, I), opts)
)	     
numericalImageDim (List, NumericalVariety) := ZZ => opts -> (F, NV) -> (
    numericalImageDim(F, ideal(components(NV)#0), numericalImageSample(F, NV), opts)
)

numericalImageInvariants = method()
numericalImageInvariants (List, Ideal) := List => (F,I) -> (
    makeBertiniInput(F, I)
)

makeGraphIdeal = method()
makeGraphIdeal (List, Ideal) := Ideal => (F, I) -> (
    --F = F/(f -> if instance(f, Sequence) then f/(h->promote(h,ring I)) else (promote(f, ring I), 1_(ring I)));
    F = rationalForm F;
    targetvariable := symbol targetvariable;
    TARGETVARIABLES := (targetvariable_1..targetvariable_#F);
    GRAPHRING := ring I ** (coefficientRing ring I)[TARGETVARIABLES];
    SOURCETOGRAPH := map(GRAPHRING, ring I);
    TARGETTOGRAPH := map(GRAPHRING, (coefficientRing ring I)[TARGETVARIABLES]);
    F = toList(1..#F)/(i -> (SOURCETOGRAPH((F#(i-1))#1))*(TARGETTOGRAPH(targetvariable_i)) - SOURCETOGRAPH((F#(i-1))#0));
    ideal(F) + SOURCETOGRAPH(I)
)

makeBertiniInput = method()
makeBertiniInput (List, Ideal) := (F, I) -> (
    GraphIdeal := makeGraphIdeal(F, I);
    inputFile := openOut temporaryFileName();
    inputFile << "CONFIG" << endl;
    inputFile << " TrackType:1;" << endl;
    inputFile << "END;" << endl;
    inputFile << "INPUT" << endl;
    inputFile << " function ";
    for i from 1 to #(GraphIdeal_*)-1 do inputFile << "equation"|toString(i)|",";
    inputFile << "equation"|toString(#(GraphIdeal_*))|";" << endl;
    inputFile << " variable_group ";
    for i from 0 to #(gens ring GraphIdeal)-2 do inputFile << toString((gens ring GraphIdeal)_(0,i))|",";
    inputFile << toString((vars ring GraphIdeal)_(0,#(gens ring GraphIdeal)-1))|";" << endl;
    inputFile << endl;
    for i from 1 to #(GraphIdeal_*) do inputFile << "equation"|toString(i)|" = "|toString((GraphIdeal_*)#(i-1))|";" << endl;
    inputFile << "END;" << endl << close;
    toString(inputFile)
)

numericalImageMembership = method()
numericalImageMembership (List, Ideal, Thing, RR) := (F,I,p,tol) -> (

)
numericalImageMembership (List, Ideal, Thing) := (F, I, p) -> numericalImageDim(F, I, p, 0.0000000001)

numericalImageHilbert = method()
numericalImageHilbert (List, Ideal, ZZ) := ZZ => (F,I,d) -> (
    local targetvariable;
    targetvariable = symbol targetvariable;
    C := coefficientRing ring I;
    TARGETRING := C[(targetvariable_1..targetvariable_#F)];
    monos := polySystem transpose basis(d, TARGETRING);
    print "Applying F to sample points ...";
    time samplePoints := toList(1..monos.NumberOfPolys)/(i -> getRandomImagePoint(F, I));
    samplePoints = samplePoints/(p -> point{(p#Coordinates)/(norm(2,p))});
    print "Creating interpolation matrix ...";
    time interpolation := matrix(samplePoints/(p -> flatten entries evaluate(monos, p)));
    print "Computing SVD ...";
    time (S, U, Vt) := SVD(interpolation, DivideConquer => true);
    S = toList S;
    print S;
    for i from 1 to #S-1 do if S#i*10^5 < S#(i-1) then return #S - i;
    return 0;
)

getRandomPoint = method()
getRandomPoint (List, Ring) := List => (L, C) -> (
    L/(l -> random C)
)

getRandomImagePoint = method()
getRandomImagePoint (List, Ideal) := Thing => (F, I) -> (
    numericalEvaluate(F, point{getRandomPoint(gens ring I, coefficientRing ring I)})
)

numericalImageDegree = method()
numericalImageDegree (List, Ideal, String) := ZZ => (F, I, software) -> (
    print "Sampling point from domain ...";
    startPoint := ((numericalImageSample(F, I))#0)#0;
    print "Finding dimension of image ...";
    imageDim := numericalImageDim(F, I, {{startPoint}});
    extra := (toList(0..<dim I - imageDim)/(n -> random(1, ring I)))/(l -> l - (evaluate(polySystem{l}, startPoint#0))_(0,0));
    traceTestAttempts := 0;
    monodromyOutput := monodromyLoop(F, I, {imageDim, startPoint, traceTestAttempts, createSystem(imageDim, F, {extra, I, startPoint#1}), extra, startPoint#1}, software);
    if not traceTest(monodromyOutput, F, software) then (
	print "Failed trace test!";
	traceTestAttempts = 1;
	while traceTestAttempts < 5 do (
            monodromyOutput = monodromyLoop(F, I, {imageDim, startPoint, traceTestAttempts, createSystem(imageDim, F, {extra, I, startPoint#1}), extra}, software);
	    if not traceTest(monodromyOutput, F, software) then ( 
		traceTestAttempts = traceTestAttempts + 1;
		print "Failed trace test!";
            ) else break;
	);
    );
    #(monodromyOutput#1)
)

createSystem = method()
createSystem (ZZ, List, List) := List => (d, F, L) -> ( -- creates random linear slice
    extra := L#0;
    I := L#1;
    startPoint := (L#2)#Coordinates;
    C := coefficientRing ring I;
    slice := random(C^d, C^(#F));
    pullbackSystem := flatten entries(slice*(transpose matrix{F}));
    pullbackSystem = pullbackSystem - flatten entries(slice*transpose matrix{startPoint});
    pullbackSystem = join(pullbackSystem, extra);
    if not I == 0 then pullbackSystem = join(pullbackSystem, I_*);
    pullbackSystem
)
createSystem (List, List, Boolean) := List => (S, F, isTraceTest) -> ( -- creates (parallel) random translation(s)
    C := coefficientRing ring first F;
    v := toList(0..<min(#S, #F))/(i -> 10*random C);
    O := toList(0..<#S)/(i -> if i < min(#S, #F) then S#i - v#i else S#i);
    gamma := {random(C), random(C)};
    if isTraceTest == true then {toList(0..<#S)/(i -> if i < min(#S, #F) then S#i - gamma#0*v#i else S#i), toList(0..<#S)/(i -> if i < min(#S, #F) then S#i - gamma#1*v#i else S#i),gamma} else O
)

monodromyLoop = method()
monodromyLoop (List, Ideal, List, String) := List => (F, I, P, software) -> (
    imageDim := P#0;
    startPoint := P#1;
    traceTestAttempts := P#2;
    startSystem := P#3;
    extra := P#4;
    C := coefficientRing ring I;
    startSols := {startPoint#0};
    intersectionPoints := {startPoint#1};
    if software == "parametricM2" then (
    	a := symbol a;
    	R := C(monoid [gens ring I, a_(1,1)..a_(#F,imageDim)]);
    	toR := map(R, ring I, take(gens R,numgens ring I));
    	slice := genericMatrix(R,R_(numgens ring I),imageDim,#F);
    	pullbackSystem := flatten entries(slice * toR transpose matrix{F});
    	pt1 := random(C^imageDim,C^#F);
    	pullbackSystem = pullbackSystem - flatten entries(pt1*transpose matrix startPoint#1);
    	pullbackSystem = join(pullbackSystem, extra/(f->toR f));
    	if I != 0 then pullbackSystem = join(pullbackSystem, (toR I)_*);
    	PS := polySystem pullbackSystem;
    	makeGateMatrix(PS,Parameters=>drop(gens R,numgens ring I));
    	PH := parametricSegmentHomotopy PS; 
	pt1 = transpose flatten pt1; 
	);
    noNewPoints := 0;
    local downstairs, local currentPointsNum, local startSystem, local mid1, local mid2, local gamma;
    local homotopyStep1, local homotopyStep2, local homotopyStep3;
    G := F;
    J := I;
    if software == "Bertini" then (
        local T;
        T = symbol T;
        HOMOTOPYRING := C[T] ** ring I;
        SOURCETOHOMOTOPYRING := map(HOMOTOPYRING, ring I);
        G = F/(f -> SOURCETOHOMOTOPYRING f);
        J = SOURCETOHOMOTOPYRING I;
        T = first gens HOMOTOPYRING;
        extra = extra/(l -> SOURCETOHOMOTOPYRING l);
	startSystem = startSystem/(f -> SOURCETOHOMOTOPYRING f);
    );
    print "Executing monodromy main loop ...";
    time while noNewPoints < 10 do (
        currentPointsNum = #intersectionPoints;
        mid1 = createSystem(if traceTestAttempts > 0 then createSystem(imageDim, G, {extra, J, startPoint#1}) else startSystem, G, false);
        mid2 = createSystem(if traceTestAttempts > 0 then createSystem(imageDim, G, {extra, J, startPoint#1}) else startSystem, G, false);
        gamma = toList(0..<3)/(i -> random C);
        if software == "Bertini" then (
            homotopyStep1 = bertiniTrackHomotopy(T, T*gamma#0*startSystem + (1-T)*mid1, startSols, Verbose => false);
            homotopyStep2 = bertiniTrackHomotopy(T, T*gamma#1*mid1 + (1-T)*mid2, homotopyStep1, Verbose => false);
            homotopyStep3 = bertiniTrackHomotopy(T, T*gamma#2*mid2 + (1-T)*startSystem, homotopyStep2, Verbose => false);
        ) else if software == "parametricM2" then (
	    pt2 := random(C^(numcols PH.Parameters//2),C^1);
	    pt3 := random(C^(numcols PH.Parameters//2),C^1);
	    homotopyStep1 = trackHomotopy(specialize (PH, pt1||pt2), startSols);
            homotopyStep2 = trackHomotopy(specialize (PH, pt2||pt3), homotopyStep1);
            homotopyStep3 = trackHomotopy(specialize (PH, pt3||pt1), homotopyStep2);
	    homotopyStep3 = select(homotopyStep3, p -> p#SolutionStatus == Regular);
	) else (
            homotopyStep1 = track(gamma#0*startSystem, mid1, startSols);
            homotopyStep2 = track(gamma#1*mid1, mid2, homotopyStep1);
            homotopyStep3 = track(gamma#2*mid2, startSystem, homotopyStep2);
	    homotopyStep3 = select(homotopyStep3, p -> p#SolutionStatus == Regular);
        );
	for p in homotopyStep3 do (
            downstairs = numericalEvaluate(F, p);
            if all(intersectionPoints, q -> not downstairs == q) then (
                intersectionPoints = append(intersectionPoints, downstairs);
                startSols = append(startSols, p);
            );
        );
        print ("Points found: " | #intersectionPoints);
        if currentPointsNum == #intersectionPoints then noNewPoints = noNewPoints + 1 else noNewPoints = 0;
    );
    {startSystem, startSols}
)

traceTest = method()
traceTest (List, List, String) := Boolean => (S, F, software) -> (
    C := coefficientRing ring first F;
    startSystem := S#0;
    startSols := S#1;
    traceSlices := createSystem(startSystem, F, true);
    gamma := traceSlices#2;
    local traceStep1, local traceStep2;
    if software == "Bertini" then (
        T := first gens ring first startSystem;
        traceStep1 = bertiniTrackHomotopy(T, T*startSystem + (1-T)*traceSlices#0, startSols, Verbose => false);
        traceStep2 = bertiniTrackHomotopy(T, T*startSystem + (1-T)*traceSlices#1, startSols, Verbose => false);
    ) else (
        traceStep1 = track(startSystem, traceSlices#0, startSols);
        traceStep2 = track(startSystem, traceSlices#1, startSols);
    );
    P0 := sumPointList(startSols/(p -> numericalEvaluate(F, p)));
    P1 := sumPointList(traceStep1/(p->numericalEvaluate(F,p)));
    P2 := sumPointList(traceStep2/(p->numericalEvaluate(F,p)));
    discrepancy := (gamma#1)*(P1#Coordinates - P0#Coordinates) - (gamma#0)*(P2#Coordinates - P0#Coordinates);
    point{{sum(toList(0..<#F)/(i->random(C)*discrepancy#i))}} == point{{0_C}}
)

Point + Point := (P, Q) -> point{P#Coordinates + Q#Coordinates}

sumPointList = method()
sumPointList List := Point => L -> point{sum(L/(p -> p#Coordinates))}

beginDocumentation()

--Documention--
--<<docTemplate
doc ///
    Key
    	Implicitization
    Headline
    	a Macaulay2 implicitization package
    Description
    	Text
	    Allows for user-friendly computation of basic invariants of the image of a rational map, based on interpolation and homotopy continuation.  Geared toward large-scale and applied problems.
///

TEST ///
d = dim ker map(QQ[x,y,z,w]/ideal(x^3 - y*z*w), QQ[a..e], {x*w + 2*x*y, x*w-3*y^2, z^2, x^2 + y^2 + z^2 - w^2, 3*y*w - 2*x^2})
R = CC[x,y,z,w]
I = ideal(x^3 - y*z*w)
F = {x*w + 2*x*y, x*w-3*y^2, z^2, x^2 + y^2 + z^2 - w^2, 3*y*w - 2*x^2}
assert(numericalImageDim(F, I) == d)
///

TEST ///
R = CC[s,t]
F = flatten entries basis(4, R) - set{s^2*t^2} -- Rational quartic curve in P^3
h5 = numcols basis(5, ker map(QQ[s,t], QQ[x,y,z,w], {s^4,s^3*t,s*t^3,t^4}))
assert(numericalImageHilbert(F, ideal(0_R), 5) == h5)
///

end--

restart
load"Implicitization.m2"
loadPackage("Implicitization", Reload => true)
uninstallPackage "Implicitization"
installPackage("Implicitization", RemakeAllDocumentation => true)
viewHelp "Implicitization"
check "Implicitization"




--- Alexander-Hirschowitz: cubic forms in 5 variables (expect generic rank 7)

R = CC[a_{0,0}..a_{6,4}, b_0..b_6]
S = CC[x_0..x_4]
cubics = (entries transpose genericMatrix(R, 5, 7))/(r -> flatten entries((map(R, S, r))(basis(3, S))))
F = sum(toList(0..<7)/(i -> (cubics#i)/(a -> b_i*a)))
numericalImageDim(F, ideal 0_R)




-- Border rank 2 3x3x3 tensors

R = CC[a_0..a_2, b_0..b_2, c_0..c_2, d_0..d_2, e_0..e_2, f_0..f_2]
F = toList({0,0,0}..{2,2,2})/(t -> a_(t#0)*b_(t#1)*c_(t#2) + d_(t#0)*e_(t#1)*f_(t#2))
I = ideal 0_R
numericalImageDim(F,I)
numericalImageHilbert(F,I,2)
numericalImageHilbert(F,I,3)




-- Rational normal curve

restart
R = CC[s,t]
F = flatten entries basis(66, R)
I = ideal(0_R)
numericalImageDegree(F, I, "M2") -- should be degree 6
F = delete(s^2*t^2, flatten entries basis(4, R))
numericalImageDegree(F, I, "M2") -- should be degree 4
 S= QQ[s,t]
 degree(ker map(S, QQ[x_1..x_366], basis(365,S)))



R = CC[x,y,z,w]
I = ideal(x^10 - y*z*w)
F = {x*w + 2*x*y, x*w-3*y^2, z^2, x^2 + y^2 + z^2 - w^2, 3*y*w - 2*x^2}
numericalImageDim(F, ideal(R_1)^2)
numericalImageDim(F, ideal(R_2)^2)
numericalImageDim(F, ideal(0_R))
numericalImageDim(F, (ideal gens R)^2)
time numericalImageHilbert(F,I,4)
numericalImageDegree(F,I, "M2")



R = CC[x_1..x_15]
F = (minors(3, genericMatrix(R, 3, 5)))_*
I = ideal(0_R)
numericalImageDegree(F, I, "M2") -- should be degree 5



R = CC[a..f]
p3 = flatten entries(matrix{{a,b,c,d}}*random(CC^4, CC^7))
p1 = flatten entries(matrix{{e,f}}*random(CC^2, CC^7))
numericalImageDegree(toList(0..<7)/(i -> (p1#i)*product(delete(p3#i, p3))), ideal 0_R, "M2") -- should be degree 10


Makegraphideal(F,I)
Makebertiniinput(F,I)
numericalImageSample(F,I)
time numericalImageDim(F,I)

R = CC_100[a1,a2,a3,a4,b1,b2,b3,b4,s1,s2,s3,t1,t2,t3];

A = a1^2 + a2^2 + a3^2 + a4^2 - 1;
B = b1^2 + b2^2 + b3^2 + b4^2 - 1;

T111 = b1^2*s1+b2^2*s1-b3^2*s1-b4^2*s1-a1^2*t1-a2^2*t1+a3^2*t1+a4^2*t1;
T112 = 2*b2*b3*s1-2*b1*b4*s1-a1^2*t2-a2^2*t2+a3^2*t2+a4^2*t2;
T113 = 2*b1*b3*s1+2*b2*b4*s1-a1^2*t3-a2^2*t3+a3^2*t3+a4^2*t3;
T121 = b1^2*s2+b2^2*s2-b3^2*s2-b4^2*s2-2*a2*a3*t1+2*a1*a4*t1;
T122 = 2*b2*b3*s2-2*b1*b4*s2-2*a2*a3*t2+2*a1*a4*t2;
T123 = 2*b1*b3*s2+2*b2*b4*s2-2*a2*a3*t3+2*a1*a4*t3;
T131 = b1^2*s3+b2^2*s3-b3^2*s3-b4^2*s3-2*a1*a3*t1-2*a2*a4*t1;
T132 = 2*b2*b3*s3-2*b1*b4*s3-2*a1*a3*t2-2*a2*a4*t2;
T133 = 2*b1*b3*s3+2*b2*b4*s3-2*a1*a3*t3-2*a2*a4*t3;
T211 = 2*b2*b3*s1+2*b1*b4*s1-2*a2*a3*t1-2*a1*a4*t1;
T212 = b1^2*s1-b2^2*s1+b3^2*s1-b4^2*s1-2*a2*a3*t2-2*a1*a4*t2;
T213 = -2*b1*b2*s1+2*b3*b4*s1-2*a2*a3*t3-2*a1*a4*t3;
T221 = 2*b2*b3*s2+2*b1*b4*s2-a1^2*t1+a2^2*t1-a3^2*t1+a4^2*t1;
T222 = b1^2*s2-b2^2*s2+b3^2*s2-b4^2*s2-a1^2*t2+a2^2*t2-a3^2*t2+a4^2*t2;
T223 = -2*b1*b2*s2+2*b3*b4*s2-a1^2*t3+a2^2*t3-a3^2*t3+a4^2*t3;
T231 = 2*b2*b3*s3+2*b1*b4*s3+2*a1*a2*t1-2*a3*a4*t1;
T232 = b1^2*s3-b2^2*s3+b3^2*s3-b4^2*s3+2*a1*a2*t2-2*a3*a4*t2;
T233 = -2*b1*b2*s3+2*b3*b4*s3+2*a1*a2*t3-2*a3*a4*t3;
T311 = -2*b1*b3*s1+2*b2*b4*s1+2*a1*a3*t1-2*a2*a4*t1;
T312 = 2*b1*b2*s1+2*b3*b4*s1+2*a1*a3*t2-2*a2*a4*t2;
T313 = b1^2*s1-b2^2*s1-b3^2*s1+b4^2*s1+2*a1*a3*t3-2*a2*a4*t3;
T321 = -2*b1*b3*s2+2*b2*b4*s2-2*a1*a2*t1-2*a3*a4*t1;
T322 = 2*b1*b2*s2+2*b3*b4*s2-2*a1*a2*t2-2*a3*a4*t2;
T323 = b1^2*s2-b2^2*s2-b3^2*s2+b4^2*s2-2*a1*a2*t3-2*a3*a4*t3;
T331 = -2*b1*b3*s3+2*b2*b4*s3-a1^2*t1+a2^2*t1+a3^2*t1-a4^2*t1;
T332 = 2*b1*b2*s3+2*b3*b4*s3-a1^2*t2+a2^2*t2+a3^2*t2-a4^2*t2;
T333 = b1^2*s3-b2^2*s3-b3^2*s3+b4^2*s3-a1^2*t3+a2^2*t3+a3^2*t3-a4^2*t3;

F = {b1^2*s1+b2^2*s1-b3^2*s1-b4^2*s1-a1^2*t1-a2^2*t1+a3^2*t1+a4^2*t1,
2*b2*b3*s1-2*b1*b4*s1-a1^2*t2-a2^2*t2+a3^2*t2+a4^2*t2,
2*b1*b3*s1+2*b2*b4*s1-a1^2*t3-a2^2*t3+a3^2*t3+a4^2*t3,
b1^2*s2+b2^2*s2-b3^2*s2-b4^2*s2-2*a2*a3*t1+2*a1*a4*t1,
2*b2*b3*s2-2*b1*b4*s2-2*a2*a3*t2+2*a1*a4*t2,
2*b1*b3*s2+2*b2*b4*s2-2*a2*a3*t3+2*a1*a4*t3,
b1^2*s3+b2^2*s3-b3^2*s3-b4^2*s3-2*a1*a3*t1-2*a2*a4*t1,
2*b2*b3*s3-2*b1*b4*s3-2*a1*a3*t2-2*a2*a4*t2,
2*b1*b3*s3+2*b2*b4*s3-2*a1*a3*t3-2*a2*a4*t3,
2*b2*b3*s1+2*b1*b4*s1-2*a2*a3*t1-2*a1*a4*t1,
b1^2*s1-b2^2*s1+b3^2*s1-b4^2*s1-2*a2*a3*t2-2*a1*a4*t2,
-2*b1*b2*s1+2*b3*b4*s1-2*a2*a3*t3-2*a1*a4*t3,
2*b2*b3*s2+2*b1*b4*s2-a1^2*t1+a2^2*t1-a3^2*t1+a4^2*t1,
b1^2*s2-b2^2*s2+b3^2*s2-b4^2*s2-a1^2*t2+a2^2*t2-a3^2*t2+a4^2*t2,
-2*b1*b2*s2+2*b3*b4*s2-a1^2*t3+a2^2*t3-a3^2*t3+a4^2*t3,
2*b2*b3*s3+2*b1*b4*s3+2*a1*a2*t1-2*a3*a4*t1,
b1^2*s3-b2^2*s3+b3^2*s3-b4^2*s3+2*a1*a2*t2-2*a3*a4*t2,
-2*b1*b2*s3+2*b3*b4*s3+2*a1*a2*t3-2*a3*a4*t3,
-2*b1*b3*s1+2*b2*b4*s1+2*a1*a3*t1-2*a2*a4*t1,
2*b1*b2*s1+2*b3*b4*s1+2*a1*a3*t2-2*a2*a4*t2,
b1^2*s1-b2^2*s1-b3^2*s1+b4^2*s1+2*a1*a3*t3-2*a2*a4*t3,
-2*b1*b3*s2+2*b2*b4*s2-2*a1*a2*t1-2*a3*a4*t1,
2*b1*b2*s2+2*b3*b4*s2-2*a1*a2*t2-2*a3*a4*t2,
b1^2*s2-b2^2*s2-b3^2*s2+b4^2*s2-2*a1*a2*t3-2*a3*a4*t3,
-2*b1*b3*s3+2*b2*b4*s3-a1^2*t1+a2^2*t1+a3^2*t1-a4^2*t1,
2*b1*b2*s3+2*b3*b4*s3-a1^2*t2+a2^2*t2+a3^2*t2-a4^2*t2,
b1^2*s3-b2^2*s3-b3^2*s3+b4^2*s3-a1^2*t3+a2^2*t3+a3^2*t3-a4^2*t3}

I = ideal(A,B)

G = {(T111,T333),
(T112,T333),
(T113,T333),
(T121,T333),
(T122,T333),
(T123,T333),
(T131,T333),
(T132,T333),
(T133,T333),
(T211,T333),
(T212,T333),
(T213,T333),
(T221,T333),
(T222,T333),
(T223,T333),
(T231,T333),
(T232,T333),
(T233,T333),
(T311,T333),
(T312,T333),
(T313,T333),
(T321,T333),
(T322,T333),
(T323,T333),
(T331,T333),
(T332,T333),
(T333,T333)}

getRandomPointTrifocal = method()
getRandomPointTrifocal (List, ZZ) := List => (variables, tol) -> (
    L := toList(0..13);
    L0 := L/(i -> random CC_tol);
    alist := L0_{0,1,2,3}/norm2(L0_{0,1,2,3});
    blist := L0_{4,5,6,7}/norm2(L0_{4,5,6,7});
    alist | blist | L0_{8..13}
)




restart
loadPackage("Implicitization", Reload => true)
upstairs = random(CC^3, CC^5)
stiefelRing = CC[x_1..x_15]
pluckerRing = CC[p_1..p_10]
--pluckerMap = map(stiefelRing, pluckerRing, (minors(3, genericMatrix(stiefelRing,3,5)))_*)
downstairs = apply((minors(3, genericMatrix(stiefelRing,3,5)))_*, m->sub(m,apply(toList(0..<15), i->(gens stiefelRing)#i=>(flatten entries(upstairs))#i)))
mySlice = random(CC^7,CC^10)
mySlice = flatten entries(mySlice*transpose matrix{gens pluckerRing} - mySlice*transpose matrix{downstairs})--make it pass through downstairs
QstiefelRing = QQ[X_1..X_15]
QpluckerRing = QQ[P_1..P_10]
basechange = map(pluckerRing, QpluckerRing, gens pluckerRing)
QpluckerRelations = (ker(map(QstiefelRing, QpluckerRing, (minors(3, genericMatrix(QstiefelRing,3,5)))_*)))_*
pluckerRelations = apply(QpluckerRelations, r->basechange(r))
zeroDimSystem = flatten join(mySlice,pluckerRelations)
square1 = squareUp(polySystem(zeroDimSystem))
square2 = squareUp(polySystem(zeroDimSystem))
square1 === square2
time solutions1 = solveSystem square1
time solutions2 = solveSystem square2
--distanceMatrix = matrix(toList(0..<12)/(i -> toList(0..<12)/(j -> norm(2,(solutions1#i)#Coordinates - (solutions2#j)#Coordinates))))
trueSolutions = flatten flatten(apply(solutions1, s -> solutions2/(t -> if s == t then s else {})))

use stiefelRing
F = (minors(3, genericMatrix(stiefelRing,3,5)))_*
pullbackSlice = mySlice/(l -> sub(l, toList(0..<#F)/(i -> (gens ring first mySlice)#i => F#i)))
I = ideal 0_stiefelRing
extra = (toList(0..<dim I - 7)/(n -> random(1, ring I)))/(l -> l - (evaluate(polySystem{l}, point upstairs))_(0,0))
monodromyLoop(F, I, {7, {point upstairs, point{downstairs}}, true, join(pullbackSlice, extra), extra}, "Bertini");


--specificMonodromy(F, ideal 0_stiefelRing, {point upstairs, point{downstairs}}, pullbackSlice);


specificMonodromy = method()
specificMonodromy (List, Ideal, Thing, List) := List => (F, I, startPoint, pullbackSlice) -> (
    imageDim := 7;
    C := coefficientRing ring I;
    passedTraceTest := true;
    startSols := {startPoint#0};
    intersectionPoints := {startPoint#1};
    startSystem := pullbackSlice;
    extra := (toList(0..<dim I - imageDim)/(n -> random(1, ring I)))/(l -> l - (evaluate(polySystem{l}, startPoint#0))_(0,0));
    startSystem = join(startSystem, extra);
    --print extra;
    G := F;
    J := I;
    noNewPoints := 0;
    local P, local currentPointsNum, local startSystem, local mid1, local mid2, local gamma;
    local homotopyStep1, local homotopyStep2, local homotopyStep3;
    while noNewPoints < 10 do (
        currentPointsNum = #intersectionPoints;
        mid1 = createSystem(if not passedTraceTest then createSystem(imageDim, G, {extra, J, startPoint#1}) else startSystem, G, false);
        mid2 = createSystem(if not passedTraceTest then createSystem(imageDim, G, {extra, J, startPoint#1}) else startSystem, G, false);
        gamma = toList(0..<3)/(i -> random C);
	homotopyStep1 = track(gamma#0*startSystem, mid1, startSols);
        homotopyStep2 = track(gamma#1*mid1, mid2, homotopyStep1);
        homotopyStep3 = track(gamma#2*mid2, startSystem, homotopyStep2);
	for p in homotopyStep3 do (
            P = numericalEvaluate(F, p);
            print(intersectionPoints/(q -> norm2(q#Coordinates - P#Coordinates)));
            if all(intersectionPoints, q -> not P == q) then (
                intersectionPoints = append(intersectionPoints, P);
                startSols = append(startSols, p);
            );
        );
        print ("Points found: " | #intersectionPoints);
        if currentPointsNum == #intersectionPoints then noNewPoints = noNewPoints + 1 else noNewPoints = 0;
    );
    {startSystem, startSols}
)


