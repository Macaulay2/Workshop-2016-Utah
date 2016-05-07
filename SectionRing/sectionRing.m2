needsPackage "Divisor"

dualToIdeal = method();
dualToIdeal(Ideal) := (I) -> (
	R := ring(I);
	M :=module(I);
	moduleToIdeal(Hom(M,R),IsGraded=>true,ReturnMap=>true)
);

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

globallyGenerated = method();
globallyGenerated(WDiv) := (D) -> (				--Finds the smallest positive number such that O_X(a*D) is globally generated
	a:=1;

	while ((1%(baseLocus(a*D)) == 0) != true) do (		--Uses a binary search to improve speed
		a =2*a;
	);

	upperbound := a;
	lowerbound := ceiling(a/2);

	while (lowerbound < upperbound-1) do (
		a = ceiling((lowerbound + upperbound)/2);
		if ((1%(baseLocus(a*D)) == 0) != true) then (
			lowerbound = a;
		)
		else if ((1%(baseLocus(a*D)) == 0) == true) then (
			upperbound = a;
		);

	);
	upperbound
);


globallyGenerated(Ideal) := (I) -> (
	globallyGenerated(divisor(I))
);

globallyGenerated(Module) := (M) -> (
	globallyGenerated(moduleToDivisor(M))
);

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

isMRegular = method();
isMRegular(CoherentSheaf,CoherentSheaf,ZZ) := (F,G,m) ->(	--Outputs whether of not F is m-regular relative to G
	V := variety(F);
	dV := dim(V);
	j=1;
	bool = true;
	while(j<(dV+1)) do (
		if (bool == true) then(
			if(m!=j) then(
				bool = (HH^j((F**(G^**(m-j)))) == 0);
			)
			else if (m==j) then (
				bool = (HH^j(F) == 0);
			);
		);
		j = j+1;
	);
	bool
);

isMRegular(CoherentSheaf,ZZ) := (F,m) ->(		--Outputs whether F is m-regular (rel O_X(1))
	V :variety(F);
	G = OO_V(1);
	mRegularParticular(F,G,m)
);

isMRegularOO = method();
isMRegularOO(CoherentSheaf,ZZ) := (F,m) -> (			--Tests if all higher cohomologies of F vanish
	V := variety(F);
	dV := dim(V);
	j:=1;
	bool = true;
	while(j<(dV+1)) do (
		if (bool == true) then(
			bool = (HH^j(F) == 0);
		);
		j = j+1;
	);
	bool
);

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

mRegular = method();
mRegular(CoherentSheaf,CoherentSheaf) := (F,G) -> (	--Computes the m-regularity of the sheaf F relative to G
	bool0 := isMRegular(F,G,0);
	m:=0;
	lowerbound:=0;
	upperbound:=0;
	a:=0;
	if (bool0 == true) then (			--Tests negative-regularity using a binary search
		m=-1;
		while (isMRegular(F,G,m)) do (
			m=2*m;
		);

		lowerbound = m;
		upperbound = ceiling(m/2);

		while (lowerbound < upperbound-1) do (
			a = ceiling((lowerbound + upperbound)/2);
			if (isMRegular(F,G,a) != true) then (
				lowerbound = a;
			)
			else if (isMRegular(F,G,a) == true) then (
				upperbound = a;
			);
		);
	)

	else if (bool0==false) then (
		m=1;
		while (isMRegular(F,G,m) != true) do (
			m=2*m;
		);

		upperbound = m;
		lowerbound = ceiling(m/2);

		while (lowerbound < upperbound-1) do (
			a = ceiling((lowerbound + upperbound)/2);
			if (isMRegular(F,G,a) != true) then (
				lowerbound = a;
			)
			else if (isMRegular(F,G,a) == true) then (
				upperbound = a;
			);
		);
	);
	upperbound
);

mRegular(CoherentSheaf) := (F) -> ( 				--m-Regularity of a sheaf (relative OO_V(1))
	V := variety(F);					
	mRegular(F,(OO_V(1)))
);

mRegular(Ideal) := (I) -> ( 					--Returns the number m for which O_X(D) is m-regular, where  I is an ideal,
	R := ring(I);						--and D is the corresponding divisor to I.
	F := sheaf(Hom(module(I),R));
	mRegular(F)
);

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

sectionRing = method();
sectionRing(Ideal) := (I) -> (
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
	R := ring(I);
	l := globallyGenerated(I);
	j := 0;
	bound := l;
	G = first entries gens I;
	J_l = Hom(ideal(apply(G, z->z^l)),R);
	while(j<l) do (
		J_j = Hom(ideal(apply(G, z->z^j)),R);
		bound = max(bound,l*(mRegular(sheaf(J_j),sheaf(J_l)))+j);
		j = j+1;
	);
	bound = bound + 1;
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
	KK:= coefficientRing(R);
	Z := dualToIdeal(I);
	Shift := (Z#1)#0;
	J_1 = reflexifyIdeal((Z#0));
	F_1 = basis(Shift,J_1);
	n_1 = numColumns(F_1);
	F_1 = map(R^(numRows(F_1)),R^(n_1),F_1);
	Map_1 = (gens J_1)*(F_1);				--Map_i Represents the map H^0(O_X(iD)) -> J^(i)
	myVars = {};						--Begins to create a list for the necessary variables
	DegreeList :={};					--and a list of their corresponding degrees.
	i:=1;
	while ( i < bound) do(
		J_i = reflexivePower(i,J_1);
		F_i = basis((Shift*i),J_i);
		n_i = numColumns((F_i));			--Rank of H^0(O_X(iD))
		F_i = map(R^(numRows((F_i))),R^(n_i),F_i);
		Map_i = (gens J_i)*(F_i);
		myVars = (myVars) | {toList(Y_{i,1}..Y_{i,n_i})};  --Creates a list of variables, grouped internally by degree for convenience
		l=0;
		while(l<n_i) do(				--Creates the degree list for the section ring
			DegreeList = DegreeList | toList({i});
			l = l+1;
		);
		i=i+1;
	);
	
	Vars := flatten myVars;

	S := KK [Vars,Degrees=>DegreeList];
	myVars = apply(myVars, z->apply(z,x->value(x)));
	numDegs = #myVars;
	myVarsData = {};
	i=0;
	while (i<numDegs) do (
		j=0;
		while(j<(#(myVars#i))) do(
			myVarsData = myVarsData | {myVars#i#j,Map_(i+1)_j};
			j=j+1;
		);
		i = i+1;
	);

	RelIdeal := ideal(0);
	Spar = S;

	c=1;
	while((c<bound) and (n_c>0)) do (
		Vect_c = transpose matrix{myVars#(c-1)};
		c=c+1;
	);
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
	b := 0;
	LengPa := 0;
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

	j=2;
	while ( (dim(Spar) >  dim(R)) or (isDomain(Spar) != true)) do ( --Create a list of relations
		
		Part_j = partitions(j);				--Uses the partitions of j to get relations coming from products of elements
		LengP = #(Part_j);   					--of lower degrees into the elements of degree j
		
		a:=0;
		AdmPart_j = {};					--Creates a list of admissable partitions, given that there are variables only
		while (a<LengP) do(					--in sufficiently small degree
			if(Part_j#a#0 < bound) then (
				AdmPart_j = AdmPart_j | {(Part_j)#a}; 
			);
			a=a+1;
		);
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 	
	
		LengP = #(AdmPart_j);

		a=0;

		TotMap_a = Map_((AdmPart_j)#a#0);			--Starts to create the map O_X(a_1 D) ** ... ** O_X(a_n D) \oplus O_X(j D) -> R		
		TotVect_a = Vect_((AdmPart_j)#a#0);
		b=1;
		LengPa = #((AdmPart_j)#a);
		while (b<LengPa) do (
			TotMap_a = TotMap_a ** Map_((AdmPart_j)#a#b);
			TotVect_a = TotVect_a ** Vect_((AdmPart_j)#a#b);
			b = b+1;
		);

		MapTot = TotMap_0;
		VectTot = TotVect_0;

		a=1;
		while(a<LengP) do ( 						--Relate each of the partitions in the form (a1>=a2>=...>=an) of j
			TotMap_a = Map_((AdmPart_j)#a#0);			--Starts to create the map O_X(a_1 D) ** ... ** O_X(a_n D) \oplus O_X(j D) -> R		
			TotVect_a = Vect_((AdmPart_j)#a#0);

			b=1;
			LengPa = #((AdmPart_j)#a);
			while (b<LengPa) do (
				TotMap_a = TotMap_a ** Map_((AdmPart_j)#a#b);
				TotVect_a = TotVect_a ** Vect_((AdmPart_j)#a#b);
				b = b+1;
			);

			MapTot = MapTot | TotMap_a;
			VectTot = VectTot || TotVect_a;
			a = a+1;
		);

		KerT = generators ker(MapTot);		

		NumCols = numColumns(KerT);
		e := 0;
			
		while (e < NumCols) do (
			L = flatten entries KerT_{e};
			if ((isVectScalar L) == true) then (
				L = convertScalarVect(S,L);
				Rel = sub((entries (matrix{L}*VectTot))#0#0,S);
				RelIdeal = trim(RelIdeal + ideal(Rel));
				Spar = S/RelIdeal;
			);
			e=e+1;
		); 

		j=j+1;
	);
	SectionRing = minimalPresentation Spar;
	SectionRing
);

isVectScalar = L -> (
	Ramb = ring (L#0); 
	all(L, z -> (degree(z) <= degree (sub(1, Ramb))) ) 
);

convertScalarVect = (newS, L) -> (apply(L, z->sub(z, newS)));

removalOfVariablesInternal = (L,n,Variables) -> (
	
);
