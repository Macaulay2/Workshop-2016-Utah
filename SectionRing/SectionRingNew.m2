newPackage( "SectionRing",
Version => "0.1", Date => "May 08 2016", Authors => {
     {Name=> "Andrew Bydlon",
     Email=> "thelongdivider@gmail.com",
     HomePage => "http://www.math.utah.edu/~bydlon/"
     }
}, --this file is in the public domain
Headline => "A package for computing the section ring of a Weil Divisor.", DebuggingMode => true, Reload=>true)

export{
	"globallyGenerated",
	"isMRegular",
	"mRegular",
	"sectionRing",
	"isVectScalar",
	"convertScalarVect"
}
needsPackage "Divisor"

dualToIdeal = method();
dualToIdeal(Ideal) := (I) -> (
	R := ring(I);
	M := module(I);
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
	j:=1;
	bool := true;
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
	local V;
	local G;
	local mRegularParticular;	
	V = variety(F);
	G = OO_V(1);
	mRegularParticular(F,G,m)
);

--isMRegularOO = method();
--isMRegularOO(CoherentSheaf,ZZ) := (F,m) -> (			--Tests if all higher cohomologies of F vanish
--	V := variety(F);
--	dV := dim(V);
--	j:=1;
--	bool = true;
--	while(j<(dV+1)) do (
--		if (bool == true) then(
--			bool = (HH^j(F)) == 0);
--		);
--		j = j+1;
--	);
--	bool
--);

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

