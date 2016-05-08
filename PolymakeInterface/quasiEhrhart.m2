restart
loadPackage("Polyhedra")




--Input: a polyhedron (perhaps rational)
--Output: a list of polynomials defining a quasi polynomail, index starts at 0 mod(period)
----method: polyhedra can already count lattice points. This program naively samples
------------the number of lattice points in a collection of dilations in order to perform
------------interpolation to create polynomials. 
------------It isn't quick. I'm fairly confident there is a quicker way to do this described
------------in a paper on the arxiv (which is implemented in latte?)
rationalEhrhart=method()
rationalEhrhart(Polyhedron) := (P) ->(
    M:=vertices(P);--pull the extreme vertices of the polytope
    R:=QQ[t];--create a ring for the quasi polynomials to live
    l=lcm(for r in flatten entries M list(denominator(r)));--pull denominators from the vertices of P
    d:=dim(P);
    print(concatenate("The rational simplex given must be stretched by a factor of ", toString(l), " in order to be integral"));
    print(concatenate("Therefore we must compute ",toString(l), " different polynomials. So this may take some time"));
    polyMatList:={};
    for k from 0 to l-1 do(     --it's quasi poly, so we need to interpolate polynomials for each k<lcm
	print(concatenate("Now computing quasipolynomial corresponding to mod(", toString(k),")"));
	Ak= new MutableMatrix from id_(QQ^(d+1));
	bk=matrix{{}};
	for j from 0 to d do ( ---I know each polynomial is degree d, so I  need d sample points to perform interpolation
	    if k+j==0 then bk=bk|(matrix{{1}});--the zeroth dilation of a polytope has one point in it.
	    if k+j!=0 then(
	    m:=k+j*l; --this is my dialation factor for the original polytope.
	    for i from 0 to d do(
		Ak_(j,i)=m^i;
		),
	    bk=bk|(matrix{{#(latticePoints(m*P))}});
	    ),
	    ),
	hardAk:=new Matrix from Ak;
    	polyMatList=join(polyMatList,{(hardAk^(-1))*transpose(bk)});
	),
    varMat:=matrix{{}};
    for i from 0 to d do(
	varMat=varMat|(matrix{{((R)_0)^i}});
	),
    ehrPolys:={};
    ehrPolys = join(ehrPolys,for A in polyMatList list((entries(varMat*A))#0#0));
    return(ehrPolys);
    )



   
P=convexHull(matrix{{0,0,1/3},{0,1/2,0}})--here's a rational polytope
ehrhart(P)--here's what polyhedra says the "ehrhart polynomial" should be.
--one easy way to see why this is inaccurate is to notice that if you dilate P by a factor of 1, you should get 1 lattice point. not 2
netList rationalEhrhart(P)

Q=convexHull(matrix{{1/2,1/2,3/2,3/2},{1/2,3/2,1/2,3/2}})
netList rationalEhrhart(Q)
