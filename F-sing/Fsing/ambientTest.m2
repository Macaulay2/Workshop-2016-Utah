findTestElementAmbient = method()
randomSubset = method()

findTestElementAmbient(Ring) := (R) ->
(
	I = ideal R;
	n = #gens R - dim R;
	M = jacobian I;
	r = rank target M;
	c = rank source N;
	testEle = ideal(sub(0,ambient R));
	while(isSubset(testEle, I))
	do(
	   testEle = minors(n,M, First =>{randomSubset(r,n),randomSubset(c,n)}, Limit =>1);
	);
	testEle
)

randomSubset(ZZ,ZZ) := (m,n) ->
(
	L = for i from 0 to m-1 list i;
	for i from 0 to m-n-1 do (L = delete(L#(random(0,m-1-i)),L));
	L
)


	