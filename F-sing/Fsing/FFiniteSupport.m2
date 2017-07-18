
--- Given matrices A, B with target R^alpha find all v\in R^alpha such that B v \in Image A
--- by finding partial syzygies
matrixColon= (A, B) ->(
assert(target(A)==target(B));
m:=rank source B;
M:=B | A;
S:=syz(M);
S^(toList(0..m-1))
)


--- Given a generating morphism U:coker(A) -> F(coker A), compute a generating root U:coker(L) -> F(coker L)
generatingRoot= (A,U) ->(
	R:=ring(A);
	L:=A;
	alpha:=rank target A;
	LL:=transpose matrix{toList(alpha:0_R)};
	while ((( L)%( LL))!=0) do
	{
		LL=L;
		L=L | matrixColon(frobenius( L ),U);
		L=mingens image L;
---		print("=================================================================");
---		print(L);
	};
	L
)


--- Given a generating morphism U:coker(A) -> F(coker A), compute the support of the F-module
FFiniteSupport= (A,U) ->(
	R:=ring(A);
	alpha:=rank source U;
	LL:=id_(R^alpha);
	L:=ethRoot(1, U*LL);
	while (((gens image LL) %( (gens image L)|A))!=0) do
	{
		LL=L;
---		L=mingens image U*L;
		L=U*L;
		L=ethRoot(1,L);
---?		L=mingens image L;
---		print("=================================================================");
---		print(L);
	};
	answer:=prune subquotient(A | L, A);
	mingens radical annihilator answer
)

