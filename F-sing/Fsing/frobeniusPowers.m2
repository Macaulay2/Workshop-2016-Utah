-- CHANGES IN ORDERS OF ARGUMENTS

--fastExp(N,f)  from (f,N)
--frobeniusPower(e,I) from (I,e)
--generalizedFrobeniusPower(e,I) from (I,e)



--*************************************************
--*************************************************
--This file is used for taking various types of 
--powers of ideals in characteristic p>0. 
--*************************************************
--*************************************************



--Computes powers of elements in char p>0, using that Frobenius is an endomorphism. If 
--N = N_0 + N_1 p + ... + N_e p^e, where 0 <= N_i < p, then this computes f^N as
--f^(N_0) (f^(N_1))^p ... (f^(N_e))^(p^e). 

fastExp = (N,f) ->
(
     p:=char ring f;
     E:=basePExp(p,N);
     product(#E, e -> (sum(terms f^(E#e), g -> g^(p^e))))
)

--------------------------------------------------------------------------------------------------------

--Outputs the p^e-th Frobenius power of an ideal, or the p^e-th (entry-wise) Frobenius power of a matrix.

frobeniusPower = method()

frobeniusPower(ZZ,Ideal) := (e,I) ->
(
     R := ring I;
     p := char R;
     G := I_*;
     if #G==0 then ideal(0_R) else ideal(apply(G, j -> fastExp(j,p^e)))
)

frobeniusPower(ZZ,Matrix) := (e,M) ->
(
    p:=char ring M;
    matrix apply(entries M,u -> apply(u, j -> fastExp(j,p^e)))
)

--------------------------------------------------------------------------------------------------------

--This is an internal function.  Given ideals I,J and a positive integer e, consider
--the following chain of ideals:
--K_1 = (I*J)^[1/p^e] and K_{s+1} = (I*K_s)^[1/p^e]
--This chain is ascending, and has the property that once two consecutive terms
--agree, the chain stabilizes.  This function outputs the stable ideal of this chain.

stableIdeal = (e,I,J) -> (
K := ideal(0_{ring I});
L := ethRoot(e,I*J);
while (isSubset(L,K)==false) do
    (
    K = L;              
    L = ethRoot(e,I*K);
    );
    trim K 
)

--------------------------------------------------------------------------------------------------------

--Outputs the generalized Frobenius power of an ideal; either the N-th Frobenius power of N/p^e-th one.

genFrobeniusPower = method(Options => {gfpStrategy => Naive});

--Computes the integral generalized Frobenius power I^[N]
genFrobeniusPower(ZZ,Ideal) := opts -> (N,I) -> 
(
     R := ring I;
     p := char R;
     G := first entries mingens I;
     if #G==0 then return ideal(0_R);
     if #G==1 then return ideal(fastExp(N,G#0));
     E := basePExp(p,N);
     product(#E, m->frobeniusPower(m,I^(E#m)))
)

--Computes the generalized Frobenius power I^[N/p^e]
genFrobeniusPower(ZZ,ZZ,Ideal) := opts -> (e,N,I) ->
(
     R := ring I;
     p := char R;
     G := first entries mingens I;
     if #G==0 then return ideal(0_R);
     rem := N % p^e;
     M := N // p^e;
     J := genFrobeniusPower(M,I);  --component when applying Skoda's theorem
     
    if opts.gfpStrategy==Safe then 
    (
	E := basePExp(p,rem);
	J*product(#E, m->ethRoot(e-m,I^(E#m)));  --step-by-step computation of generalized Frobenius power of I^[rem/p^e]
                                                                            --using the base p expansion of rem/p^e < 1
    )
    else J*ethRoot(e,genFrobeniusPower(rem,I))  --Skoda to compute I^[N/p^e] from I^[rem/p^e] 
 )

--Computes the generalized Frobenius power I^[t] for a rational number t 
genFrobeniusPower(QQ,Ideal) := (t,I) ->
(
    p = char ring I;
    L := divideFraction(p,t); 
    a := L#0; b := L#1; c := L#2;     --write t = a/(p^b*(p^c-1))
    if c==0 then return genFrobeniusPower(b,a,I)  --if c = 0, call simpler function
    	else 
	(
	    rem := a % (p^c-1);      
	    quot := a // (p^c-1);     
	    J := stableIdeal(c, genFrobeniusPower(rem,I), I);
	    ethRoot(b,genFrobeniusPower(quot,I)*J)
        );
)