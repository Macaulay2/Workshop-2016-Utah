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

--Outputs the generalized Frobenius power of an ideal; either the N-th Frobenius power of N/p^e-th one.

genFrobeniusPower = method()

genFrobeniusPower(ZZ,Ideal) := (N,I) ->
(
     p := char ring I;
     E := basePExp(p,N);
     product(#E, m->frobeniusPower(m,I^(E#m)))
)

genFrobeniusPower(ZZ,MonomialIdeal) := (N,I) ->
(
     p := char ring I;
     E := basePExp(p,N);
     product(#E, m->frobeniusPower(m,I^(E#m)))
)

genFrobeniusPower(ZZ,ZZ,Ideal) := (e,N,I) ->
(
    ethRoot(e,genFrobeniusPower(N,I))
)

genFrobeniusPower(ZZ,ZZ,MonomialIdeal) := (e,N,I) ->
(
    ethRoot(e,genFrobeniusPower(N,I))
)
