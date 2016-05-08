newPackage(
	"GLmnReps",
    	Version => "1.0", 
    	Date => "May 8, 2016",
    	Authors => {
	     {Name => "Michael Perlman", Email => "mperlman@nd.edu", HomePage => "https://www3.nd.edu/~mperlman/"},
	     {Name => "Claudiu Raicu", Email => "craicu@nd.edu", HomePage => "https://www3.nd.edu/~craicu/"}
	     },
    	Headline => "Representations of gl(mn) and syzygies",
	PackageExports => {"SchurRings","BGG"},
    	DebuggingMode => false
    	)
    
export{"hhh","charL", "charK", "bettiPart"}

charE = method()
charE(ZZ,ZZ) := (m,n) -> (
    parsmn := flatten for i from 0 to m*n list select(partitions i,x->#x <= n and (try(x#0 <= m) else true));
    apply(parsmn,par -> {toList conjugate par,toList par})   
    )
--charE = memoize charE

charK = method()
charK(List,List,SchurRing,SchurRing) := (mn,ll,S,T) -> (
    cE := charE(mn#0,mn#1);
    chrE := sum apply(cE,x -> S_(x#0) * T_(x#1));
    S_ll * T_ll * chrE
    )
charK(ZZ,ZZ,List) := (m,n,ll) -> (
    T := schurRing(getSymbol"t000t000t",n);
    S := schurRing(T,getSymbol"s000s000s",m);
    chrK := charK({m,n},ll,S,T);
    flatten apply(listForm(chrK),ter -> apply(listForm(ter#1), x -> {ter#0,x#0,x#1}))
    )
--charE = memoize charE

delta := k -> for i from 0 to k-1 list k-1-i
delta = memoize delta
dotact := (per,mu,k) -> (dk := delta k; (for i from 0 to k-1 list mu#(per#i) + dk#(per#i))-dk)
dotact = memoize dotact
leng := (per,k) -> (lng := 0; for i from 0 to k-2 do for j from i+1 to k-1 do if per#i>per#j then lng = lng + 1;lng)
leng = memoize leng

generateMus := (lam,k) ->
(
    prs := {};
    for i from 0 to k-2 do
    (
    	noexit := true;
    	for j from i+1 to k-1 do if (lam#i - lam#j < j - i and noexit) then prs = prs | {(i,j)} else noexit = false;
    	);
    perms := select(permutations k,per -> (include := true; for c in prs do if per#(c#0)>per#(c#1) then include = false; include));
    perms = apply(perms, per -> apply(sort(apply(per,i->(per#i,i))),j->j#1));
    musleng := for per in perms list (pl := dotact(per,lam,k); lst := last pl; (leng(per,k),reverse(for i from 0 to k-1 list (lst = max(lst,pl#(k-1-i));lst))));
    musleng / (x->x#1)
    )
--generateMus = memoize generateMus

charL = method()
charL(ZZ,ZZ,List) := (m,n,lam) -> (
    lam = lam | splice{n-#lam:0};
    mus := generateMus(lam,n);
    MAX := lam + splice{n:m};
    rez := 0;
    S := schurRing(getSymbol"s000s000s",m);
    T := schurRing(S,getSymbol"t000t000t",n);
    for mu in mus do 
        for ll in mu..MAX do if ll == reverse sort(ll) then 
	    rez = rez + (-1)^(sum(ll-lam)) * charK({m,n},ll,S,T);
    selsmall := select(listForm rez,x->all(splice{0..#(x#0)-1},i-> (x#0#i <= MAX#i)));
    flatten apply(selsmall,x -> apply(listForm(x#1),ter -> {ter#0,x#0,ter#1}))
    )
--charL = memoize charL

bettiPart = method()
bettiPart(ZZ,ZZ,ZZ,List) := (n,m,p,lam) -> (
     r := local r;
     s := local s;
     x := local x;
     R := schurRing(r,n);
     T := schurRing(s,m);
     conjlam := toList conjugate( new Partition from lam);
     d := dim r_lam;
     e := dim s_lam;
     kk := ZZ/p;
     S := kk[x_(1,1)..x_(n,m)];
     M := genericMatrix(S,m,n);
     lis := for i from 0 to d*e-1 list
     (
    A := random(kk^m,kk^m);
    B := random(kk^n,kk^n);
    N := A * M * B;
    product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1})
     );
    J := ideal lis;
    I := ideal mingens J;
    if (numgens I != e*d) then error"wrong number of generators"
    else betti res I --needs to be replaced by next two lines
    --else F:= res(I, FastNonminimal => true);
   -- betti(F, Minimize => true)
    )

bettiPart(ZZ,ZZ,List) := (n,m,lam) -> (
    bettiPart(n,m,32003,lam)
    )

g = method()
g(ZZ,ZZ,List,ZZ) := (n,m,lam,d) -> (
    F:=select(charL(n,m,lam), x-> sum(x#0) === d);
    t:= local t;
    s:= local s;
    T1:= schurRing(t,n);
    T2:= schurRing(s,m);
    tempList:=apply(F, x-> dim(T1_(x#0))*dim(T2_(x#1))*x#2);
    sum tempList
    )

h = method()
h(ZZ,ZZ,List) := (n,m,lam) -> (
    G:= apply(charL(n,m,lam),x->sum(x#0));
    hlist:= for i from sum lam to max G list
    (
    g(n,m,lam,i)
    );
    hlist
    )

  
end

restart
uninstallPackage"GLmnReps"
installPackage"GLmnReps"
debug needsPackage"GLmnReps"

charL(3,2,{2,1})
charK(3,2,{1})
charE(4,3)
h(3,3,{3,1})
g(3,3,{3,1},5)

N=matrix{{1,2,3},{4,5,6}}
N_{0,1}^{0,1}




g := method()
g(ZZ,ZZ,List,ZZ) := (n,m,lam,d) -> (
    F:=select(charL(n,m,lam), x-> sum(x#0) === d);
    t:= local t;
    s:= local s;
    T1:= schurRing(t,n);
    T2:= schurRing(s,m);
    tempList:=apply(F, x-> dim(T1_(x#0))*dim(T2_(x#1))*x#2);
    sum tempList
    )

h := method()
h(ZZ,ZZ,List) := (n,m,lam) -> (
    G:= apply(charL(n,m,lam),x->sum(x#0));
    hlist:= for i from sum lam to max G list
    (
    g(n,m,lam,i)
    );
    hlist
    )


I = ideal mingens J;
if (numgens I != d^2) then error"wrong number of generators"

--F = resolution(I,DegreeLimit => 5)
--F = resolution I
F = res(I, FastNonminimal => true)
betti(F, Minimize => true)
end

restart
a = 4
b = 2
time load"syz12.m2"
betti F

select(charL(3,3,{3,1}), x -> sum(x#0) == 5)
T = schurRing(t,3)
sum apply(o22, x -> dim(T_(x#0)) * dim(T_(x#1)) * x#2)
apply(charL(3,3,{3,1}),x-> sum(x#0))
max oo
min ooo
viewHelp List
