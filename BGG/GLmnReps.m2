newPackage(
	"GLmnReps",
    	Version => "1.0", 
    	Date => "May 8, 2016",
    	Authors => {
	     {Name => "Michael Perlman", Email => "a", HomePage => "https://www3.nd.edu/~mperlman/"},
	     {Name => "Claudiu Raicu", Email => "craicu@nd.edu", HomePage => "https://www3.nd.edu/~craicu/"}
	     },
    	Headline => "Representations of gl(mn) and syzygies",
	PackageExports => {"SchurRings","BGG"},
    	DebuggingMode => false
    	)
    
export{"charE", "charL", "charK"}

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
  
end

restart
uninstallPackage"GLmnReps"
installPackage"GLmnReps"
charL(3,2,{2,1})
charK(3,2,{1})
charE(4,3)
