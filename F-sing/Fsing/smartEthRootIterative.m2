--we should make an iterative implementation of this function at some point.
--here's a start...


--smartEthRoot(ZZ, List, List) := (e, exponentList, idealList) -> (
----the idealList, can take a list of ideals, a list of elements, or a mix of lists of ideals or elements
--    -- make a list of the number of generators of each ideal in idealList
--    minGensList = apply(idealList, jj -> (if (class jj === Ideal) then #(first entries mingens (jj)) else 1 ));
--
--    -- see what's the biggest power of p smaller than a/m where 'a' is in exponentList
--    -- and m is in minGensList. I.e. we want to find the largest e such that a >= mp^e
--    minGensLog = apply(minGensList, exponentList, (mm, aa) -> (
--        n = floorLog(p, aa//mm);
--        if (n > e) then e else n
--    ));
--
--    tripleList = sort apply(minGensLog, idealList, exponentList, (a,b,c) -> {a,b,c});
--
--
--    R := ring(idealList#0);
--    answer :=  ideal(1_R);
--    p := char(R);
--
--    for i from 0 to length(idealList) - 1 do (
--        answer = answer*(idealList#j)^(exponentList#j - p^minGensLog#j)
--    );
--    
--    j := 0;
--    for i from 0 to e do (
--        if i == tripleList#j#0 then (
--            answer = answer*
--            j = j+1;
--        );
--    );
--
--);


