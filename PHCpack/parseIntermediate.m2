f = (output,numsols) -> (
    L := get output;
    rgx = "\\*{5} +path([[:digit:]])\s+\\*+$";
    start := (regex(rgx,L))#0#0;
    linesL := lines substring(start,L);
    solsize := 0;
    for i from 1 to #linesL do (
        if match("^\\*",linesL#i) then (
            solsize = i;
            break;
            ););
    chunksize := numsols*solsize;
    chunks := {};
    while match(rgx,linesL#0) do (
        chunks = append(chunks,take(linesL,chunksize));
        linesL = drop(linesL,chunksize);
    );
    chunks = for chunk in chunks list (
        for l in chunk list (
            replace(rgx,"solution \\1 :",l);
        );
    );
    
)