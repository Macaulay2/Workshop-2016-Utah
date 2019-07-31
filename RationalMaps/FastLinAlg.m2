newPackage( "FastLinAlg",
Version => "0.1", Date => "July 24th, 2019", Authors => {
    {Name => "Boyana Martinova",
    Email=> "u1056124@utah.edu"
    },
    {Name => "Marcus Robinson",
    Email => "robinson@math.utah.edu",
    HomePage => "http://www.math.utah.edu/~robinson"
    },
    {Name => "Karl Schwede",
    Email=> "schwede@math.utah.edu",
    HomePage=> "http://www.math.utah.edu/~schwede"
    },
    {Name => "Yuhui (Wei) Yao",
    Email=> "yuhuiyao4ever@gmail.com"
    }
}, --this file is in the public domain
Headline => "A package for faster linear algebra operations.", DebuggingMode => true, Reload=>true)
export{
  "selectSmallestTerms",
  "chooseSubmatrixSmallestDegree",
  "chooseSubmatrixLargestDegree",
--  "chooseMinorLargestDegree",
--  "chooseMinorSmallestDegree",
  --"randomMaxPosition",
  --"randomMinPosition",
  --"randomMinPositions",
  --"replaceSmallestTerm",
  --"replaceLargestTerm",
  "chooseRandomSubmatrix",
  --"locationToSubmatrix",
  "reorderPolynomialRing",
  --"getDetOfSubmatrix",
--  "nonzeroEntries",
  "Rn",
  "RnReductionP",
  "chooseGoodMinors",
  "projDim",
  "isRankAtLeast",
  "getSubmatrixOfRank",
  "recursiveMinors",
  --options to export
  "DetStrategy", --to pass on to calls of determinant
  "MinDimension", --option for projDim
  "Random",
  "LexLargest",
  "LexSmallest",
  "LexSmallestTerm",
  "RandomNonzero",
  "GRevLexLargest",
  "GRevLexSmallest",
  "GRevLexSmallestTerm",
  "MaxMinors",
  "Rank", --a value for Strategy in isRankAtLeast
  "MutableSmallest", --these should be protected but maybe not exported
  "MutableLargest", --these should be protected but maybe not exported
  "Threads",
  "MinorsCache",
  --"RecursiveMinors",
  --premade stratgies
  "Recursive",
  "StrategyDefault",
  "StrategyGRevLexSmallest",
  "StrategyLexSmallest",
  "StrategyRandom",
  "StrategyCurrent"
}


--***********************************
--*****Control options universally***
--***********************************
StrategyDefault := new HashTable from {
    LexLargest => 0,
    LexSmallestTerm => 40,
    LexSmallest => 40,
    GRevLexSmallestTerm => 0,
    GRevLexSmallest => 0,
    GRevLexLargest => 0,
    Random => 20,
    RandomNonzero => 0
};

StrategyCurrent = new MutableHashTable from StrategyDefault;

StrategyGRevLexSmallest = new HashTable from {LexLargest=>0, LexSmallestTerm => 0, LexSmallest=>0, GRevLexSmallestTerm => 50, GRevLexSmallest => 50, GRevLexLargest=>0,Random=>0,RandomNonzero=>0};

StrategyLexSmallest = new HashTable from {LexLargest=>0, LexSmallestTerm => 50, LexSmallest=>50, GRevLexSmallestTerm=>0, GRevLexSmallest=>0, GRevLexLargest=>0,Random=>0,RandomNonzero=>0};

StrategyRandom = new HashTable from {LexLargest=>0, LexSmallestTerm => 0, LexSmallest=>0, GRevLexSmallestTerm=>0, GRevLexSmallest=>0, GRevLexLargest=>0,Random=>100,RandomNonzero=>0};

optRn := {
    Verbose => false,
    MaxMinors => null,
    Strategy => StrategyDefault,
    DetStrategy => Cofactor
};

optInternalChooseMinor := {
    MutableSmallest => null,
    MutableLargest => null,
    Strategy => StrategyDefault,
    Verbose => false
};

optProjDim := {
    MinDimension => 0,
    Verbose => false,
    Strategy => StrategyDefault,
    DetStrategy => Cofactor,
    MaxMinors => null
};

optIsRankAtLeast :=  {
    Verbose => false,
    DetStrategy => Rank,
    MaxMinors => null,
    Strategy => StrategyDefault
};

optChooseGoodMinors := {
    Verbose => false,
    Strategy => StrategyDefault,
    DetStrategy=>Cofactor
};


--************************************
-------------------------------------
--Choosing submatrices---------------
-------------------------------------
--************************************

chooseSubmatrixLargestDegree = method(Options => {});
--------------------------------------
--This command tries to find a n1xn1 submatrix with
--determinant of biggest possible degree.
--It does it by picking the largest element in the matrix,
--then removing that row and column, and repeating the process
--n1 times.
--It returns the list of rows and columns that determine the submatrix.
--------------------------------------


chooseSubmatrixLargestDegree(ZZ, Matrix) := opts -> (n1, M1) -> (
    rCt := numRows M1;
    cCt := numColumns M1;
    i := 0;
    j := 0;
    curM1 := M1;
    --degreeMatrx := matrix apply(entries M1, l1 -> apply(l1, i->degree i));
    curList := null;
    curMax := null;
    curCol := null;
    curRow := null;
    newCurCol := null;
    newCurRow := null;

    returnRowList := {};--maybe these should be mutable lists, we'll do that later...
    returnColList := {};
    --keepRowList := {};
    --keepColList := {};

    while (i < n1) do (
        --print concatenate("in loop, i =", toString(i));
        curList = flatten entries curM1;
        curMax = randomMaxPosition curList;
        --print (curList#curMax);
        curRow = curMax // cCt;
        curCol = curMax % cCt;
        curM1 = submatrix'(curM1, {curRow}, {curCol});
        rCt = rCt - 1;
        cCt = cCt - 1;
        --we need to adjust based on submatrices
        newCurRow = curRow;
        j = (#returnRowList) - 1;
        --and for columns
        newCurCol = curCol;
        j = (#returnColList) - 1;
        returnRowList = append(returnRowList, newCurRow);
        returnColList = append(returnColList, newCurCol);
        i = i + 1;
    );
    --print returnRowList;
    returnRowList = new MutableList from returnRowList;
    returnColList = new MutableList from returnColList;
    j = n1-1;
    --when we removed rows and columns, the numbering probably got off.
    --We are fixing that in the following two loops.
    while (j >= 0) do ( --compare entry j
        newCurRow = returnRowList#j;
        i = j-1;
        --print "here1";
        while  (i >= 0) do (--compare entry j to previous terms in the list
            --print "here2";
            if (newCurRow >= returnRowList#i) then (newCurRow = newCurRow + 1;);
            i = i-1;
        );
        returnRowList#j = newCurRow;
        j = j-1;
    );
    j = n1-1;
    while (j >= 0) do ( --compare entry j
        newCurCol = returnColList#j;
        i = j-1;
        while  (i >= 0) do (--compare entry j to previous terms in the list
            if (newCurCol >= returnColList#i) then (newCurCol = newCurCol + 1;);
            i = i-1;
        );
        returnColList#j = newCurCol;
        j = j-1;
    );
    return {new List from returnRowList, new List from returnColList};
);

--the following commands use chooseSubmatrixLargestDegree
--to return the actual submatrix.

chooseMinorLargestDegree = method(Options => {});
chooseMinorLargestDegree(ZZ, Matrix) := o -> (n1, M1) -> (
    bestDegree := chooseSubmatrixLargestDegree(n1, M1);
    minorRowList := bestDegree#0;
    minorColList := bestDegree#1;
    return (M1^minorRowList)_minorColList;
)
chooseMinorLargestDegree(ZZ, MutableMatrix) := o -> (n1, M1) -> (
    return new Matrix from chooseMinorLargestDegree(n1, matrix(M1), o);
  );

  --------------------------------------
  --This command finds an n1xn1 submatrix
  --This selects it randomly, but it makes sure there is a nonzero entry in each row and column.
  --------------------------------------


chooseRandomNonzeroSubmatrix = method(Options => {});

chooseRandomNonzeroSubmatrix(ZZ, Matrix) := opts -> (n1, M1) -> (
  rCt := numRows M1;
  cCt := numColumns M1;
  i := 0;
  j := 0;
  curM1 := M1;
  curList := null;
  curEntry := null;
  entryList := {};
  curCol := null;
  curRow := null;
  newCurCol := null;
  newCurRow := null;

  returnRowList := {};--maybe these should be mutable lists, we'll do that later...
  returnColList := {};

  while (i < n1) do (
      --print concatenate("in loop, i =", toString(i));
      --curList = flatten entries curM1;
      entryList = entries transpose matrix nonzeroEntries(curM1);
      curEntry = entryList#(random(#entryList));
      --print (curList#curMax);
      curRow = curEntry#0;
      curCol = curEntry#1;
      curM1 = submatrix'(curM1, {curRow}, {curCol});
      rCt = rCt - 1;
      cCt = cCt - 1;
      --we need to adjust based on submatrices
      newCurRow = curRow;
      --print curRow;
      j = (#returnRowList) - 1;
      --and for columns
      newCurCol = curCol;
      j = (#returnColList) - 1;
      returnRowList = append(returnRowList, newCurRow);
      returnColList = append(returnColList, newCurCol);
      i = i + 1;
  );
  returnRowList = new MutableList from returnRowList;
  returnColList = new MutableList from returnColList;
  j = n1-1;
  --when we removed rows and columns, the numbering probably got off.
  --We are fixing that in the following two loops.
  while (j >= 0) do ( --compare entry j
      newCurRow = returnRowList#j;
      i = j-1;
      --print "here1";
      while  (i >= 0) do (--compare entry j to previous terms in the list
          --print "here2";
          if (newCurRow >= returnRowList#i) then (newCurRow = newCurRow + 1;);
          i = i-1;
      );
      returnRowList#j = newCurRow;
      j = j-1;
  );
  j = n1-1;
  while (j >= 0) do ( --compare entry j
      newCurCol = returnColList#j;
      i = j-1;
      while  (i >= 0) do (--compare entry j to previous terms in the list
          if (newCurCol >= returnColList#i) then (newCurCol = newCurCol + 1;);
          i = i-1;
      );
      returnColList#j = newCurCol;
      j = j-1;
  );
  return {new List from returnRowList, new List from returnColList};
);

--------------------------------------
--This command takes in a matrix, and replaces the zeros with high degree polynomials
--------------------------------------


replaceZeros= method(Options=>{});

replaceZeros(Matrix):= Matrix => o->(M2) -> (
    Mute := mutableMatrix M2;
    m := numRows M2;
    n := numColumns M2;
    M2ent := flatten entries M2;
    largeDeg := 1;
    if (#M2ent > 0) then largeDeg = 2*max(largeDeg, max(flatten apply(flatten entries M2, z->degree z)));

    largeGen := (last gens ambient ring M2)^(2*largeDeg+2);

    i := 0;
    while (i<n*m) do (
        Row := i//n;
        Col := i % n;
        if ((flatten entries(M2_{Col}))#Row==0)
        then Mute_(Row, Col)=largeGen;
        i=i+1;
    );
    unMute := matrix Mute;
    return unMute;
);

selectSmallestTerms = method(Options=>{});
--this function takes a matrix and replaces each entry with the smallest monomial term
--this is useful when we want to find the minor with the smallest term

selectSmallestTerms(Matrix) := Matrix => o->(M2) -> (
    entryList := entries M2;
    newEntryList := apply(entryList, myRow -> apply(myRow, f2 -> min terms f2));
    matrix newEntryList
);

chooseSubmatrixSmallestDegree = method(Options=>{});
    --------------------------------------
    --This command tries to find a n1xn1 submatrix with
    --determinant with a smallest term
    --Frequently this is used in concert with selectSmallestTerms which replaces
    --the entries of a matrix with their smallest term.
    --It does this by choosing the smallest (with respect to the monomial order)
    --term of the matrix, removing that row and column, and repeating
    --n1 times.
    --It returns the list of rows and columns that determine the submatrix
    --------------------------------------

chooseSubmatrixSmallestDegree(ZZ, Matrix) := o -> (n1, M3) -> (
          --M1 := new Matrix replaceZeros(M3);
          M1 := new Matrix from M3;
          rCt := numRows M1;
          cCt := numColumns M1;
          i := 0;
          j := 0;
          curM1 := M1;
          curList := null;
          curMax := null;
          curCol := null;
          curRow := null;
          newCurCol := null;
          newCurRow := null;

          returnRowList := {};
          returnColList := {};

          while (i < n1) do (
              --print concatenate("in loop, i =", toString(i));
              curList = flatten entries curM1;
              curMax = randomMinPosition curList;
              curRow = curMax // cCt;
              curCol = curMax % cCt;
              curM1 = submatrix'(curM1, {curRow}, {curCol});
              rCt = rCt - 1;
              cCt = cCt - 1;
              newCurRow = curRow;
              j = (#returnRowList) - 1;
              newCurCol = curCol;
              j = (#returnColList) - 1;
              returnRowList = append(returnRowList, newCurRow);
              returnColList = append(returnColList, newCurCol);
              i = i + 1;
          );
          returnRowList = new MutableList from returnRowList;
          returnColList = new MutableList from returnColList;
          j = n1-1;
          while (j >= 0) do (
              newCurRow = returnRowList#j;
              i = j-1;
              while  (i >= 0) do (
                  if (newCurRow >= returnRowList#i) then (newCurRow = newCurRow + 1;);
                  i = i-1;
              );
              returnRowList#j = newCurRow;
              j = j-1;
          );
          j = n1-1;
          while (j >= 0) do (
              newCurCol = returnColList#j;
              i = j-1;
              while  (i >= 0) do (
                  if (newCurCol >= returnColList#i) then (newCurCol = newCurCol + 1;);
                  i = i-1;
              );
              returnColList#j = newCurCol;
              j = j-1;
          );
          return {new List from returnRowList, new List from returnColList};
      );

--these functions just return the submatrix the previous function identifies.

chooseMinorSmallestDegree = method(Options => {});
chooseMinorSmallestDegree(ZZ, MutableMatrix) := o -> (n1, M1) -> (
    return chooseMinorSmallestDegree(n1, matrix(M1), o);
  );

  chooseMinorSmallestDegree(ZZ, Matrix) := o -> (n1, M1) -> (
      bestDegree := chooseSubmatrixSmallestDegree(n1, M1);
      minorRowList := bestDegree#0;
      minorColList := bestDegree#1;
      return (M1^minorRowList)_minorColList;
  );

--the following command chooses a *random* maximum element of a list
randomMaxPosition = method(Options=>{});
randomMaxPosition(List) := o -> (L1) -> (
    newList := apply(#L1, i -> {L1#i, random((#L1)^2), i});
    newList2 := random(newList);
    j := maxPosition(newList2);
    return ((newList2#j)#2);
);

--the following command chooses a *random* minimum element of a list
randomMinPosition = method(Options=>{});
randomMinPosition(List) := o -> (L1) -> (
    newList := apply(#L1, i -> {L1#i, random((#L1)^2), i});
    newList2 := random(newList);
    j := minPosition(newList2);
    return ((newList2#j)#2);
);

--the following command chooses n1 *random* minumum element of a list
randomMinPositions = method(Options=>{});
randomMinPositions(ZZ, List) := o -> (n1, L1) -> (
    newList := apply(#L1, i -> {L1#i, random((#L1)^2), i});
    newList2 := random(newList);
    newList3 := sort(newList2);
    return apply(take(n1, newList3), z -> z#0);
);

--this function replaces one of the smallest terms in the matrix by a larger term,
-- hopefully then allowing us to identify the next smallest term
-- in fact, it chooses one of the entries we picked when
--identifying our smallest matrix and randomly increases it.
replaceSmallestTerm= method(Options=>{});
replaceSmallestTerm(List, MutableMatrix) := opts -> (submatrixS, M1) -> (
  ambR:= ambient ring(M1);
  rowListS := submatrixS#0;
  colListS := submatrixS#1;
  mutedSM := M1;
  myRand := random(#rowListS);
  moddedRow := rowListS#(myRand);
  moddedCol := colListS#(myRand);
  val := (M1_(moddedRow, moddedCol))*((gens ambR)#(random(#(gens ambR))));
  mutedSM_(moddedRow, moddedCol) = val;
  return mutedSM;
  );

--this works similarly to the above, but it doesn't choose a random term among
--the largest, it really picks the biggest one.   First we drop the extra terms,
--then we replace the term by a smaller degree term.
replaceLargestTerm= method(Options=>{});
replaceLargestTerm(List, MutableMatrix):= opts-> (submatrixL, M1)->(
  --submatrixL := chooseSubmatrixLargestDegree(FR, M1);
  rowListL := submatrixL#0;
  colListL := submatrixL#1;
  largeSubmatrix :=(M1^rowListL)_colListL;
  termLists := flatten entries monomials(M1_(rowListL#0, colListL#0));
  local myRand;

  muted := M1;
  if (#termLists > 1) then ( --drop the smallest degree term basically, and sum what's left
      minTerm := randomMinPosition(termLists);
      muted_(rowListL#0, colListL#0) = sum(drop(termLists, {minTerm, minTerm}));
  )
  else if (#termLists == 1) then (
      minTerm = new MutableList from factor(termLists#0);
      if (#minTerm > 0) then (
          myRand = random(#minTerm);
          minTerm#myRand = new Power from {(minTerm#myRand)#0, (minTerm#myRand#1) - 1};
          muted_(rowListL#0, colListL#0) = value (new Product from minTerm);
      )
      else (
          muted_(rowListL#0, colListL#0) = sub(0, ring M1);
      );
      j := minPosition
  )
  else(
      muted_(rowListL#0, colListL#0) = 0;
      );
  --replaced := new Matrix from muted;
  return muted;
  );

--this just chooses a random submatrix
chooseRandomSubmatrix = method(Options=>{});

chooseRandomSubmatrix(ZZ, Matrix) := opts -> (n1, M1) ->
(
    rowL := random(toList(0..(numRows M1 - 1)));
    colL := random(toList(0..(numColumns M1 - 1)));
    rowL = sort take(rowL, n1);
    colL = take(colL, n1);
    return {rowL, colL};
    )


--this takes in a list of locations in a matrix determining a submatrix
--and returns the sorted list of rows and columns.
locationToSubmatrix = method();

locationToSubmatrix(List) := (L1) ->
(  {sort(L1#0), sort(L1#1)}  );


--in our functions, we choose smallest elements in our matrices, with
--respect to lots of different monomial orders, so we need code to
--change the polynomial ring order, randomly
reorderPolynomialRing = method(Options=>{});

reorderPolynomialRing(Symbol, PolynomialRing) := opts-> (myOrder, R1) -> (
        coeff := coefficientRing R1;
        genList := generators R1;
        newGenList := random(genList);
        return coeff[newGenList, MonomialOrder => myOrder];
);

--this just gets the det of a submatrix, it has some strategy options
getDetOfSubmatrix = method(Options=>{DetStrategy=>Cofactor});

getDetOfSubmatrix(Matrix, List) := opts -> (M2, submat) -> (
    mRowListS := submat#0;
    mColListS := submat#1;
    smallestSubmatrix :=(M2^mRowListS)_mColListS;

    if (opts.DetStrategy === Recursive) then (
        J := recursiveMinors(#(submat#0), smallestSubmatrix);
        return first first entries gens J;
    );

    return determinant(smallestSubmatrix, Strategy=>opts.DetStrategy);-- take determinant
);

--this function is used in the choose nonzero matrix bit
nonzeroEntries = method(Options=>{});
nonzeroEntries (Matrix):= opts -> (M1) ->(
  entryList := flatten entries M1;
  numTerms := #entryList;

  n := numColumns M1;

  i:=0;
  thatTerm:= null;
  nonzeroRowList := {};
  nonzeroColList := {};

  nzRowNum := null;
  nzColNum := null;
  while (i<numTerms) do (
    if not (entryList#i ==0) then (
      nzRowNum = i//n;
      nzColNum = i%n;
      nonzeroRowList = append(nonzeroRowList, nzRowNum);
      nonzeroColList = append(nonzeroColList, nzColNum);
      );
    i=i+1;
    );
    return {nonzeroRowList, nonzeroColList};
  );



--this function checks Rn via reduction mod p
RnReductionP = method(Options=>optRn);
RnReductionP(ZZ, Ring, ZZ):= opts -> (n1, R1, p)-> (
    ambR := ambient R1;
    genList := generators(ambR);
    ambRing := ZZ/p[genList];

    Id := sub(ideal R1, ambRing);
    return Rn(n1, ambRing/Id, opts);
);

--an internal function which chooses a minor, based on chance
internalChooseMinor = method(Options=>optInternalChooseMinor);

internalChooseMinor(ZZ, Ring, Matrix, Matrix) := opts -> (minorSize, ambR, nonzeroM, M1) -> (
    --if (opts.Verbose or (debugLevel > 0)) then print "internalChooseMinor: starting.";
    totalPercent :=  opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest + opts.Strategy#GRevLexSmallestTerm + opts.Strategy#GRevLexLargest + opts.Strategy#Random + opts.Strategy#RandomNonzero;
    myRandom := random(totalPercent);
    submatrixS1 := null;
    local R2;
    local f;
    mutM2 := opts.MutableSmallest;
    mutM1 := opts.MutableLargest;
    local M2;
    if (myRandom < opts.Strategy#LexSmallest) then (
        R2 = reorderPolynomialRing(Lex, ambR); --do the same with respect to a Lex ordering
        f = map(R2, ambR);
        M2 = f(nonzeroM);
        submatrixS1 = chooseSubmatrixSmallestDegree(minorSize, M2);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm) then (
        R2 = reorderPolynomialRing(Lex, ambR); --do the same with respect to a Lex ordering
        f = map(R2, ambR);
        M2 = f(nonzeroM);
        submatrixS1 = chooseSubmatrixSmallestDegree(minorSize, selectSmallestTerms M2);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest) then
    (
        R2 = reorderPolynomialRing(Lex, ambR); --do the same with respect to a Lex ordering
        f = map(R2, ambR);
        M2 = f(nonzeroM);
        submatrixS1 = chooseSubmatrixLargestDegree(minorSize, M2);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest) then
    (
        R2 = reorderPolynomialRing(GRevLex, ambR);
        f = map(R2, ambR);
        M2 = f(matrix mutM2); --put the matrix in the ring with the new order
        submatrixS1 = chooseSubmatrixSmallestDegree(minorSize, M2);
        mutM2 = replaceSmallestTerm(submatrixS1, mutM2);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest + opts.Strategy#GRevLexSmallestTerm) then
    (
        R2 = reorderPolynomialRing(GRevLex, ambR);
        f = map(R2, ambR);
        M2 = f(matrix mutM2); --put the matrix in the ring with the new order
        submatrixS1 = chooseSubmatrixSmallestDegree(minorSize, selectSmallestTerms M2);
        mutM2 = replaceSmallestTerm(submatrixS1, mutM2);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest + opts.Strategy#GRevLexSmallestTerm + opts.Strategy#GRevLexLargest ) then (
        R2 = reorderPolynomialRing(GRevLex, ambR);
        f = map(R2, ambR);
        M1 = f(matrix mutM1); --put the matrix in the ring with the new order
        submatrixS1 = chooseSubmatrixLargestDegree(minorSize, M1);
        mutM1 = replaceLargestTerm(submatrixS1, mutM1);
        --this needs to be written
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest + opts.Strategy#GRevLexSmallestTerm + opts.Strategy#GRevLexLargest + opts.Strategy#Random) then (
        submatrixS1 = chooseRandomSubmatrix(minorSize, M1);
    )
    else if (myRandom < opts.Strategy#LexSmallest + opts.Strategy#LexSmallestTerm + opts.Strategy#LexLargest + opts.Strategy#GRevLexSmallest + opts.Strategy#GRevLexSmallestTerm + opts.Strategy#GRevLexLargest + opts.Strategy#Random + opts.Strategy#RandomNonzero) then (
        submatrixS1 = chooseRandomNonzeroSubmatrix(minorSize, M1);
    );
    --if (opts.Verbose or (debugLevel > 0)) then print "internalChooseMinor: finished.";
    return submatrixS1;
)



Rn = method(Options=>optRn);

Rn(ZZ, Ring) := opts -> (n1, R1) -> (
    ambR := ambient R1;
    Id := ideal R1;
    M1 := sub(jacobian Id, ambR);
    numberRelations := numColumns(M1);
    n := numgens ambR;
    r := dim R1;
    fullRank := n-r;
    myRand := 0;
    Q := null; C := null; C2 := null; g := null;
    --1/0;

    possibleMinors := binomial(n, fullRank)*binomial(numberRelations, fullRank);
    minNumberToCutDown := n1+1;
    local numberOfMinorsCompute;
    if (opts.MaxMinors === null) then (
        numberOfMinorsCompute = 10*minNumberToCutDown + 8*ceiling(log_1.3(possibleMinors));)
    else (
        numberOfMinorsCompute = opts.MaxMinors;);

    minTerm := sub(0, R1);
    mutM1 := mutableMatrix(M1);
    nonzeroM := replaceZeros(M1); --
    mutM2 := mutableMatrix(nonzeroM); --for smallest grevlex computations

    searchedSet := new MutableHashTable from {}; --used to store which determinants have already been computed

    sumMinors := ideal(sub(0, ambR));

    quotient1 := ambR/(Id+sumMinors);
    d := dim(quotient1);

    if (opts.Verbose or debugLevel > 0) then print concatenate("Rn: ring dimension =", toString(d), ", there are ", toString(possibleMinors), " possible minors, we will compute up to ", toString(numberOfMinorsCompute), " of them.");

    i := 0;
    k := 0;
    j := 2*minNumberToCutDown+3;
    R2 := reorderPolynomialRing(GRevLex, ambR);
    f := map(R2, ambR);
    myRandom := 0;
    local M2;
    local submatrixS1;
    if (opts.Verbose or debugLevel > 0) then print concatenate("Rn: About to enter loop");
    while ( (r-d <= n1) and (i < numberOfMinorsCompute) and (#searchedSet < possibleMinors)) do (
        while (i <= 1.3^j) and (i < numberOfMinorsCompute) do (
            submatrixS1 = internalChooseMinor(fullRank, ambR, nonzeroM, M1, Strategy=>opts.Strategy, MutableSmallest=>mutM2, MutableLargest=>mutM1);
            if  (not (submatrixS1 === null)) and (not (searchedSet#?(locationToSubmatrix(submatrixS1)))) then (
                searchedSet#(locationToSubmatrix(submatrixS1)) = true;
                sumMinors = sumMinors + ideal(getDetOfSubmatrix(M1, submatrixS1, DetStrategy=>opts.DetStrategy));
                k = k+1;
            );
            i = i+1;
        );
        if (opts.Verbose or debugLevel > 0) then print concatenate("Rn:  Loop step, about to compute dimension.  Submatrices considered: ", toString(i), ", and computed = ", toString(# keys searchedSet) );
        mutM2 = mutableMatrix(nonzeroM); --reset this matrix periodically
        mutM1 = mutableMatrix(M1);
        quotient1 = ambR/(Id+sumMinors);
        d = dim(quotient1);
        if (opts.Verbose or debugLevel > 0) then print concatenate("Rn:  partial singular locus dimension computed, = ", toString(d));
        j = j+1;
    );
    if (opts.Verbose or debugLevel > 0) then print concatenate("Rn:  Loop completed, submatrices considered = ", toString(i), ", and computed = ", toString(# keys searchedSet), ".  singular locus dimension appears to be = ", toString(d));
    if ( (r-d) > n1) then (return true);
    if  (r-d<n+1) and (#searchedSet <= possibleMinors) then (return null);
    if (#searchedSet == possibleMinors) and ((r-d)<n1+1) then (return false);
);



chooseGoodMinors = method(Options=>optChooseGoodMinors);

chooseGoodMinors(ZZ, ZZ, Matrix) := opts -> (howMany, minorSize, M1) -> (
    R1 := ring M1;
    ambR := ambient R1;
    Id := ideal R1;
    possibleMinors := binomial(numColumns M1, minorSize)*binomial(numRows M1, minorSize);
    M1 = sub(M1, ambR);
    mutM1 := mutableMatrix(M1);
    nonzeroM := replaceZeros(M1); --
    mutM2 := mutableMatrix(nonzeroM); --for smallest grevlex computations

    searchedSet := new MutableHashTable from {}; --used to store which determinants have already been computed
    i := 0; -- how many times through the loop we go
    k := 0; -- how many actual minors we found
    j := 3; -- controlling when to reset the grevlex matrix
    maxAttempts := 5*log_2(possibleMinors) + 10; --the absolute max number of attempts
    R2 := reorderPolynomialRing(GRevLex, ambR);
    f := map(R2, ambR);
    sumMinors := ideal(sub(0, ambR));
    local M2;
    local submatrixS1;
--first we try smallest submatrices with respect to several different monomial orders
    while ( (k < howMany) and (i < maxAttempts) and (#searchedSet < possibleMinors)) do (
        while (i <= j^1.5) do (
            submatrixS1 = internalChooseMinor(minorSize, ambR, nonzeroM, M1, Strategy=>opts.Strategy, MutableSmallest=>mutM2, MutableLargest=>mutM1);

            if (not (submatrixS1 === null)) and (not searchedSet#?(locationToSubmatrix(submatrixS1))) then (
                searchedSet#(locationToSubmatrix(submatrixS1)) = true;
                sumMinors = sumMinors + ideal(getDetOfSubmatrix(M1, submatrixS1, DetStrategy=>opts.DetStrategy));
                k = k+1;
            );
            i = i+1;
        );
        mutM2 = mutableMatrix(nonzeroM); --reset this matrix periodically
        mutM1 = mutableMatrix(M1);
        j = j+1;
      --reset the matrix before randomly traversing it again
    );
    sub(sumMinors, R1)
)



projDim = method(Options=>optProjDim);

projDim(Module) := opts -> (N1) -> (
    if (isHomogeneous N1) then return pdim N1;
    ambRing := ring (N1);
    if (not instance(ambRing, PolynomialRing)) then error "projDim: currently this only works for modules over polynomial rings";
    myRes := resolution minimalPresentation N1;
    myDiffs := myRes.dd;
    myLength := length myRes;
    firstRank := rank myRes_myLength;
    goodMinorsOptions := new OptionTable from {Strategy=>opts.Strategy, Verbose=>opts.Verbose, DetStrategy=>opts.DetStrategy}; --just grab the options relevant to chooseGoodMinors
    if (debugLevel > 0) or opts.Verbose then print concatenate("projDim: resolution computed!  length =", toString myLength, " rank =", toString firstRank);
    firstDiff := myDiffs_myLength;

    local minorsCount; --how many minors?
    if instance(opts.MaxMinors, BasicList) then (
        minorsCount = (opts.MaxMinors)#0;)
    else if instance(opts.MaxMinors, ZZ) then (
        minorsCount = opts.MaxMinors;)
    else ( minorsCount = firstRank*4+5; );

    theseMinors := chooseGoodMinors(minorsCount, firstRank, firstDiff, goodMinorsOptions);
    if (debugLevel > 0) or opts.Verbose then print concatenate("projDim: first minors computed!  minors found =", toString(#first entries gens theseMinors));
    curDim := dim theseMinors;
    if (curDim >= 0) then (return myLength);

    i := myLength-1;
    if (debugLevel > 0) or opts.Verbose then print "projDim: computed dim, now starting loop.";
    while (i>opts.MinDimension) and (curDim <= -1) do (
        --print concatenate("in loop: ", toString(i));
        firstRank = (rank myRes_i)-firstRank;
        if instance(opts.MaxMinors, BasicList) then (
            minorsCount = (opts.MaxMinors)#(myLength - i); )
        else if instance(opts.MaxMinors, ZZ) then (
            minorsCount = opts.MaxMinors; )
        else ( minorsCount = firstRank*4+5; );
        theseMinors = chooseGoodMinors(minorsCount, firstRank, myDiffs_i, goodMinorsOptions);
        if (debugLevel > 0) or opts.Verbose then print concatenate("projDim: in loop, about to compute dim.  current length =", toString myLength, " rank =", toString firstRank);
        curDim = dim theseMinors;
        if (curDim >= 0) then (return i);
        i=i-1;
    );
    return i;
);

--checks if the rank is at least something.  First it calls getSubmatrixOfRank, if it finds one, then great, otherwise it calls rank
isRankAtLeast = method(Options => optIsRankAtLeast);

isRankAtLeast(ZZ, Matrix) := opts -> (n1, M0) -> (
    if (n1 > numRows M0) or (n1 > numColumns M0) then return false;
    if (n1 == numRows M0) and (n1 == numColumns M0) then return (rank M0 == n1);
    val := getSubmatrixOfRank(n1, M0, opts);
    if (val === null) then ( return (rank M0 >= n1); ) else return true;
);

getSubmatrixOfRank = method(Options => optIsRankAtLeast);

getSubmatrixOfRank(ZZ, Matrix) := opts -> (n1, M0) -> (
    --print opts;
    if (n1 > numRows M0) or (n1 > numColumns M0) then return null;
    R1 := ring M0;
    ambRing := ambient R1;
    i := 0;
    M1 := sub(M0, ambRing);
    possibleMinors := binomial(numColumns M1, n1)*binomial(numRows M1, n1);
    attempts := min(possibleMinors, 2+log_10(possibleMinors));
    if not (opts.MaxMinors === null) then attempts = opts.MaxMinors;
    mutM1 := mutableMatrix(M1); --for largest grevlex computations
    nonzeroM := replaceZeros(M1); --
    mutM2 := mutableMatrix(nonzeroM); --for smallest grevlex computations

    internalMinorsOptions := new OptionTable from {Strategy=>opts.Strategy}; --just grab the options relevant to chooseGoodMinors

    searchedSet := new MutableHashTable from {}; --used to store which ranks have already been computed

    if not (opts.MaxMinors === null) then (
        attempts = opts.MaxMinors;
    );
    subMatrix := null;
    val := null;
    if (debugLevel > 0) or opts.Verbose then print ("getSubmatrixOfRank: Trying to find a submatrix of rank at least: " | toString(n1) | " with attempts = " | toString(attempts) | ".  DetStrategy=>" | toString(opts.DetStrategy));
    while (i < attempts)  do (
        subMatrix = internalChooseMinor(n1, ambRing, nonzeroM, M1, internalMinorsOptions++{MutableSmallest=>mutM2, MutableLargest=>mutM1});
        --if (debugLevel > 0) or opts.Verbose then print ("getSubmatrixOfRank: found subMatrix " | toString(subMatrix));
        if (not (subMatrix === null)) and (not (searchedSet#?(locationToSubmatrix(subMatrix)))) then (
            searchedSet#(locationToSubmatrix(subMatrix)) = true;
            if (opts.DetStrategy === Rank) then (
                mRowListS := subMatrix#0;
                mColListS := subMatrix#1;
                if (rank sub((M1^mRowListS)_mColListS, R1) >= n1) then (
                    if (debugLevel > 0) or opts.Verbose then print ("getSubmatrixOfRank: found one, in " | toString(i+1) | " attempts");
                    return {mRowListS, mColListS};
                );
            )
            else if (opts.DetStrategy === null) then ( --we use a default strategy
                val = sub(getDetOfSubmatrix(M1, subMatrix, DetStrategy=>opts.DetStrategy), R1);
                if (not (val == 0)) then (if (debugLevel > 0) or opts.Verbose then print ("getSubmatrixOfRank: found one, in " | toString(i+1) | " attempts"); return {subMatrix#0, subMatrix#1};);
            )
            else ( --we use the given strategy
                val = sub(getDetOfSubmatrix(M1, subMatrix, DetStrategy=>opts.DetStrategy), R1);
                if (not (val == 0)) then (if (debugLevel > 0) or opts.Verbose then print ("getSubmatrixOfRank: found one, in " | toString(i+1) | " attempts"); return {subMatrix#0, subMatrix#1};);
            );
        );
        i = i+1;
    );
    return null;
);

--**********************************************
--**Alternate minors command, multithreaded  ***
--**computes smaller minors before larger    ***
--**minors.  It stores them.  This is faster ***
--**for matrices where multiplications are   ***
--**expensive.  It uses more ram though.     ***
--**********************************************

recursiveMinors =method(Options=>{MinorsCache => true, Threads=>0, Verbose=>false});

recursiveMinors(ZZ, Matrix) := Ideal => opts -> (n, M1) -> (
    if not ((M1#cache)#?MinorsCache) then (M1.cache)#MinorsCache = new MutableHashTable from {};
    --if not ((M1#cache)#?MinorsCache) then print "hello?";
    --if (n == numColumns M1) then M1 = transpose M1;
    H1 := null;
    if ((M1#cache)#MinorsCache)#?n then (
        if (debugLevel > 0) or (opts.Verbose) then print "recursiveMinors: we found stored data, using it.";
        H1 = ((M1#cache)#MinorsCache)#n; )
    else (
        H1 = recursiveMinorsTableP(n, M1, opts);
        if (opts.MinorsCache == true) then (((M1.cache).MinorsCache)#n = H1;); --store data only if told to
    );
    ideal values H1
);


--this is the multithreaded internal function
recursiveMinorsTableP = method(Options=>{MinorsCache => true, Threads=>0, Verbose=>true});

recursiveMinorsTableP(ZZ, Matrix) := opts -> (n, M1) -> (
    if not ((M1#cache)#?MinorsCache) then (M1.cache)#MinorsCache = new MutableHashTable from {};
    if (opts.Threads <= 1) then (
        if (debugLevel > 0) then print "recursiveMinorsP: Going to single threaded version.";
        return recursiveMinorsTable(n, M1);
    );
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    lenL := length L;
    R1 := {};
    if (debugLevel > 0) then    print "recursiveMinorsTableP: about to recurse";
    if (n > 2) then (
        if ((M1.cache).MinorsCache)#?(n-1) then (
            if (debugLevel > 0) then print "recursiveMinorsTable: we found stored data, using it to avoid recursion.";
            R1 = ((M1.cache).MinorsCache)#(n-1); )
        else (
            R1 = recursiveMinorsTableP(n-1, submatrix'(M1, {rowCt-1}, ),opts);
        );
    )
    else (
        if ((M1.cache).MinorsCache)#?n then (
            if (debugLevel > 0) then print "recursiveMinorsTableP: we found stored data, using it.";
            return ((M1.cache).MinorsCache)#n; )
        else (
            return new HashTable from  apply(L, i -> (i => det (M1^(i#0)_(i#1))) ); --this should be multithreaded too presumably...
        );
    );
    tempL := null;
    if (debugLevel > 0) then print "3TableP: did recursion, going to do tasks";
    taskList := apply(opts.Threads, i -> (tempL =  take(L, {(i*(lenL-1))//(opts.Threads), ((i+1)*(lenL-1))//(opts.Threads)});
        return createTask(temp, (tempL, n, M1, R1) ); ) );
    apply(taskList, t -> schedule t);
    while true do (
        sleep 0;--this should be replaced by a usleep or nanosleep command.
        if (all(taskList, t->isReady(t))) then break;
        );
    myList := flatten flatten apply(taskList, t -> taskResult(t));
    H := new HashTable from myList;
    apply(taskList, tt -> cancelTask(tt));
    return H;
)

--this is just a temporary function
temp = (Lis,n, M1, R1) -> (
    --print "temp started";
	K := apply(Lis, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R1#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol})))));
	return K
);

--this is the single threaded function
recursiveMinorsTable = (n, M1) -> (
    if not ((M1#cache)#?MinorsCache) then (M1.cache)#MinorsCache = new MutableHashTable from {};
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    R := {};
    if (n > 2) then (
        if ((M1.cache).MinorsCache)#?(n-1) then (
            if (debugLevel > 0) then print "recursiveMinorsTable: we found stored data, using it to avoid recursion.";
            R = ((M1.cache).MinorsCache)#(n-1); )
        else (
            R = recursiveMinorsTable(n-1, submatrix'(M1, {rowCt-1}, )))
        )
    else (
        if ((M1.cache).MinorsCache)#?n then (
            if (debugLevel > 0) then print "recursiveMinorsTable: we found stored data, using it.";
            return ((M1.cache).MinorsCache)#n; )
        else (
            return new HashTable from  apply(L, i -> (i => det (M1^(i#0)_(i#1))) )
        );
    );
    H := new HashTable from  apply(L, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol}))))); --det (M^(i#0)_(i#1))

    return H;
)



beginDocumentation();

document {
    Key => FastLinAlg,
    Headline => "faster linear algebra in certain cases",
    EM "FastLinAlg", " is a for computing (applications of) function field linear algebra more quickly in certain settings.",
    BR{},BR{},
    "This package does... ", TO "chooseGoodMinors", " which ",
    BR{},BR{},
    BOLD "Useful functions:",BR{},
    UL {
	  {TO "chooseGoodMinors", " can let you determine if a divisor is Cartier or if a power is Cartier."},
	  {TO "Rn", " checks whether a ring is regular in codimension n in the affirmative." },
	  {TO "projDim", " checks the projective dimension of a module and will give better answers than " },
	},
	BR{},
	BOLD "Acknowledgements:",BR{},BR{},
	"The authors would like to thank David Eisenbud and Srikanth Iyengar for useful conversations and comments on the development of this package."
}

doc ///
    Key
        chooseGoodMinors
        (chooseGoodMinors, ZZ, ZZ, Matrix)
        [chooseGoodMinors, Verbose]
    Headline
        returns an ideal generated by interesting minors in a matrix
    Usage
        chooseGoodMinors(count, minorSize, M)
    Inputs
        count: ZZ
        minorSize: ZZ
        M: Matrix
    Outputs
        : Ideal
    Description
        Text
            This returns an ideal generated by approximately {\tt count} minors of size {\tt minorSize} the matrix M.
        Example
            R = QQ[x, y, z];
///

doc ///
    Key
        chooseSubmatrixSmallestDegree
        (chooseSubmatrixSmallestDegree, ZZ, Matrix)
    Headline
        returns coordinates for low degree submatrix of a matrix
    Usage
        chooseSubmatrixSmallestDegree(n1, M3)
    Inputs
        n1: ZZ
        M3: Matrix
    Outputs
        : List
    Description
        Text
            returns the row and column coordinates in the original matrix for a specified sized submatrix of terms lower in degree than other terms of the matrix.
        Example
            R = QQ[x,y,z];
            M=matrix{{x^2,x^3,x^4},{y^4,y^2,y^3},{z^3,z^4,z^2}};
            chooseSubmatrixSmallestDegree(2, M)
///

doc ///
    Key
        chooseSubmatrixLargestDegree
        (chooseSubmatrixLargestDegree, ZZ, Matrix)
    Headline
        returns coordinates for higher degree submatrix of a matrix
    Usage
        chooseSubmatrixLargestDegree(n1, M1)
    Inputs
        n1: ZZ
        M1: Matrix
    Outputs
        : List
    Description
        Text
            returns the row and column coordinates in the original matrix for a specified sized submatrix of
            terms likely higher in degree than other terms of the matrix.
        Example
            R = QQ[x,y,z];
            M=matrix{{x^2,x^3,x^4},{y^4,y^2,y^3},{z^3,z^4,z^2}};
            chooseSubmatrixLargestDegree(2, M)
///

doc ///
    Key
        chooseRandomSubmatrix
        (chooseRandomSubmatrix, ZZ, Matrix)
    Headline
        returns coordinates for a random submatrix
    Usage
        chooseRandomSubmatrix(n1, M1)
    Inputs
        n1: ZZ
        M1: Matrix
    Outputs
        : List
    Description
        Text
            Returns the row and column coordinates in the original matrix for a specified sized random submatrix.
        Example
            R = QQ[x,y,z];
            M=matrix{{x^2,x^3,x^4},{y^4,y^2,y^3},{z^3,z^4,z^2}};
            chooseRandomSubmatrix(2, M)
///

doc ///
    Key
        isRankAtLeast
        (isRankAtLeast, ZZ, Matrix)
        [isRankAtLeast, Verbose]
    Headline
        determines if the matrix has rank at least a number
    Usage
        isRankAtLeast(n1, M1)
    Inputs
        n1: ZZ
        M1: Matrix
    Outputs
        : List
    Description
        Text
            This function tries to quickly determine if the matrix has a given rank.
            {\tt isRankAtLeast} calls @TO getSubmatrixOfRank@.  If that function finds
            a submatrix of a certain rank, this returns true.  If that function
            fails to find a submatrix of a certain rank, this simply calls
            @TO rank@.  To control the number of times {\tt getSubmatrixOfRank} considers
            submatrices, use the option {\tt MaxMinors}.
        Example
            R = QQ[x,y,z]
    SeeAlso
        getSubmatrixOfRank
        [isRankAtLeast, Strategy]
///

doc ///
    Key
        getSubmatrixOfRank
        (getSubmatrixOfRank, ZZ, Matrix)
        [getSubmatrixOfRank, Verbose]
    Headline
        tries to find a submatrix of the given rank
    Usage
        getSubmatrixOfRank(n1, M1)
    Inputs
        n1: ZZ
        M1: Matrix
    Outputs
        : List
    Description
        Text
            This function looks at subtrices of the given matrix, and tries to find
            one of the specified rank.  If it fails to find such a matrix,
            the function returns null.  The option {\tt MaxMinors} is used to
            control how many minors to consider.  If left {\tt null}, the number
            considered is based on the size of the matrix.
        Example
            R = QQ[x,y,z]
        Text
            The option {\tt Strategy} is used to used to control how the function computes
            the rank of the submatrices considered.  See @TO [getSubmatrixOfRank, Strategy]@.
    SeeAlso
        isRankAtLeast
        [getSubmatrixOfRank, Strategy]
///

doc ///
    Key
        Rn
        (Rn, ZZ, Ring)
        [Rn, Verbose]
    Headline
        attempts to show that the ring is regular in codimension n
    Usage
        Rn(n, R)
    Inputs
        n: ZZ
        R: Ring
    Outputs
        : Boolean
            true, if the ring is Rn, false if it determines it is not
        : null
            if no determination is made
    Description
        Text
            This function looks at interesting minors of the jacobian matrix to try to verify that the ring is Rn.
            It is frequently much faster at giving an affirmative answer than computing the dimension of the ideal of all minors of the Jacobian.
            We begin with a simple example which is R1, but not R2.
        Example
            R = QQ[x, y, z]/ideal(x*y-z^2);
            Rn(1, R)
            Rn(2, R)
        Text
            Next we consider a more interesting example that is R1 but not R2, and highlight the speed differences.  Note that {\tt Rn(2, R)} returns nothing, as the function did not determine if the ring was Rn.
        Example
            T = ZZ/101[x1,x2,x3,x4,x5,x6,x7];
            I =  ideal(x5*x6-x4*x7,x1*x6-x2*x7,x5^2-x1*x7,x4*x5-x2*x7,x4^2-x2*x6,x1*x4-x2*x5,x2*x3^3*x5+3*x2*x3^2*x7+8*x2^2*x5+3*x3*x4*x7-8*x4*x7+x6*x7,x1*x3^3*x5+3*x1*x3^2*x7+8*x1*x2*x5+3*x3*x5*x7-8*x5*x7+x7^2,x2*x3^3*x4+3*x2*x3^2*x6+8*x2^2*x4+3*x3*x4*x6-8*x4*x6+x6^2,x2^2*x3^3+3*x2*x3^2*x4+8*x2^3+3*x2*x3*x6-8*x2*x6+x4*x6,x1*x2*x3^3+3*x2*x3^2*x5+8*x1*x2^2+3*x2*x3*x7-8*x2*x7+x4*x7,x1^2*x3^3+3*x1*x3^2*x5+8*x1^2*x2+3*x1*x3*x7-8*x1*x7+x5*x7);
            S = T/I;
            dim S
            time Rn(1, S)
            time Rn(2, S)
        Text
            On one machine, a call of {\tt dim singularLocus R} took {\tt 32.4468} seconds and returned {\tt 1}.  On the same machine, {\tt Rn(1, R)} took only {\tt 0.153042} seconds.
        Text
            The following is a (pruned) affine chart on Abelian surface (the product of elliptic curves).  It is nonsingular, as our function verifies.
            If one does not prune it, then the {\tt dim singularLocus} call takes an enormous amount of time, otherwise {\tt dim singularLocus} frequently takes a similar amount of time to our function.
        Example
            R = QQ[c, f, g, h]/ideal(g^3+h^3+1,f*g^3+f*h^3+f,c*g^3+c*h^3+c,f^2*g^3+f^2*h^3+f^2,c*f*g^3+c*f*h^3+c*f,c^2*g^3+c^2*h^3+c^2,f^3*g^3+f^3*h^3+f^3,c*f^2*g^3+c*f^2*h^3+c*f^2,c^2*f*g^3+c^2*f*h^3+c^2*f,c^3-f^2-c,c^3*h-f^2*h-c*h,c^3*g-f^2*g-c*g,c^3*h^2-f^2*h^2-c*h^2,c^3*g*h-f^2*g*h-c*g*h,c^3*g^2-f^2*g^2-c*g^2,c^3*h^3-f^2*h^3-c*h^3,c^3*g*h^2-f^2*g*h^2-c*g*h^2,c^3*g^2*h-f^2*g^2*h-c*g^2*h,c^3*g^3+f^2*h^3+c*h^3+f^2+c);
            dim(R)
            time (dim singularLocus (R))
            time Rn(2, R)
            time Rn(2, R)
            time Rn(2, R)
        Text
            The function works by choosing interesting looking submatrices, computing their determinants, and periodically (based on a logarithmic growth setting), computing the dimension of a subideal of the Jacobian.
            The option {\tt Verbose} can be used to see this in action.
        Example
            time Rn(2, S, Verbose=>true)
        Text
            The maximum number of minors considered can be controlled by the option {\tt MaxMinors}.
            If the user does not specify this, then the default value is given by the formula {\tt TODO}.
            The function does not recompute determinants, so {\tt MaxMinors} is only an upper bound on the number of
            minors computed.
        Example
            time Rn(2, S, Verbose=>true, MaxMinors=>30)
        Text
            This function has many options which allow you to fine tune the strategy used to find interesting minors.
            You can pass it a {\tt HashTable} specifying the strategy via the option {\tt Strategy}.  See @TO LexSmallest@ for how to construct this {\tt HashTable}.
            The default strategy is {\tt StrategyDefault}, which seems to work well on the examples we have explored.  However, caution must be taken, even in the above examples,
            certain strategies work well while others do not.  In the Abelian surface example, {\tt LexSmallest} works very well,
            while {\tt LexSmallestTerm} does not even typically correctly identify the ring as nonsingular
            (this is because of the fact that there are a small number of entries with nonzero constant terms, which are selected repeatedly).
            However, in our first example, the {\tt LexSmallestTerm} is much faster, and {\tt Random} does not perform well at all.
        Example
            StrategyCurrent#Random = 0;
            StrategyCurrent#LexSmallest = 100;
            StrategyCurrent#LexSmallestTerm = 0;
            time Rn(2, R, Strategy=>StrategyCurrent)
            time Rn(2, R, Strategy=>StrategyCurrent)
            time Rn(1, S, Strategy=>StrategyCurrent)
            time Rn(1, S, Strategy=>StrategyCurrent)
            StrategyCurrent#LexSmallest = 0;
            StrategyCurrent#LexSmallestTerm = 100;
            time Rn(2, R, Strategy=>StrategyCurrent)
            time Rn(2, R, Strategy=>StrategyCurrent)
            time Rn(1, S, Strategy=>StrategyCurrent)
            time Rn(1, S, Strategy=>StrategyCurrent)
            time Rn(1, S, Strategy=>StrategyRandom)
            time Rn(1, S, Strategy=>StrategyRandom)
///



document {
    Key => {GRevLexLargest, GRevLexSmallest, GRevLexSmallestTerm, LexLargest, LexSmallest, LexSmallestTerm, Random, RandomNonzero, [chooseGoodMinors, Strategy], [getSubmatrixOfRank, Strategy], [Rn, Strategy], [isRankAtLeast, Strategy], "StrategyCurrent", "StrategyDefault", "StrategyLexSmallest", "StrategyGRevLexSmallest", "StrategyRandom"},
    Headline => "strategies for choosing submatrices",
    "Many of the core functions of this package allow the user to fine tune the strategy of how submatrices are selected.  Different strategies yield markedly different performance or results on these examples.
    These are controlled by specifying a {\tt Strategy => } option, pointing to a {\tt HashTable}.  This HashTable should have the following keys.",
    UL {
        {TT "GRevLexLargest", ": try to find submatrices where each row and column has a large entry with respect to a random ", TT "GRevLex", "order."},
        {TT "GRevLexSmallest", ": try to find submatrices where each row and column has a small entry with respect to a random ", TT "GRevLex", "order."},
        {TT "GRevLexSmallestTerm", ": find submatrices where each row and column has an entry with a small term with respect to a random ", TT "GRevLex", "order."},
        {TT "LexLargest", ": try to find submatrices where each row and column has a large entry with respect to a random ", TT "Lex", "order."},
        {TT "LexSmallest", ": try to find submatrices where each row and column has a small entry with respect to a random ", TT "Lex", "order."},
        {TT "LexSmallestTerm", ": find submatrices where each row and column has an entry with a small term with respect to a random ", TT "Lex", "order."},
        {TT "Random", ": find random submatrices "},
        {TT "RandomNonzero", ": find random submatrices that have nonzero rows and columns"}
    },
    "Each such key should point to an integer.",
	BR{},BR{},
    "Functions such as ", TO chooseGoodMinors, " will select a number of random submatrices based on the values of those keys.  For example, if ",
    TT "LexSmallest", " and ", TT "LexLargest", " are set to ", TT "50", " approximately the submatrics will be smallest with respect to ", TT "Lex",
    " and the other half will be largest with respect to ", TT "Lex.", "The values do not need to add up to 100.",
    BR{},BR{},
    "These functions all work by finding the optimal entry with respect to the given strategy, removing that row and column, and then choosing the next optimal entry.  ",
    "This is done until a submatrix of the desired size has been found.",
    BR{},BR{},
    "In some functions, the ", TT "GRevLex", " versions of this strategy will modify the working matrix in a loop, repeatedly lowering/raising the degree of elements",
    "so as to ensure that different choices are made.",
    BR{},BR{},
    "This package comes with several default strategies exported to the user.",
    UL {
        {TT "StrategyDefault", ": 40% of the matrices are ", TT "LexSmallest", " and ", TT "LexSmallestTerm", " each, the remaining 20% are ", TT "Random" },
        {TT "StrategyLexSmallest", ": 50% of the matrices are ", TT "LexSmallest", " and 50% are ", TT "LexSmallestTerm"},
        {TT "StrategyGRevLexSmallest", ": 50% of the matrices are ", TT "GRevLexSmallest", " and 50% are ", TT "GRevLexLargest"},
        {TT "StrategyRandom", ": chooses 100% random submatrices."},
    },
    "Additionally, a ", TT "MutableHashTable", " named ", TT "StrategyCurrent", " is also exported.  It begins as the default strategy, but the user can modify it."
}


doc ///
    Key
        recursiveMinors
        (recursiveMinors, ZZ, Matrix)
        [recursiveMinors, MinorsCache]
        [recursiveMinors, Threads]
        MinorsCache
    Headline
        uses a recursive cofactor algorithm to compute the ideal of minors of a matrix
    Usage
        I = recursiveMinors(n, M, Threads=>t, MinorsCache=>b)
    Inputs
        n: ZZ
           the size of minors to compute
        M: Matrix
           the matrix which to compute the ideal of minors
        t: ZZ
           an optional input, which describes the number of threads to uses
        b: Boolean
           an optional input, which says whether to cache in input
    Outputs
        I: Ideal
           the ideal of minors
    Description
        Text
           Given a matrix $M$, this computes the ideal of determinants of size $n \times n$ submatrices.
           This is a recursive function.  It computes the minors via cofactor.  This function is {\it currently}
           implemented in Macaulay2, whereas the {\tt minors} function uses
        Example
           R = QQ[x,y];
           M = random(R^{5,5,5,5,5,5}, R^7);
           time I2 = recursiveMinors(4, M, Threads=>0);
           time I1 = minors(4, M, Strategy=>Cofactor);
           I1 == I2
///

doc ///
    Key
        projDim
        (projDim, Module)
        [projDim, Verbose]
        [projDim, MaxMinors]
    Headline
        finds an upper bound for the projective dimension of a module
    Usage
        n = projDim(N, MinDimension=>d)
    Inputs
        N: Module
            a module over a polynomial ring
        d: ZZ
            the minimum projective dimension of the module
    Outputs
        n: ZZ
           an upper bound for the projective dimension of N
    Description
        Text
            The function {\tt pdim} returns the length of a projective resolution.
            If the module passed is not Homogeneous, then the projective resolution may not be minimal and so {\tt pdim} can
            give the wrong answer.  This function {\tt projDim} tries to improve this bound
            by considering ideals of appropriately sized minors of the resolution (starting from the end of the resolution and working backwards).
            Using the option {\tt MinDimension} (default value 0)
            gives a lower bound on the projective dimension, increasing it can thus improve the speed of computation.
        Example
            R = QQ[x,y];
            I = ideal((x^3+y)^2, (x^2+y^2)^2, (x+y^3)^2, (x*y)^2);
            pdim(module I)
            time projDim(module I)
            time projDim(module I, MinDimension => 1)
        Text
            The option MaxMinors can be used to control how many minors are computed at each step.
            If this is not specified, the number of minors is a function of the rank of the matrix at the
            current step of the resolution.
            Otherwise the user can set the option {\tt MaxMinors => ZZ} to specify that a fixed integer is used for
            each step.  Alternately, the user can set {\tt MaxMinors => List} passing a list of how many minors
            should be computed at each step (working backwards).  TODO, give an example. 
    SeeAlso
        pdim
///

doc ///
    Key
        MinDimension
    Headline
        an option for projDim
    Description
        Text
            This option is used to tell the function {\tt projDim} not to look for projective dimension below the option value.
    SeeAlso
        projDim
///

doc ///
    Key
        reorderPolynomialRing
        (reorderPolynomialRing, Symbol, PolynomialRing)
    Headline
        produces an isomorphic polynomial ring with a different, randomized, monomial order
    Usage
        R1 = reorderPolynomialRing(orderType, R)
    Inputs
        R: Ring
            a polynomial ring
        orderType: Symbol
            a valid monomial order, such as {\tt GRevLex}
    Outputs
        R1: ZZ
           a polynomial ring
    Description
        Text
            This function takes a polynomial ring and produces a new polynomial ring with {\tt MonomialOrder} of type {\tt orderType}.
            The order of the variables are randomized.
        Example
            R = QQ[x,y,z,w];
            x > y and y > z and z > w
            use reorderPolynomialRing(GRevLex, R)
            x > y
            y > z
            z > w
///

doc ///
    Key
        DetStrategy
        Recursive
        Rank
        [chooseGoodMinors, DetStrategy]
        [getSubmatrixOfRank, DetStrategy]
        [isRankAtLeast, DetStrategy]
        [projDim, DetStrategy]
        [Rn, DetStrategy]
        [RnReductionP, DetStrategy]
    Headline
        DetStrateg is a strategy for allowing the user to choose how determinants (or rank), is computed
    Description
        Text
            Passing the option {\tt DetStrategy => Symbol} controls how certain functions compute determinants (or in some cases, rank).
            For all methods, {\tt Bareiss}, {\tt Cofactor} and {\tt Recursive} are valid options.  The first two come included with Macaulay2,
            the third is just a call to @TO recursiveMinors@.
        Text
            Additionally, for the methods {\tt getSubmatrixOfRank} and {\tt isRankAtLeast},
            one can also pass {\tt DetStrategy} the option {\tt Rank} (the default).  Not using {\tt Rank} tells
            the functions to check rank by computing the determinant, instead of using the internal
            {\tt rank} function.  This is not usually recommended, but soemtimes it can be effective.
    SeeAlso
        recursiveMinors
        Bareiss
        Cofactor
        rank
///

TEST /// --check #0 (Rn)
R = QQ[x,y,z]/ideal(x^2-y*z);
assert(Rn(1, R)=== true);
///

TEST /// --check #1 (Rn)
T = ZZ/101[x1,x2,x3,x4,x5,x6,x7];
I =  ideal(x5*x6-x4*x7,x1*x6-x2*x7,x5^2-x1*x7,x4*x5-x2*x7,x4^2-x2*x6,x1*x4-x2*x5,x2*x3^3*x5+3*x2*x3^2*x7+8*x2^2*x5+3*x3*x4*x7-8*x4*x7+x6*x7,x1*x3^3*x5+3*x1*x3^2*x7+8*x1*x2*x5+3*x3*x5*x7-8*x5*x7+x7^2,x2*x3^3*x4+3*x2*x3^2*x6+8*x2^2*x4+3*x3*x4*x6-8*x4*x6+x6^2,x2^2*x3^3+3*x2*x3^2*x4+8*x2^3+3*x2*x3*x6-8*x2*x6+x4*x6,x1*x2*x3^3+3*x2*x3^2*x5+8*x1*x2^2+3*x2*x3*x7-8*x2*x7+x4*x7,x1^2*x3^3+3*x1*x3^2*x5+8*x1^2*x2+3*x1*x3*x7-8*x1*x7+x5*x7);
assert((Rn(1,T/I) === true) or (Rn(1,T/I) === true));
///

TEST /// --check #2 (ProjDim)
R = QQ[x,y];
I = ideal((x^3+y)^2, (x^2+y^2)^2, (x+y^3)^2, (x*y)^2);
assert(projDim(module I)==1);
///

TEST /// --check #3 (Rn, RnReductionP)
R = (QQ)[YY_1, YY_2, YY_3, YY_4, YY_5, YY_6, YY_7, YY_8, YY_9];
J =  ideal(YY_8^2-YY_7*YY_9,YY_6*YY_8-YY_5*YY_9,YY_3*YY_8-YY_2*YY_9,YY_2*YY_8-YY_1*YY_9,YY_6*YY_7-YY_5*YY_8,YY_3*YY_7-YY_1*YY_9,YY_2*YY_7-YY_1*YY_8,YY_6^2-YY_4*YY_9,YY_5*YY_6-YY_4*YY_8,YY_4*YY_6+YY_1*YY_8-10*YY_1*YY_9-YY_
     2*YY_9+10*YY_3*YY_9,YY_3*YY_6-YY_8*YY_9-10*YY_9^2,YY_2*YY_6-YY_7*YY_9-10*YY_8*YY_9,YY_1*YY_6-YY_7*YY_8-10*YY_7*YY_9,YY_5^2-YY_4*YY_7,YY_4*YY_5+YY_1*YY_7-10*YY_1*YY_8-YY_1*YY_9+10*YY_2*YY_9,YY_3*YY_5-YY_7*YY_9-10*YY_8
     *YY_9,YY_2*YY_5-YY_7*YY_8-10*YY_7*YY_9,YY_1*YY_5-YY_7^2-10*YY_7*YY_8,YY_4^2+YY_7^2-YY_9^2,YY_3*YY_4-YY_5*YY_9-10*YY_6*YY_9,YY_2*YY_4-YY_5*YY_8-10*YY_5*YY_9,YY_1*YY_4-YY_5*YY_7-10*YY_5*YY_8,YY_2^2-YY_1*YY_3,YY_1*YY_2-
     10*YY_1*YY_3-YY_2*YY_3+10*YY_3^2+YY_4*YY_8+10*YY_4*YY_9,YY_1^2-YY_3^2+YY_4*YY_7+20*YY_4*YY_8-YY_4*YY_9);
assert((Rn(1, R/J)===null) and (RnReductionP(1, R/J, 7)===null));
///

TEST/// --check #4 (Rn)
T = (ZZ/101)[YY_1, YY_2, YY_3, YY_4, YY_5];
J = ideal(YY_2*YY_3+3*YY_3^2+43*YY_1*YY_4+YY_3*YY_4-43*YY_2*YY_5+50*YY_3*YY_5,YY_2^2-18*YY_3^2+YY_1*YY_4-18*YY_2*YY_4-39*YY_3*YY_4-YY_1*YY_5+8*YY_2*YY_5-47*YY_3*YY_5,YY_1^2+16*YY_1*YY_2-3*YY_1*YY_3-32*YY_3^2+23*YY_1*YY_4-42*YY_2*YY_4-43*YY_3*YY_4+19*YY_1*YY_5+20*YY_2*YY_5-34*YY_3*YY_5,YY_3^3+16*YY_1*YY_2*YY_4+8*YY_1*YY_3*YY_4-36*YY_3^2*YY_4-YY_1*YY_4^2-47*YY_3*YY_4^2+45*YY_1*YY_3*YY_5-31*YY_3^2*YY_5+7*YY_1*YY_4*YY_5+16*YY_2*YY_4*YY_5+37*YY_3*YY_4*YY_5-16*YY_1*YY_5^2+36*YY_2*YY_5^2+15*YY_3*YY_5^2,YY_1*YY_3^2-48*YY_1*YY_2*YY_4+21*YY_1*YY_3*YY_4-45*YY_3^2*YY_4+47*YY_1*YY_4^2+10*YY_2*YY_4^2-47*YY_3*YY_4^2+13*YY_4^3-25*YY_1*YY_2*YY_5-33*YY_1*YY_3*YY_5+45*YY_3^2*YY_5+24*YY_1*YY_4*YY_5-36*YY_2*YY_4*YY_5-41*YY_3*YY_4*YY_5+26*YY_4^2*YY_5-27*YY_1*YY_5^2+30*YY_2*YY_5^2-13*YY_3*YY_5^2-24*YY_4*YY_5^2-17*YY_5^3);
assert(Rn(1, T/J)===true);
--we shoudl change these a bit
--assert(Rn(1, T/J, LexSmallest=>0, Random=>0)===true);
-- note: if we set LexLargest=>0 running just on LexSmallest then Rn returns null
-- same happens for all =>0, GRevLexSmallest=> nonzero
-- GRevLexLargest works
-- assert(Rn(1, T/J, LexLargest=>0, LexSmallest=>0, GRevLexLargest=>20, Random=>0)===true);
-- assert(Rn(1, T/J, LexLargest=>0, LexSmallest=>0)===true);
-- assert(Rn(1, T/J, LexLargest=>0, LexSmallest=>0, Random=>0, RandomNonzero=>20)===true);
///

TEST/// --check #5, this is an affine chart on an Abelian surface, it is nonsingular
S = QQ[c, f, g, h];
J = ideal(g^3+h^3+1,f*g^3+f*h^3+f,c*g^3+c*h^3+c,f^2*g^3+f^2*h^3+f^2,c*f*g^3+c*f*h^3+c*f,c^2*g^3+c^2*h^3+c^2,f^3*g^3+f^3*h^3+f^3,c*f^2*g^3+c*f^2*h^3+c*f^2,c^2*f*g^3+c^2*f*h^3+c^2*f,c^3-f^2-c,c^3*h-f^2*h-c*h,c^3*g-f^2*g-c*g,c^3*h^2-f^2*h^2-c*h^2,c^3*g*h-f^2*g*h-c*g*h,c^3*g^2-f^2*g^2-c*g^2,c^3*h^3-f^2*h^3-c*h^3,c^3*g*h^2-f^2*g*h^2-c*g*h^2,c^3*g^2*h-f^2*g^2*h-c*g^2*h,c^3*g^3+f^2*h^3+c*h^3+f^2+c);
assert((Rn(2, S/J) === true) or (Rn(2, S/J) === true));
///

TEST /// --check #6, we found this example by dehomogenizing a homongeneous example, pdim does not provide the correct answer (of course, it does if you rehomogenize)
S = QQ[t_0, t_1, t_2, t_3, t_4, t_5];
J = ideal(-t_2^3+2*t_1*t_2*t_3-t_0*t_3^2-t_1^2*t_4+t_0*t_2*t_4,-t_2^2*t_3+t_1*t_3^2+t_1*t_2*t_4-t_0*t_3*t_4-t_1^2*t_5+t_0*t_2*t_5,-t_2*t_3^2+t_2^2*t_4+t_1*t_3*t_4-t_0*t_4^2-t_1*t_2*t_5+t_0*t_3*t_5,-t_3^3+2*t_2*t_3*t_4-t_1*t_4^2-t_2^2*t_5+t_1*t_3*t_5,-t_2^2*t_4+t_1*t_3*t_4+t_1*t_2*t_5-t_0*t_3*t_5-t_1^2+t_0*t_2,-t_2*t_3*t_4+t_1*t_4^2+t_2^2*t_5-t_0*t_4*t_5-t_1*t_2+t_0*t_3,-t_3^2*t_4+t_2*t_4^2+t_2*t_3*t_5-t_1*t_4*t_5-t_2^2+t_1*t_3,-t_2*t_4^2+t_2*t_3*t_5+t_1*t_4*t_5-t_0*t_5^2-t_1*t_3+t_0*t_4,
    -t_3*t_4^2+t_3^2*t_5+t_2*t_4*t_5-t_1*t_5^2-t_2*t_3+t_1*t_4,-t_4^3+2*t_3*t_4*t_5-t_2*t_5^2-t_3^2+t_2*t_4);
assert(projDim(module J, MinDimension=>2) == 2)
///




end
