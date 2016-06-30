
+ M2 --no-readline --print-width 79
Macaulay2, version 1.9.0.1
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : load "EthRoots"
stdio:1:1:(3): error: file not found on path: "EthRoots"

i2 : load "EthRoots.m2"

i3 : R = QQ[x,y,z]

o3 = R

o3 : PolynomialRing

i4 : I = ideal(x^3, x*y*z, x*y^2, z^3)

             3            2   3
o4 = ideal (x , x*y*z, x*y , z )

o4 : Ideal of R

i5 : kk = ZZ/7

o5 = kk

o5 : QuotientRing

i6 : R = kk[x,y,z]

o6 = R

o6 : PolynomialRing

i7 : I = ideal(x^3, x*y*z, x*y^2, z^3)

             3            2   3
o7 = ideal (x , x*y*z, x*y , z )

o7 : Ideal of R

i8 : ethRoot(4, I)

o8 = ideal 1

o8 : Ideal of R

i9 : ethRoot(1, I)

o9 = ideal 1

o9 : Ideal of R

i10 : R= kk[x,y,z]/(x^5 - y^3, x^7 - z^3, y^7 - z^5)

o10 = R

o10 : QuotientRing

i11 : I = (x^2 + 1)

       2
o11 = x  + 1

o11 : R

i12 : ethRoot(2, I)
stdio:12:1:(3): error: no method found for applying ethRoot to:
     argument 1 :  2 (of class ZZ)
                    2
     argument 2 :  x  + 1 (of class R)

i13 : I = ideal(x^2 + 1)

             2
o13 = ideal(x  + 1)

o13 : Ideal of R

i14 : ethRoot(2, I)
EthRoots.m2:61:43:(3):[3]: error: ethRoot: Expected an ideal in a PolynomialRing.
EthRoots.m2:61:43:(3):[3]: --entering debugger (type help to see debugger commands)
EthRoots.m2:61:43-61:43: --source code:
    if (class R =!= PolynomialRing) then (error "ethRoot: Expected an ideal in a PolynomialRing.");

ii15 : break

i16 : R = kk[x,y,z]

o16 = R

o16 : PolynomialRing

i17 : I = ideal(x^21, x^7*y^7*z^7, y^5) 

              21   7 7 7   5
o17 = ideal (x  , x y z , y )

o17 : Ideal of R

i18 : ethRoot(2, I)

o18 = ideal 1

o18 : Ideal of R

i19 : ethRoot(1, I)

o19 = ideal 1

o19 : Ideal of R

i20 : I = ideal(x^8, y^8, z^8)

              8   8   8
o20 = ideal (x , y , z )

o20 : Ideal of R

i21 : ethRoot(1, I)

o21 = ideal (z, y, x)

o21 : Ideal of R

i22 : benchmark "ethRoot(1, I)"

o22 = .00557421034782609

o22 : RR (of precision 53)

i23 : benchmark "ethRoot(1, I, EthRootStrategy=>MonomialBasis)"

o23 = .00130697908047337

o23 : RR (of precision 53)

i24 : ethRoot(1, ideal(x^7))

o24 = ideal x

o24 : Ideal of R

i25 : ethRoot(1, ideal(x^8))

o25 = ideal x

o25 : Ideal of R

i26 : ethRoot(1, ideal(x^13))

o26 = ideal x

o26 : Ideal of R

i27 : ethRoot(1, ideal(x^14))

             2
o27 = ideal x

o27 : Ideal of R

i28 : I = ideal(x^2)

             2
o28 = ideal x

o28 : Ideal of R

i29 : I_*

        2
o29 = {x }

o29 : List

i30 : I = ideal(x^2, y^2, x*y, x+y, x-y)

              2   2
o30 = ideal (x , y , x*y, x + y, x - y)

o30 : Ideal of R

i31 : I_*

        2   2
o31 = {x , y , x*y, x + y, x - y}

o31 : List

i32 : I = ideal(x^2, y^2, y^2*x, x^2*y, x+y, x^2 + y, 2*x + 2*x)

              2   2     2   2           2
o32 = ideal (x , y , x*y , x y, x + y, x  + y, -3x)

o32 : Ideal of R

i33 : I_*

        2   2     2   2           2
o33 = {x , y , x*y , x y, x + y, x  + y, -3x}

o33 : List

i34 : load "EthRoots.m2"

i35 : I = ideal(x^8, y^8, z^8)

              8   8   8
o35 = ideal (x , y , z )

o35 : Ideal of R

i36 : benchmark "ethRoot(1, I)"

o36 = .00138777625047801

o36 : RR (of precision 53)

i37 : R = kk[x_1..x_{10}]
stdio:36:9:(3): error: no method for binary operator _ applied to objects:
--            x (of class R)
--      _     1 (of class ZZ)

i38 : R = kk[x_1..x_9]
stdio:37:9:(3): error: no method for binary operator _ applied to objects:
--            x (of class R)
--      _     1 (of class ZZ)

i39 : R = kk[x_1,..,x_9]
stdio:38:12:(3): error: syntax error at '..'

i39 : l = {a..f}

o39 = {(a, b, c, d, e, f)}

o39 : List

i40 : kk[a..f]

o40 = kk[a, b, c, d, e, f]

o40 : PolynomialRing

i41 : R = kk[x_{1..10}]
stdio:40:9:(3): error: no method for binary operator _ applied to objects:
--            x (of class R)
--      _     {(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)} (of class List)

i42 : {1..100}

o42 = {(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      -------------------------------------------------------------------------
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
      -------------------------------------------------------------------------
      39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
      -------------------------------------------------------------------------
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
      -------------------------------------------------------------------------
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
      -------------------------------------------------------------------------
      93, 94, 95, 96, 97, 98, 99, 100)}

o42 : List

i43 : R = kk[x_{1} .. x_{100}]
stdio:42:9:(3): error: no method for binary operator _ applied to objects:
--            x (of class R)
--      _     {1} (of class List)

i44 : R = kk[x_0..x_3]
stdio:43:9:(3): error: no method for binary operator _ applied to objects:
--            x (of class R)
--      _     0 (of class ZZ)

i45 : x

o45 = x

o45 : R

i46 : R = kk[q_0..q_{100}]
stdio:45:11:(3): error: no method for binary operator .. applied to objects:
--            0 (of class ZZ)
--     ..     {100} (of class List)

i47 : R = kk[q_0 .. q_{100}]
stdio:46:12:(3): error: no method for binary operator .. applied to objects:
--            0 (of class ZZ)
--     ..     {100} (of class List)

i48 : R = kk[q_{0} .. q_{100}]

o48 = R

o48 : PolynomialRing

i49 : benchmark "ethRoot(1, ideal(q_0^8, q_1^8, q_2^8))"
currentString:1:51:(3):[5]: error: no method for binary operator ^ applied to objects:
--            q  (of class IndexedVariable)
--             0
--      ^     8 (of class ZZ)
currentString:1:51:(3):[5]: --entering debugger (type help to see debugger commands)
currentString:1:48-1:52: --source code:
timing scan(1, iBenchmark -> (ethRoot(1, ideal(q_0^8, q_1^8, q_2^8));null;null))

ii50 : benchmark "ethRoot(1, ideal((q_0)^8, (q_1)^8, (q_2)^8))"
currentString:1:53:(3):[5]: error: no method for binary operator ^ applied to objects:
--            q  (of class IndexedVariable)
--             0
--      ^     8 (of class ZZ)
currentString:1:8:(3):[3]: --back trace--

ii51 : break
stdio:48:1:(3): error: expected a list, hash table, or sequence

i52 : benchmark "ethRoot(1, ideal((q_0)^8, (q_1)^8, (q_2)^8))"
currentString:1:53:(3):[5]: error: no method for binary operator ^ applied to objects:
--            q  (of class IndexedVariable)
--             0
--      ^     8 (of class ZZ)
currentString:1:53:(3):[5]: --entering debugger (type help to see debugger commands)
currentString:1:49-1:54: --source code:
timing scan(1, iBenchmark -> (ethRoot(1, ideal((q_0)^8, (q_1)^8, (q_2)^8));null;null))

ii53 : q_0

oo53 = q
        0

oo53 : IndexedVariable

ii54 : q_0^8
stdio:2:4:(3): error: no method for binary operator ^ applied to objects:
--            q  (of class IndexedVariable)
--             0
--      ^     8 (of class ZZ)

ii55 : q^8_0
stdio:3:2:(3): error: no method for binary operator ^ applied to objects:
--            q (of class IndexedVariableTable)
--      ^     8 (of class ZZ)

ii56 : (q_0)

oo56 = q
        0

oo56 : IndexedVariable

ii57 : (q_0)^8
stdio:5:6:(3): error: no method for binary operator ^ applied to objects:
--            q  (of class IndexedVariable)
--             0
--      ^     8 (of class ZZ)

ii58 : x

oo58 = x

oo58 : kk[x, y, z]

ii59 : x^8

        8
oo59 = x

oo59 : kk[x, y, z]

ii60 : break
stdio:49:1:(3): error: expected a list, hash table, or sequence

i61 : R

o61 = R

o61 : PolynomialRing

i62 : describe R

o62 = kk[q   , q   , q   , q   , q   , q   , q   , q   , q   , q   , q    ,
          {0}   {1}   {2}   {3}   {4}   {5}   {6}   {7}   {8}   {9}   {10} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {11}   {12}   {13}   {14}   {15}   {16}   {17}   {18}   {19}   {20} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {21}   {22}   {23}   {24}   {25}   {26}   {27}   {28}   {29}   {30} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {31}   {32}   {33}   {34}   {35}   {36}   {37}   {38}   {39}   {40} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {41}   {42}   {43}   {44}   {45}   {46}   {47}   {48}   {49}   {50} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {51}   {52}   {53}   {54}   {55}   {56}   {57}   {58}   {59}   {60} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {61}   {62}   {63}   {64}   {65}   {66}   {67}   {68}   {69}   {70} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {71}   {72}   {73}   {74}   {75}   {76}   {77}   {78}   {79}   {80} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q    ,
       {81}   {82}   {83}   {84}   {85}   {86}   {87}   {88}   {89}   {90} 
                                                                           
      -------------------------------------------------------------------------
      q    , q    , q    , q    , q    , q    , q    , q    , q    , q     ,
       {91}   {92}   {93}   {94}   {95}   {96}   {97}   {98}   {99}   {100} 
                                                                            
      -------------------------------------------------------------------------
      Degrees => {101:1}, Heft => {1}, MonomialOrder => {MonomialSize => 32},
                                                        {GRevLex => {101:1}}
                                                        {Position => Up    }
      -------------------------------------------------------------------------
      DegreeRank => 1]

i63 : R_3

o63 = q
       {3}

o63 : R

i64 : R_3^8

       8
o64 = q
       {3}

o64 : R

i65 : benchmark "ethRoot(1, ideal(R_0^8, R_1^8, R_2^8, R_3^8))"

o65 = .0095925485804878

o65 : RR (of precision 53)

i66 : benchmark "ethRoot(1, ideal(R_0^8, R_1^8, R_2^8, R_3^8), EthRootStrategy => Substitution)"

o66 = .0543927473809525

o66 : RR (of precision 53)

i67 : apply({1..10}, i -> benchmark "ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution)"
        C-c C-c
i67 : apply({1..10}, i -> benchmark "ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution)"
        C-c C-c
i67 : apply({1..10}, i -> benchmark "ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution)")
currentString:1:54:(3):[6]: error: no method for binary operator * applied to objects:
--            7 (of class ZZ)
--      *     i (of class Symbol)
currentString:1:54:(3):[6]: --entering debugger (type help to see debugger commands)
currentString:1:53-1:55: --source code:
timing scan(1, iBenchmark -> (ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution);null;null))

ii68 : break
stdio:58:21:(3):[1]: error: expected a list, hash table, or sequence

i69 
      
      i

o69 = i

o69 : Symbol

i70 : {1..10}

o70 = {(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)}

o70 : List

i71 : apply({1..10}, i -> print i)
(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

o71 = {}

o71 : List

i72 : list{1..10}
stdio:64:1:(3): error: syntax error at 'list'

i72 : 1..10

o72 = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

o72 : Sequence

i73 : apply(1..10, i -> benchmark "ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution)")
currentString:1:54:(3):[6]: error: no method for binary operator * applied to objects:
--            7 (of class ZZ)
--      *     i (of class Symbol)
currentString:1:54:(3):[6]: --entering debugger (type help to see debugger commands)
currentString:1:53-1:55: --source code:
timing scan(1, iBenchmark -> (ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution);null;null))

ii74 : break
stdio:65:19:(3):[1]: error: expected a list, hash table, or sequence

i75 : apply(1..10, i -> print i)
1
2
3
4
5
6
7
8
9
10

o75 = (, , , , , , , , , )

o75 : Sequence

i76 : apply(1..10, i -> 2*i)

o76 = (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

o76 : Sequence

i77 : apply(1..10, i -> ethRoot(1, ideal(R_0^(7*i+1), R_1^(7*i+1), R_2^(7*i+1), R_3^(7*i+1)), EthRootStrategy => Substitution))

                                               2     2     2     2          
o77 = (ideal (q   , q   , q   , q   ), ideal (q   , q   , q   , q   ), ideal
               {3}   {2}   {1}   {0}           {3}   {2}   {1}   {0}        
      -------------------------------------------------------------------------
        3     3     3     3             4     4     4     4             5   
      (q   , q   , q   , q   ), ideal (q   , q   , q   , q   ), ideal (q   ,
        {3}   {2}   {1}   {0}           {3}   {2}   {1}   {0}           {3} 
      -------------------------------------------------------------------------
       5     5     5             6     6     6     6             7     7   
      q   , q   , q   ), ideal (q   , q   , q   , q   ), ideal (q   , q   ,
       {2}   {1}   {0}           {3}   {2}   {1}   {0}           {3}   {2} 
      -------------------------------------------------------------------------
       7     7             8     8     8     8             9     9     9   
      q   , q   ), ideal (q   , q   , q   , q   ), ideal (q   , q   , q   ,
       {1}   {0}           {3}   {2}   {1}   {0}           {3}   {2}   {1} 
      -------------------------------------------------------------------------
       9             10    10    10    10
      q   ), ideal (q   , q   , q   , q   ))
       {0}           {3}   {2}   {1}   {0}

o77 : Sequence

i78 : fn ( ZZ ) := n -> fn (n)
stdio:69:11:(3): error: expected left hand parameter to be a function, type, or a hash table

i79 : myfn ( ZZ ) := n -> myfn (n)
stdio:70:11:(3): error: expected left hand parameter to be a function, type, or a hash table

i80 : myfn (ZZ) := n -> (
      benchmark "ethRoot(1, ideal(R_0^(7*n+1), R_1^(7*n+1), R_2^(7*n+1), R_3^(7*n+1)), EthRootStrategy => Substitution)"
      )
stdio:71:11:(3): error: expected left hand parameter to be a function, type, or a hash table

i81 : myfn = method();

i82 : myfn (ZZ) := n -> (
      benchmark "ethRoot(1, ideal(R_0^(7*n+1), R_1^(7*n+1), R_2^(7*n+1), R_3^(7*n+1)), EthRootStrategy => Substitution)"
      )

o82 = {*Function[stdio:75:16-76:1]*}

o82 : FunctionClosure

i83 : apply(1..10, i -> myfn(i))
currentString:1:54:(3):[7]: error: no method for binary operator * applied to objects:
--            7 (of class ZZ)
--      *     n (of class Symbol)
currentString:1:54:(3):[7]: --entering debugger (type help to see debugger commands)
currentString:1:53-1:55: --source code:
timing scan(1, iBenchmark -> (ethRoot(1, ideal(R_0^(7*n+1), R_1^(7*n+1), R_2^(7*n+1), R_3^(7*n+1)), EthRootStrategy => Substitution);null;null))

ii84 : 