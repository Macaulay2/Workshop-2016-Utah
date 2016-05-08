
+ /Applications/Macaulay2-1.9/bin/M2 --no-readline --print-width 90
Macaulay2, version 1.9
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : f = () -> (
        local a;
        local b;
        a = 4;
        5)

o1 = f

o1 : FunctionClosure

i2 : f()

o2 = 5

i3 : a

o3 = a

o3 : Symbol

i4 : f = () -> ( local x ; QQ[x] )

o4 = f

o4 : FunctionClosure

i5 : f()
stdio:10:23:(3):[1]: error: encountered object not usable as variable at position 0 in list:
        null (of class Nothing)

i6 : f = () -> ( local x = local x ; QQ[x] )
stdio:12:21:(3): error: left hand side of assignment inappropriate

i7 : f = () -> ( local x; x = symbol x ; QQ[x] )

o7 = f

o7 : FunctionClosure

i8 : f()

o8 = QQ[x]

o8 : PolynomialRing

i9 : x

o9 = x

o9 : Symbol

i10 : use o8

o10 = QQ[x]

o10 : PolynomialRing

i11 : x

o11 = x

o11 : Symbol

i12 : f = () -> ( local x; x = symbol x ; QQ[x] ; x^3 )

o12 = f

o12 : FunctionClosure

i13 : f()

       3
o13 = x

o13 : QQ[x]

i14 : oo + x
stdio:20:4:(3): error: no method for binary operator + applied to objects:
--             3
--            x  (of class QQ[x])
--      +     x (of class Symbol)

i16 : x = (ring o13) _ 0

o16 = x

o16 : QQ[x]

i17 : o13 + x

       3
o17 = x  + x

o17 : QQ[x]

i18 : help genericMatrix 

o18 = genericMatrix -- make a generic matrix of variables
      ***************************************************

      Synopsis
      ========

        * Usage:genericMatrix(R,r,m,n)
        * Inputs:
            * R, a ring
            * r, a ring element, which is a variable in the ring R (this input is
              optional)
            * m, an integer
            * n, an integer
        * Outputs:
            * a matrix, with m rows and n columns whose entries are variables in the ring
              R starting with r

      Description
      ===========

      +---------------------------+
      |i1 : R = ZZ[a..z];         |
      +---------------------------+
      |i2 : genericMatrix(R,a,2,4)|
      |                           |
      |o2 = | a c e g |           |
      |     | b d f h |           |
      |                           |
      |             2       4     |
      |o2 : Matrix R  <--- R      |
      +---------------------------+
      |i3 : genericMatrix(R,i,3,2)|
      |                           |
      |o3 = | i l |               |
      |     | j m |               |
      |     | k n |               |
      |                           |
      |             3       2     |
      |o3 : Matrix R  <--- R      |
      +---------------------------+


      Omitting the input r is the same as having r be the first variable in R.

      +-------------------------+
      |i4 : genericMatrix(R,2,4)|
      |                         |
      |o4 = | a c e g |         |
      |     | b d f h |         |
      |                         |
      |             2       4   |
      |o4 : Matrix R  <--- R    |
      +-------------------------+
      |i5 : genericMatrix(R,3,2)|
      |                         |
      |o5 = | a d |             |
      |     | b e |             |
      |     | c f |             |
      |                         |
      |             3       2   |
      |o5 : Matrix R  <--- R    |
      +-------------------------+

      See also
      ========

        * "vars(Ring)" -- row matrix of the variables
        * "genericSkewMatrix" -- make a generic skew symmetric matrix of variables
        * "genericSymmetricMatrix" -- make a generic symmetric matrix

      Ways to use genericMatrix :
      ===========================

        * genericMatrix(Ring,RingElement,ZZ,ZZ)
        * genericMatrix(Ring,ZZ,ZZ)

o18 : DIV

i19 : ideal

o19 = ideal

o19 : CompiledFunctionClosure

i20 : methods ideal

o20 = {(ideal, List)         }
      {(ideal, Matrix)       }
      {(ideal, Module)       }
      {(ideal, MonomialIdeal)}
      {(ideal, Number)       }
      {(ideal, QuotientRing) }
      {(ideal, Ring)         }
      {(ideal, RingElement)  }
      {(ideal, Sequence)     }
      {(ideal, String)       }
      {(ideal, Variety)      }

o20 : VerticalList

i21 : ideal Ideal := I -> I

o21 = {*Function[stdio:27:18-27:21]*}

o21 : FunctionClosure

i22 : ideal ideal 3

o22 = ideal 3

o22 : Ideal of ZZ

i23 : code methods hilbertSeries
*** output flushed ***
i24 : hold 3

o24 = 3

o24 : Expression of class Holder

i25 : oo^4

       4
o25 = 3

o25 : Expression of class Power

i26 : hilbertSeries ( QQ[x,y] )

          1
o26 = --------
             2
      (1 - T)

o26 : Expression of class Divide

i27 : peek oo

                       2
o27 = Divide{1, (1 - T) }

i28 : peek o25

o28 = Power{3, 4}

i29 : peek'_5 o26

o29 = Divide{ZZ[T]{1}, Product{Power{ZZ[T]{1-T}, 2}}}

i30 : peek'_2 o26

                                      2
o30 = Divide{ZZ[T]{1}, Product{(1 - T) }}

i31 : peek'_3 o26

o31 = Divide{ZZ[T]{1}, Product{Power{1 - T, 2}}}

i32 : value o26
stdio:38:1:(3): error: not implemented : fraction fields of rings with inverses

i33 : m = matrix {{2,3},{4,5}}

o33 = | 2 3 |
      | 4 5 |

               2        2
o33 : Matrix ZZ  <--- ZZ

i34 : m_1

o34 = | 3 |
      | 5 |

        2
o34 : ZZ

i35 : help (symbol -, Set)
stdio:41:1:(3): error: documentation key for '- Set' encountered, but no method installed

i36 : viewHelp Set

i37 : {3}-{4}

o37 = {-1}

o37 : List

i38 : {3}-set {4}

o38 = {3}

o38 : List

i39 : set {3}-set {4}

o39 = set {3}

o39 : Set

i40 : set {3}- {4}

o40 = set {3}

o40 : Set

i41 : viewHelp "Depth::Depth"

i42 : apropos "oaded"

o42 = {loadedFiles, loadedPackages}

o42 : List

i43 : loadedFiles

o43 = MutableHashTable{...106...}

o43 : MutableHashTable

i44 : peek oo

o44 = MutableHashTable{0 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/command.m2                                     }
                       1 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/classes.m2
                       2 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/option.m2
                       3 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/methods.m2
                       4 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/profile.m2
                       5 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/autoload.m2
                       6 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/system.m2
                       7 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/debugging.m2
                       8 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/remember.m2
                       9 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/set.m2
                       10 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/files.m2
                       11 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/fold.m2
                       12 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/max.m2
                       13 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/structure.m2
                       14 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/combinatorics.m2
                       15 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/lists.m2
                       16 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/nets.m2
                       17 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/robust.m2
                       18 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/content.m2
                       19 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/html0.m2
                       20 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/validate.m2
                       21 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/expressions.m2
                       22 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/peek.m2
                       23 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/printing.m2
                       24 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/gateway.m2
                       25 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/rings.m2
                       26 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/integers.m2
                       27 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/engine.m2
                       28 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/enginering.m2
                       29 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/rationals.m2
                       30 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/reals.m2
                       31 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/quotient.m2
                       32 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/powers.m2
                       33 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/orderedmonoidrings.m2
                       34 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/variables.m2
                       35 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/indeterminates.m2
                       36 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/ofcm.m2
                       37 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/tables.m2
                       38 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/modules.m2
                       39 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/matrix.m2
                       40 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/matrix1.m2
                       41 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/mutablemat.m2
                       42 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/quotring.m2
                       43 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/multilin.m2
                       44 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/flint.m2
                       45 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/genmat.m2
                       46 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/modules2.m2
                       47 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/gb.m2
                       48 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/matrix2.m2
                       49 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/colon.m2
                       50 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/galois.m2
                       51 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/ringmap.m2
                       52 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/newring.m2
                       53 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/matrix3.m2
                       54 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/ext.m2
                       55 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/tor.m2
                       56 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/duals.m2
                       57 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/gradedmodules.m2
                       58 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/chaincomplexes.m2
                       59 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/res.m2
                       60 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/monideal.m2
                       61 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/radical.m2
                       62 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/factor.m2
                       63 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/integrate.m2
                       64 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/http.m2
                       65 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/minPres.m2
                       66 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/monomcurve.m2
                       67 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/fano.m2
                       68 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/schubert.m2
                       69 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/code.m2
                       70 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/dotdot.m2
                       71 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/local.m2
                       72 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/packages.m2
                       73 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/document.m2
                       74 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/hypertext.m2
                       75 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/texhtml.m2
                       76 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/html.m2
                       77 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/varieties.m2
                       78 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/mathml.m2
                       79 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/texmacs.m2
                       80 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/pretty.m2
                       81 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/undoc.m2
                       82 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/obsolete.m2
                       83 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/exports.m2
                       84 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/tvalues.m2
                       85 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/typicalvalues.m2
                       86 => /Applications/Macaulay2-1.9/share/Macaulay2/Text.m2
                       87 => /Applications/Macaulay2-1.9/share/Macaulay2/SimpleDoc.m2
                       88 => /Applications/Macaulay2-1.9/share/Macaulay2/Elimination.m2
                       89 => /Applications/Macaulay2-1.9/share/Macaulay2/LLLBases.m2
                       90 => /Applications/Macaulay2-1.9/share/Macaulay2/PrimaryDecomposition/GTZ.m2
                       91 => /Applications/Macaulay2-1.9/share/Macaulay2/PrimaryDecomposition/Shimoyama-Yokoyama.m2
                       92 => /Applications/Macaulay2-1.9/share/Macaulay2/PrimaryDecomposition/Eisenbud-Huneke-Vasconcelos.m2
                       93 => /Applications/Macaulay2-1.9/share/Macaulay2/PrimaryDecomposition/radical.m2
                       94 => /Applications/Macaulay2-1.9/share/Macaulay2/PrimaryDecomposition.m2
                       95 => /Applications/Macaulay2-1.9/share/Macaulay2/ReesAlgebra.m2
                       96 => /Applications/Macaulay2-1.9/share/Macaulay2/IntegralClosure.m2
                       97 => /Applications/Macaulay2-1.9/share/Macaulay2/Parsing.m2
                       98 => /Applications/Macaulay2-1.9/share/Macaulay2/Classic.m2
                       99 => /Applications/Macaulay2-1.9/share/Macaulay2/TangentCone.m2
                       100 => /Applications/Macaulay2-1.9/share/Macaulay2/ConwayPolynomials.m2
                       101 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/last.m2
                       102 => /Applications/Macaulay2-1.9/share/Macaulay2/Core/setup.m2
                       103 => /Users/dan/Library/Application Support/Macaulay2/init.m2
                       104 => /Applications/Macaulay2-1.9/share/Macaulay2/Macaulay2Doc.m2
                       105 => /Applications/Macaulay2-1.9/share/Macaulay2/Depth.m2

i45 : uninstallAllPackages()
-- uninstalling package NumericalAlgebraicGeometry
