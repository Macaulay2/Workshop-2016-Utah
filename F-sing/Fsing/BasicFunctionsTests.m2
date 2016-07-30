-- denom test#1 - (Rational)
TEST ///
time J = denom(17/23);
assert(J == 23)
///

-- denom test#2 - (Integer)
TEST ///
time J = denom(17);
assert(J == 1)
///

-- num test#1 - (Rational)
TEST ///
time J = num(23/17);
assert(J == 23)
///

-- num test#2 - (Integer)
TEST ///
time J = num(17);
assert(J == 17)
///

-- fracPart test#1
TEST ///
time J = fracPart(15/7);
assert(J == 1/7)
///

-- floorLog test#1 
TEST ///
time J = floorLog(2,3);
assert(J == 1)
///

-- multOrder test#1 
TEST ///
time J = multOrder(10, 7);
assert(J == 6)
///

-- divideFraction test#1 - (denominator a power of p)
TEST ///
time J = divideFraction(7, 19/49);
assert(J == {19, 2, 0} )
///

-- divideFraction test#2 - (denominator not power of p)
TEST ///
time J = divideFraction(2, 5/56);
assert(J == {5, 3, 3} )
///

-- divideFraction test#3 - (denominator not power of p)
TEST ///
time J = divideFraction(2, 5/24);
assert(J == {5, 3, 2} )
///


-- divideFraction test#4 - (negative)
TEST ///
time J = divideFraction(7, -19/49);
assert(J == { -19, 2, 0} )
///

-- findNearPthPowerAbove test#1 
TEST ///
time J = findNearPthPowerAbove(5,2,1/25);
assert(J == 1/25 )
///

-- findNearPthPowerAbove test#2 
TEST ///
time J = findNearPthPowerAbove(5,2,1/27);
assert(J == 1/25 )
///

-- findNearPthPowerBelow test#1 
TEST ///
time J = findNearPthPowerBelow(7,2,1/49);
assert(J == 1/49 )
///

-- findNearPthPowerBelow test#2 
TEST ///
time J = findNearPthPowerBelow(7,2,1/47);
assert(J == 1/49 )
///

-- nontrivialPowerSet test#1 
TEST ///
L = {a, b, c};
time J = nontrivialPowerSet(L);
assert(J == {{a}, {b}, {a, b}, {c}, {a, c}, {b, c}, {a, b, c}})
///

-- numberToPrimeFactorList test#1 - (no multiplicity)
TEST ///
time J = numberToPrimeFactorList(210);
assert(J == {2, 3, 5, 7})
///

-- numberToPrimeFactorList test#2 - (multiplicity)
TEST ///
time J = numberToPrimeFactorList(136);
assert(J == {2, 2, 2, 17})
///

-- getFactorList test#1 - (prim)
TEST ///
time J = getFactorList(13);
assert(J == {13})
///

-- getFactorList test#2 - (composite)
TEST ///
time J = getFactorList(12);
assert(J == {2, 3, 12, 4, 6})
///

-- findNumberBetweenWithDenom test1
TEST ///
time J = findNumberBetweenWithDenom(3, 0, 1);
assert(J == {0, 1/3, 2/3, 1})
///

-- findNumberBetweenWithDenom test2
TEST ///
time J = findNumberBetweenWithDenom(4, 1, 2);
assert(J == {1, 5/4, 3/2, 7/4, 2})
///

-- findNumberBetween test1
TEST ///
time J = findNumberBetween(3, 1, 2);
assert(J == {0, 4/3, 3/2, 5/3, 1})
///

-- findNumberBetween test2
TEST ///
time J = findNumberBetween(4, 1, 2);
assert(J == {1, 5/4, 4/3, 3/2, 5/3, 7/4, 2})
///

-- digit test#1
TEST ///
time J = digit(3, 2, 3/4);
assert(J == 0)
///

-- digit test#2
TEST ///
time J = digit(3, 1, 3/4);
assert(J == 2)
///

-- digit test#3
TEST ///
time J = digit(5, 3, 1/13);
assert(J == 4)
///

-- digit test#4
TEST ///
L = {3/4, 1/13};
time J = digit(5, 3, L);
assert(J == {3,4})
///

-- basePExp test#1 
TEST ///
time J = basePExp(2, 22);
assert(J == {0, 1, 1, 0, 1})
///

-- basePExp test#2 
TEST ///
time J = basePExp(5, 399);
assert(J == {4, 4, 0, 3})
///

-- basePExp test#3 
TEST ///
time J = basePExp(2, 4, 1/5);
assert(J == {0, 0, 1, 1})
///

-- basePExp test#4 
TEST ///
time J = basePExp(7, 7, 1/19);
assert(J == {0, 2, 4, 0, 2, 4,0})
///

-- truncatedBasePExp test#1
TEST ///
time J = truncatedBasePExp(7, 4, 1/19);
assert(J == 18/343)
///

-- truncatedBasePExp test#2
TEST ///
time J = truncatedBasePExp(7, 4, 1/29);
assert(J == 82/2401)
///

-- truncatedBasePExp test#3
TEST ///
time J = truncatedBasePExp(7, 4, {1/19, 1/29});
assert(J == {18/343, 82/2401)
///

-- baseP1 - NEEED

-- carryTest - test#1
TEST ///
time J = carryTest(2, {1/4, 1/4});
assert(J == 3)
///


-- firstCarry test#1
TEST ///
time J = firstCarry(2, {1/4, 1/4});
assert(J == 3)
///

-- firstCarry test#2
TEST ///
time J = firstCarry(2, {1/2, 1/2});
assert(J == 2)
///

-- reciprocal test#1 - (integer)
TEST ///
V = {2, 3, 19};
time J = reciprocal(V);
assert(J == {1/2, 1/3, 1/19})
///

-- reciprocal test#2 - (rational)
TEST ///
V = {2, 1/3, 19};
time J = reciprocal(V);
assert(J == {1/2, 3, 1/19})
///

-- getCanVector test#1 - (list)
TEST ///
time J = getCanVector(3, 7);
assert(J == {0, 0, 0, 1, 0, 0, 0})
///

-- getNumAndDenom test#1
TEST ///
V = {2, -1/12, 1/5};
time J = getNumAndDenom(V);
assert(J == ({120, -5, 12}, 60))
///

-- taxicabNorm test#1 - (integer)
TEST ///
V = {2, -2, 3};
time J = taxicabNorm(V);
assert(J == 7)
///

-- taxicabNorm test#2 - (rational)
TEST ///
V = {1, 1/2, 3};
time J = taxicabNorm(V);
assert(J == 9/2)
///

-- xInt test#1 
TEST ///
time J = xInt(1, 1, -2, -2);
assert(J == 0)
///

-- allPartitions test#1
TEST ///
time J = allPartitions(3, 2);
assert(J == {{2, 1}, {1, 2}})
///

-- maxIdeal test#1
TEST ///
R = QQ[x, y, z]
time J = maxIdeal(R);
assert(J == monomialIdeal(x, y, z))
///
