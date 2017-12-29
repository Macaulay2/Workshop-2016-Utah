loadPackage "NormalToricVarieties"

-- The following example arises from the Kreuzer-Skarke database, 
-- as h11=11 examples number 2 (indexed from 0).
-- Here is the entry:

{*
(4 11  M:22 11 N:16 8 H:11,19 [-16],    1   0   0  -1   0  -1  -1   1   1   1  -1)
                                        0   1   1   2   0  -1   0  -3  -2   1   3
                                        0   0   2   2   0   0   1  -4  -3   2   2
                                        0   0   0   0   1  -1  -1  -1  -1   1   1
*}
M = matrix {
    {1, 0, 0, -1, 0, -1, -1, 1, 1, 1, -1}, 
    {0, 1, 1, 2, 0, -1, 0, -3, -2, 1, 3}, 
    {0, 0, 2, 2, 0, 0, 1, -4, -3, 2, 2}, 
    {0, 0, 0, 0, 1, -1, -1, -1, -1, 1, 1}
    }

raylist = {
    {-1, 0, 0, 0}, 
    {0, 0, 0, -1}, 
    {1, 1, -1, -1}, 
    {0, 1, -1, 0}, 
    {-1, 1, -1, 1}, 
    {1, 1, 0, -1}, 
    {1, -1, 1, 1}, 
    {0, -1, 1, 0}, 
    {-1, -1, 1, -1}, 
    {-1, -1, 0, 2}, 
    {-1, -1, 0, 1}, 
    {1, 0, 0, 0}, 
    {0, 0, 0, 1}, 
    {0, 2, -1, -1}, 
    {2, 0, 1, -1}
    }

maxcones = {
    {0, 1, 3, 10},
    {0, 1, 3, 13},
    {0, 1, 8, 10},
    {0, 1, 8, 13},
    {0, 3, 4, 10},
    {0, 3, 4, 13},
    {0, 4, 9, 10},
    {0, 4, 9, 12},
    {0, 4, 12, 13},
    {0, 5, 7, 8},
    {0, 5, 7, 12},
    {0, 5, 8, 13},
    {0, 5, 12, 13},
    {0, 7, 8, 9},
    {0, 7, 9, 12},
    {0, 8, 9, 10},
    {1, 2, 3, 10},
    {1, 2, 3, 13},
    {1, 2, 5, 13},
    {1, 2, 5, 14},
    {1, 2, 10, 11},
    {1, 2, 11, 14},
    {1, 5, 8, 13},
    {1, 5, 8, 14},
    {1, 7, 8, 10},
    {1, 7, 8, 14},
    {1, 7, 10, 11},
    {1, 7, 11, 14},
    {2, 3, 5, 11},
    {2, 3, 5, 13},
    {2, 3, 9, 10},
    {2, 3, 9, 11},
    {2, 5, 11, 14},
    {2, 9, 10, 11},
    {3, 4, 9, 10},
    {3, 4, 9, 12},
    {3, 4, 12, 13},
    {3, 5, 11, 12},
    {3, 5, 12, 13},
    {3, 9, 11, 12},
    {5, 7, 8, 14},
    {5, 7, 12, 14},
    {5, 11, 12, 14},
    {6, 7, 9, 10},
    {6, 7, 9, 12},
    {6, 7, 10, 11},
    {6, 7, 11, 14},
    {6, 7, 12, 14},
    {6, 9, 10, 11},
    {6, 9, 11, 12},
    {6, 11, 12, 14},
    {7, 8, 9, 10}
    }

V = normalToricVariety(raylist, maxcones)

end--
restart
load "slow-cohom-example.m2"

-- IV = intersectionRing V
-- gens IV
ring V -- this works fine
cl V -- ZZ^11 is very quick.
wDiv V -- ZZ^15, immediate
---fromCDivToWDiv V -- this is slow, and returns something that looks (I hope!) wrong: huge numbers in it.
---elapsedTime pic V -- takes 56 seconds on my MBP, answer is ZZ^11 (as it should be)

orbits(V,0)
orbits(V,1)
orbits(V,2)
orbits(V,3)

F = OO V_1
---elapsedTime HH^1(V, F) -- takes awhile, how long?

elapsedTime assert isProjective V -- takes 10 seconds
assert isSimplicial V -- immediate
assert not isSmooth V -- immediate

