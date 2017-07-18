TEST /// 
ZZ/5[x,y,z];
F = 2*x^7*y^3*z^8+2*x^4*z^9+2*x*y^7*z^4;
-- old version 
time assert( nu(6,F) == 2968 )  
time assert( nuList(6,F) == {0, 4, 23, 118, 593, 2968} ) 
-- using powers
time assert( newNu(6,F,ComputePreviousNus=>false,TestFunction =>testPower) == 2968 )
time assert( newNuList(6,F,TestFunction =>testPower) == {0, 0, 4, 23, 118, 593, 2968} )
time assert( newNu(6,F,TestFunction =>testPower,SearchFunction=>binarySearchRecursive) == 2968 )
time assert( newNu(6,F,TestFunction =>testPower,SearchFunction=>linearSearch) == 2968 )
-- using roots (better)
time assert( newNu(6,F,ComputePreviousNus=>false) == 2968 )  
time assert( newNuList(6,F) == {0, 0, 4, 23, 118, 593, 2968} )  
time assert( newNu(6,F,SearchFunction=>binarySearchRecursive) == 2968 )  
time assert( newNu(6,F,SearchFunction=>linearSearch) == 2968 )  
-- using colon ideals (best)
time assert( newNuList(6,F,UseColonIdeals=>true) == {0, 0, 4, 23, 118, 593, 2968} )  
time assert( newNu(6,F,UseColonIdeals=>true,TestFunction=>testPower) == 2968 )
time assert( newNu(6,ideal F,UseColonIdeals=>true,TestFunction=>testPower,SearchFunction=>binarySearchRecursive) == 2968 )
time assert( newNu(6,ideal F,UseColonIdeals=>true,TestFunction=>testPower,SearchFunction=>linearSearch) == 2968 )

 
ZZ/17[x,y,z,w];
F = -5*x*y^4*z^3-2*x^4*z^3*w+6*y*z^3*w^4+7*z*w^3-6*w^2;
-- [slow] time assert( nu(2,F) == 220 ) -- old version 
-- [slow] time assert( newNu(2,F,Strategy=>Powers) == 220 )
time assert( newNu(2,F) == 220 ) 
time assert( newNu(2,F,SearchFunction=>binarySearchRecursive) == 220) 
time assert( newNu(2,F,SearchFunction=>linearSearch) == 220 ) 
time assert( newNu(2,F,UseColonIdeals=>true) == 220 )
time assert( newNu(2,F,UseColonIdeals=>true,SearchFunction=>binarySearchRecursive) == 220 )
time assert( newNu(2,ideal F,UseColonIdeals=>true,SearchFunction=>linearSearch) == 220 )

time assert( newNu(3,F) == 3756 )  

ZZ/3[x,y,z];
I=ideal(x^2+y^3,y^2+z^3,z^2+x^3);
-- old version
time assert( nu(3,I) == 39 ) -- lucky
time assert( nuList(3,I) == {3, 12, 39} ) 
-- powers 
time assert( newNu(3,I,ComputePreviousNus=>false) == 39 ) 
time assert( newNuList(3,I) == {0, 3, 12, 39} ) 
time assert( newNuList(3,I,UseColonIdeals=>true) == {0, 3, 12, 39} ) -- colon ignored
time assert( newNu(3,I,SearchFunction=>binarySearchRecursive) == 39 ) 
time assert( newNu(3,I,SearchFunction=>linearSearch) == 39 ) -- lucky

-- old version
time assert( nuHat(5,I) == 242 ) 
time assert( nuHatList(5,I) == {2, 8, 26, 80, 242} ) 
-- new version
time assert( newNuHat(5,I,ComputePreviousNus=>false) == 242 ) 
time assert( newNuHatList(5,I) == {0, 2, 8, 26, 80, 242} ) 
time assert( newNuHat(5,I,UseColonIdeals=>true) == 242 ) -- UseColonIdeals is ignored

ZZ/13[x,y,z];
I=ideal(x^3+y^4,y^6+z^3,z^5+x^7);
time assert( nu(1,I) == 11 )  -- old version
time assert( newNu(1,I) == 11 ) 
time assert( newNu(1,I,SearchFunction=>linearSearch) == 11 ) -- lucky
time assert( newNu(1,I,ComputePreviousNus=>false) == 11 ) 
time assert( newNu(1,I,TestFunction=>testRoot) == 11 ) -- bad choice

time assert( nuHat(2,I) == 154 ) -- slow
time assert( newNuHat(2,I) == 154 )  -- slow

ZZ/7[x,y,z];
I=ideal(x+y^2,y+z^2,z+x^2);
J=ideal(x^3,y^3,z^3);
time assert( nu(1,I,J) == 60 ) -- old version
time assert( newNu(1,I,J) == 60 ) -- slow
time assert( newNu(1,I,J,TestFunction=>testRoot) == 60 ) -- better
time assert( newNu(1,I,J,TestFunction=>testRoot,SearchFunction=>linearSearch) == 60 ) 

time assert( nuHat(1,I,J) == 48 ) -- old version
time assert( newNuHat(1,I,J,ComputePreviousNus=>false) == 48 )  
time assert( newNuHatList(1,I,J) == {6, 48} )  
time assert( newNuList(1,I,J,TestFunction=>testGenFrobeniusPower) == {6, 48} )  -- same as previous
///
