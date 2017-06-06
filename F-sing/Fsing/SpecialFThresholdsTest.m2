-- binaryFormFPT test 1
TEST /// 
R = ZZ/2[x,y]
assert( binaryFormFPT( x ) == 1 )
assert( binaryFormFPT( x^10*y^15 ) == 1/15 )
assert( binaryFormFPT( x^10*y^3+x^9*y^4+x^6*y^7+x^4*y^9+x^3*y^10+x*y^12+y^13 ) == 315/2048 )
assert( binaryFormFPT( x^10*y+x^9*y^2+x^8*y^3+x^7*y^4+x^4*y^7+x^3*y^8+x*y^10 ) == 93/512 )
assert( binaryFormFPT( x^10+x^8*y^2+x^5*y^5+x^4*y^6+x^3*y^7 ) == 1/5 )
///

-- binaryFormFPT test 2
TEST /// 
R = ZZ/31[x,y]
L = { x, y, x+3*y, x+10*y }
assert( binaryFormFPT( L, { 3, 6, 6, 10 } ) == 221645/2770563 )
assert( binaryFormFPT( L, { 5, 17, 20, 11 } ) == 160922837253/4264455187205 )
assert( binaryFormFPT( L, { 5, 12, 8, 11 } ) == 1/18 )
assert( binaryFormFPT( L, { 19, 14, 17, 18 } ) == 14198854736226229/482761061031691789 )
assert( binaryFormFPT( L, { 17, 10, 17, 20 } ) == 1/32 )
///

-- binaryFormFPT test 3
TEST /// 
kk = GF( ZZ/5[a]/ideal( a^3+a+1 ) )
R = kk[x,y]
L = { x, y, x+y, x-y, x+a*y, x+a^2*y }
assert( binaryFormFPT( L, { 4, 19, 2, 11, 5, 20 } ) == 30535166380835361168/931322574615478515625 )
assert( binaryFormFPT( L, { 22, 76, 46, 30, 92, 88 } ) == 18466082398285089576311704129149/3268496584496460855007171630859375 )
assert( binaryFormFPT( L, { 3, 32, 2, 32, 73, 84 } ) == 18765116046319289672094923747611764472226361925425255193119555448/2120458113234079732946726383480129385361578897573053836822509765625 )
///



