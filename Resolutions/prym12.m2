R = (ZZ/27901)[a, b, c, d, e, f, g, h, i, j, k]
I = ideal(
    e*h-5132*f*h-13016*g*h-8896*h^2-9659*a*i-7916*b*i+11904*c*i+7166*d*i+12515*e*i-6936*f*i-11695*g*i+8052*h*i+439*i^2-5256*a*j-7608*b*j+3622*c*j-11212*d*j+7503*e*j-9970*f*j-2380*g*j-4553*h*j+1825*i*j+7709*j^2-4296*a*k-10789*b*k+3829*c*k+12449*d*k+2660*e*k-5672*f*k+4114*g*k-1825*h*k-7709*i*k,
    d*h-5139*f*h+7362*g*h-468*h^2-10524*a*i+7692*b*i+11832*c*i+13283*d*i-4347*e*i+6700*f*i+4364*g*i-2751*h*i+7510*i^2-5752*a*j-11756*b*j-4453*c*j+7930*d*j+10761*e*j-178*f*j-5727*g*j+12562*h*j-877*i*j+9745*j^2-77*a*k-8830*b*k+1556*c*k+3078*d*k-3718*e*k+8478*f*k+7829*g*k+877*h*k-9745*i*k,
    c*h+5566*f*h+7750*g*h-11229*h^2-10136*a*i-10974*b*i-12781*c*i+8152*d*i+4989*e*i-8005*f*i-11263*g*i-8889*h*i+4649*i^2+13147*a*j-11552*b*j+240*c*j-2541*d*j+166*e*j-11708*f*j+925*g*j+6765*h*j-3131*i*j+5010*j^2-3568*a*k-8392*b*k-8014*c*k+89*d*k+6299*e*k+7964*f*k-11414*g*k+3131*h*k-5010*i*k,
    b*h+13255*f*h+11575*g*h+11481*h^2-9025*a*i+497*b*i+13891*c*i-1017*d*i-10492*e*i-8007*f*i-1058*g*i-3790*h*i-10111*i^2+8526*a*j-5915*b*j+6238*c*j-9504*d*j-3367*e*j+2426*f*j-622*g*j+11514*h*j-10867*i*j+6070*j^2-7976*a*k-5221*b*k+6741*c*k-201*d*k-12849*e*k+4412*f*k-1403*g*k+10867*h*k-6070*i*k,
    a*h-4726*f*h+1269*g*h-13716*h^2-168*a*i-8618*b*i-12378*c*i-4174*d*i+10193*e*i+436*f*i-8525*g*i-9848*h*i+8850*i^2-6555*a*j+700*b*j+9050*c*j-6055*d*j-5300*e*j+4307*f*j-12539*g*j+8329*h*j-6255*i*j+6213*j^2+11678*a*k-4876*b*k+588*c*k+3595*d*k-9967*e*k-5514*f*k+10722*g*k+6255*h*k-6213*i*k,
    g^2-3537*f*h+10135*g*h+1188*h^2+13412*a*i+10766*b*i-2880*c*i+11466*d*i+691*e*i-10871*f*i+13398*g*i+10968*h*i-11273*i^2-1608*a*j-5549*b*j-4877*c*j+6789*d*j+11495*e*j-8489*f*j+2451*g*j+1065*h*j+4938*i*j-4220*j^2+8429*a*k-6589*b*k-3944*c*k-10759*d*k-6097*e*k-13419*f*k+10208*g*k-4938*h*k+4220*i*k,
    f*g+5356*f*h+9025*g*h+11215*h^2+8422*a*i-6947*b*i+12104*c*i-5980*d*i-8429*e*i-1876*f*i+10057*g*i+6767*h*i-6979*i^2+6138*a*j+13658*b*j+11921*c*j-6559*d*j-2142*e*j-10233*f*j+9527*g*j+11640*h*j+6785*i*j-9304*j^2+2139*a*k-5942*b*k+9632*c*k-5007*d*k-11039*e*k+11607*f*k-4661*g*k-6785*h*k+9304*i*k,
    e*g+3214*f*h-10526*g*h-4683*h^2+12001*a*i-5425*b*i+2684*c*i-11624*d*i+5336*e*i-1850*f*i+488*g*i-7985*h*i+6860*i^2-1847*a*j-11450*b*j-13708*c*j+8894*d*j-8323*e*j+10876*f*j+1498*g*j-10803*h*j+12170*i*j+6866*j^2+8765*a*k-2569*b*k+10457*c*k-7202*d*k-6681*e*k+6487*f*k+3943*g*k-12170*h*k-6866*i*k,
    d*g-7723*f*h+4538*g*h-10829*h^2+12556*a*i-418*b*i+5894*c*i-5473*d*i-10709*e*i+3994*f*i+5870*g*i-5814*h*i-9730*i^2-3083*a*j+6737*b*j+299*c*j+2165*d*j+11036*e*j-10789*f*j-2930*g*j+6330*h*j-6510*i*j+6658*j^2-12631*a*k+5174*b*k-11634*c*k+8333*d*k-12153*e*k+8744*f*k+3400*g*k+6510*h*k-6658*i*k,
    c*g+6596*f*h+782*g*h+4279*h^2-10060*a*i-8444*b*i+3326*c*i-8591*d*i-7778*e*i-2065*f*i+6709*g*i+1642*h*i-11019*i^2+878*a*j-4426*b*j-11275*c*j+1418*d*j-143*e*j-12319*f*j-3465*g*j-4693*h*j+6241*i*j-6070*j^2+1100*a*k-8035*b*k-236*c*k+1426*d*k+1331*e*k+1823*f*k-12189*g*k-6241*h*k+6070*i*k,
    b*g+13255*f*h-10653*g*h+8325*h^2-4466*a*i+6544*b*i-1059*c*i+2777*d*i+3101*e*i-12929*f*i-2274*g*i-4205*h*i-10454*i^2-11963*a*j-10967*b*j+10598*c*j-11031*d*j-5112*e*j+146*f*j-12796*g*j-2438*h*j-678*i*j-4924*j^2+12026*a*k-13375*b*k-5325*c*k+793*d*k-6197*e*k-10900*f*k+12892*g*k+678*h*k+4924*i*k,
    a*g+1741*f*h-1836*g*h+6290*h^2+12853*a*i-13501*b*i-7836*c*i-9557*d*i+1164*e*i-6709*f*i+4621*g*i-1201*h*i-9988*i^2+4310*a*j+2621*b*j+12149*c*j-7550*d*j-6891*e*j-4291*f*j+6049*g*j+117*h*j-57*i*j-207*j^2+5215*a*k-2592*b*k+4645*c*k-12465*d*k-6620*e*k-4848*f*k+9871*g*k+57*h*k+207*i*k,
    f^2+11990*f*h+11882*g*h+129*h^2+8899*a*i+9786*b*i+5295*c*i+5429*d*i+10938*e*i+1286*f*i+4551*g*i+6353*h*i+8074*i^2+6269*a*j-10928*b*j-6130*c*j-5714*d*j-6794*e*j-7038*f*j+11854*g*j-9469*h*j-2792*i*j+3258*j^2+5632*a*k+701*b*k+10687*c*k-6374*d*k+2358*e*k+9694*f*k+1395*g*k+2792*h*k-3258*i*k,
    e*f-10719*f*h+226*g*h-2830*h^2-13261*a*i+9538*b*i+12536*c*i-616*d*i+12738*e*i+9566*f*i+11533*g*i+7011*h*i-9696*i^2-4780*a*j+10125*b*j-5537*c*j+9455*d*j+3043*e*j+6408*f*j+8845*g*j-356*h*j+12102*i*j-3179*j^2+5240*a*k+6153*b*k-11474*c*k-12835*d*k+12790*e*k+12045*f*k+10052*g*k-12102*h*k+3179*i*k,
    d*f-6829*f*h+13737*g*h+11019*h^2-1776*a*i+1446*b*i+2186*c*i-7709*d*i+5185*e*i+7678*f*i+13697*g*i+549*h*i-13198*i^2-2042*a*j+6960*b*j+370*c*j+11926*d*j+4647*e*j-10608*f*j-13201*g*j-6404*h*j+13526*i*j-5428*j^2-9146*a*k+7339*b*k-10282*c*k+1839*d*k+13793*e*k+12652*f*k-8299*g*k-13526*h*k+5428*i*k,
    c*f-5480*f*h+9800*g*h-8773*h^2-8715*a*i+12623*b*i+10989*c*i-12999*d*i+10925*e*i-2158*f*i-2497*g*i+11351*h*i+3718*i^2+8330*a*j+5609*b*j-1827*c*j+11303*d*j+7534*e*j-4753*f*j+8651*g*j+4174*h*j-7925*i*j+6149*j^2+11303*a*k-13075*b*k+11153*c*k+12725*d*k-11878*e*k+7899*f*k-7892*g*k+7925*h*k-6149*i*k,
    b*f+8437*f*h-1516*g*h+5061*h^2+12824*a*i+5189*b*i+3498*c*i-9442*d*i-11485*e*i-3901*f*i+1609*g*i+1115*h*i-7059*i^2+7350*a*j-1757*b*j-1550*c*j-4391*d*j-10425*e*j+881*f*j+9865*g*j-4950*h*j+8441*i*j-4803*j^2-1741*a*k+10992*b*k+7439*c*k-12059*d*k-7551*e*k-10980*f*k+12009*g*k-8441*h*k+4803*i*k,
    a*f+9788*f*h-3127*g*h-11921*h^2-7837*a*i-249*b*i+9314*c*i+12677*d*i+1569*e*i-2841*f*i+4918*g*i-13535*h*i+4659*i^2+9843*a*j+489*b*j+9436*c*j-4608*d*j+3695*e*j+850*f*j-12192*g*j+852*h*j-7014*i*j+11518*j^2-9803*a*k+5788*b*k-6749*c*k+2273*d*k+6153*e*k-2174*f*k-5511*g*k+7014*h*k-11518*i*k,
    e^2-880*f*h+8553*g*h-3715*h^2+6547*a*i+5513*b*i+9225*c*i+13033*d*i+459*e*i+13281*f*i+13243*g*i-10517*h*i-10141*i^2-159*a*j+3706*b*j-5818*c*j+7006*d*j+211*e*j-1222*f*j+2835*g*j-5258*h*j-9402*i*j+1058*j^2-12931*a*k-7215*b*k-6585*c*k+5856*d*k-8306*e*k+7682*f*k-12502*g*k+9402*h*k-1058*i*k,
    d*e+10937*f*h+5870*g*h+4838*h^2-2495*a*i+5740*b*i-8689*c*i-13672*d*i-2581*e*i-437*f*i-7966*g*i-2872*h*i+9998*i^2-929*a*j+5924*b*j-7334*c*j+11392*d*j-11801*e*j-9334*f*j+6719*g*j-12368*h*j+7802*i*j+5380*j^2+2765*a*k-6895*b*k+8153*c*k+6368*d*k+12462*e*k-3847*f*k+2370*g*k-7802*h*k-5380*i*k,
    c*e+1191*f*h+11342*g*h-3791*h^2-10634*a*i+1394*b*i-5950*c*i+5495*d*i-11427*e*i+3623*f*i+7209*g*i-7943*h*i-4322*i^2+10773*a*j-596*b*j+3820*c*j+12893*d*j+10281*e*j+5015*f*j+5525*g*j-2611*h*j+4505*i*j+13279*j^2+6546*a*k-9315*b*k-2657*c*k+2655*d*k-8433*e*k+2418*f*k+6933*g*k-4505*h*k-13279*i*k,
    b*e+5887*f*h+4669*g*h+1105*h^2+11457*a*i+5807*b*i-8987*c*i+6795*d*i-309*e*i+7686*f*i-9966*g*i-7586*h*i-10006*i^2+6615*a*j+8669*b*j+4835*c*j-2724*d*j-11695*e*j-929*f*j+12621*g*j-8240*h*j+1816*i*j+8220*j^2+318*a*k-11630*b*k-2854*c*k-660*d*k+9790*e*k-5035*f*k-9655*g*k-1816*h*k-8220*i*k,
    a*e-3282*f*h-11686*g*h-3168*h^2-5583*a*i-2321*b*i-8267*c*i+3086*d*i-11003*e*i-9621*f*i-2000*g*i+10611*h*i-4281*i^2-9304*a*j-11217*b*j+9357*c*j-11521*d*j-10643*e*j+12699*f*j+7776*g*j+8816*h*j+630*i*j+4285*j^2-8417*a*k-12443*b*k-2095*c*k+4049*d*k-7531*e*k+9514*f*k-4535*g*k-630*h*k-4285*i*k,
    d^2+4712*f*h+13657*g*h+11942*h^2+4904*a*i+10402*b*i+5345*c*i+13911*d*i-9862*e*i+5338*f*i-13423*g*i-13495*h*i-13012*i^2-6510*a*j-5453*b*j+9243*c*j+3932*d*j+10801*e*j-7252*f*j-10892*g*j-2861*h*j+12031*i*j+13878*j^2+108*a*k+4747*b*k+1218*c*k-1895*d*k+8733*e*k-3514*f*k-12028*g*k-12031*h*k-13878*i*k,
    c*d-3881*f*h+5326*g*h+5633*h^2+12853*a*i-6371*b*i-2600*c*i-12262*d*i+6659*e*i-13052*f*i-8608*g*i-11508*h*i-4985*i^2+9391*a*j+9873*b*j+12124*c*j+13204*d*j+9863*e*j-1667*f*j-10637*g*j-6649*h*j+12008*i*j-4044*j^2-7273*a*k+138*b*k+11919*c*k-2137*d*k+4642*e*k-5756*f*k+11634*g*k-12008*h*k+4044*i*k,
    b*d-10844*f*h-3738*g*h+11110*h^2-11126*a*i+6464*b*i-3697*c*i+3377*d*i+11285*e*i-10392*f*i-11678*g*i+8870*h*i-8992*i^2-2327*a*j+1403*b*j-13553*c*j+13614*d*j+622*e*j-13117*f*j+7718*g*j+6743*h*j+4273*i*j+10758*j^2+2294*a*k+10176*b*k+13846*c*k+13508*d*k+13685*e*k+11313*f*k+2249*g*k-4273*h*k-10758*i*k,
    a*d+1604*f*h-7198*g*h-2478*h^2+13214*a*i+13817*b*i-8290*c*i-5575*d*i+1565*e*i-13147*f*i-4885*g*i+4498*h*i-352*i^2-12875*a*j-8171*b*j-1413*c*j+9564*d*j+6371*e*j+2031*f*j+4894*g*j-12004*h*j+2125*i*j-9321*j^2-11440*a*k+6988*b*k-12733*c*k-13927*d*k+5332*e*k-9392*f*k+12356*g*k-2125*h*k+9321*i*k,
    c^2-3048*f*h-9313*g*h-13872*h^2-2613*a*i-1630*b*i-6900*c*i+8251*d*i+13440*e*i-9262*f*i+2809*g*i-7573*h*i+7742*i^2-12354*a*j-12249*b*j+1021*c*j+6520*d*j-2608*e*j+53*f*j-5407*g*j-9869*h*j+10481*i*j-1297*j^2-8752*a*k-9272*b*k+10989*c*k-6718*d*k+11010*e*k+12980*f*k+2127*g*k-10481*h*k+1297*i*k,
    b*c-4588*f*h-10958*g*h+918*h^2+1622*a*i-2965*b*i+12023*c*i-9826*d*i+8461*e*i-12095*f*i-832*g*i-4941*h*i-12853*i^2-4607*a*j-3516*b*j-3136*c*j+10445*d*j-6699*e*j-10564*f*j+5285*g*j+3351*h*j+2470*i*j-2350*j^2-8507*a*k+12962*b*k+13583*c*k+1851*d*k+10478*e*k-344*f*k+9502*g*k-2470*h*k+2350*i*k,
    a*c-11023*f*h-6716*g*h+5378*h^2+3917*a*i+11961*b*i+3536*c*i-2032*d*i-1719*e*i+10929*f*i-5049*g*i+3748*h*i-12613*i^2-11575*a*j+11828*b*j+8481*c*j+10845*d*j-4674*e*j+9911*f*j+12137*g*j+2622*h*j-4434*i*j+12471*j^2+12537*a*k-6449*b*k+1897*c*k+461*d*k-10240*e*k+12016*f*k+9991*g*k+4434*h*k-12471*i*k,
    b^2-623*f*h+10999*g*h+12565*h^2-5843*a*i-4734*b*i+11914*c*i+13322*d*i-10143*e*i-7448*f*i+6066*g*i+886*h*i-5350*i^2+13625*a*j-7960*b*j-2609*c*j+5796*d*j+11744*e*j-9916*f*j-9281*g*j+13245*h*j-1963*i*j+11551*j^2-3954*a*k-10713*b*k+4970*c*k+12606*d*k-8715*e*k+8395*f*k-7895*g*k+1963*h*k-11551*i*k,
    a*b-1095*f*h-13537*g*h+5098*h^2-11872*a*i-1097*b*i-3172*c*i+1021*d*i-33*e*i-151*f*i+12984*g*i-8176*h*i+5333*i^2-5485*a*j+1988*b*j+13243*c*j+11980*d*j-5504*e*j+13131*f*j-4426*g*j-10459*h*j-8717*i*j-6505*j^2+1184*a*k+13637*b*k-10852*c*k-8709*d*k-3312*e*k+12602*f*k+5126*g*k+8717*h*k+6505*i*k,
    a^2+12546*f*h-12024*g*h+11614*h^2-4412*a*i+402*b*i-897*c*i+11820*d*i+1847*e*i-11944*f*i+3796*g*i-13819*h*i-10673*i^2+336*a*j+3354*b*j+5476*c*j-6879*d*j-5535*e*j-12383*f*j+4203*g*j-13313*h*j+10234*i*j-7223*j^2-2457*a*k+10605*b*k-7514*c*k+1602*d*k-3027*e*k+9616*f*k-3915*g*k-10234*h*k+7223*i*k,
    i^3+9581*f*h*j-83*g*h*j+6467*h^2*j-9425*a*i*j+5634*b*i*j+9879*c*i*j+2888*d*i*j-4570*e*i*j-7722*f*i*j+4679*g*i*j-3397*h*i*j-7849*i^2*j+3292*a*j^2+13030*b*j^2+2232*c*j^2+8943*d*j^2-4548*e*j^2-10278*f*j^2-11572*g*j^2-10707*h*j^2-9295*i*j^2-7822*j^3+105*f*h*k+10954*g*h*k+10342*h^2*k-12399*a*i*k+4043*b*i*k+10596*c*i*k-6047*d*i*k-3041*e*i*k+8392*f*i*k-7848*g*i*k-10924*h*i*k-10866*i^2*k-1520*a*j*k-4109*b*j*k-7759*c*j*k+3086*d*j*k-5173*e*j*k+8123*f*j*k-9877*g*j*k-6928*h*j*k+3300*i*j*k-10759*j^2*k-11607*a*k^2-148*b*k^2+12203*c*k^2+12860*d*k^2+4351*e*k^2+11456*f*k^2-812*g*k^2+4522*h*k^2+10759*i*k^2,
    h*i^2+3986*f*h*j-5005*g*h*j+10877*h^2*j-10534*a*i*j+3023*b*i*j-3535*c*i*j-3908*d*i*j+12889*e*i*j-8686*f*i*j-11302*g*i*j+791*h*i*j+2317*i^2*j-1231*a*j^2+11208*b*j^2-7885*c*j^2+12838*d*j^2-5063*e*j^2-4039*f*j^2+8757*g*j^2-3296*h*j^2-10180*i*j^2+12199*j^3+10098*f*h*k+7916*g*h*k+12053*h^2*k-4692*a*i*k+5794*b*i*k+11187*c*i*k-2800*d*i*k-2098*e*i*k+2363*f*i*k+10117*g*i*k+8716*h*i*k+833*i^2*k-2437*a*j*k-8196*b*j*k-2166*c*j*k+13249*d*j*k+13311*e*j*k+13354*f*j*k+5427*g*j*k-3697*h*j*k-12234*i*j*k+7945*j^2*k+8802*a*k^2+3154*b*k^2-2495*c*k^2+8774*d*k^2+10730*e*k^2-13164*f*k^2+13044*g*k^2+35*h*k^2-7945*i*k^2,
    g*i^2-1733*f*h*j+2153*g*h*j+8946*h^2*j-9141*a*i*j-7511*b*i*j-12962*c*i*j+2470*d*i*j+4962*e*i*j+2887*f*i*j+5307*g*i*j+4010*h*i*j+1041*i^2*j+4187*a*j^2-3757*b*j^2+10667*c*j^2+5964*d*j^2-8404*e*j^2+9723*f*j^2+13821*g*j^2+11065*h*j^2-3739*i*j^2-1235*j^3+6875*f*h*k+12884*g*h*k-3188*h^2*k-6241*a*i*k+3449*b*i*k-12780*c*i*k-8693*d*i*k-2002*e*i*k-11107*f*i*k-12399*g*i*k+2878*h*i*k-11728*i^2*k-11332*a*j*k+11055*b*j*k-8791*c*j*k-6614*d*j*k+6083*e*j*k-13507*f*j*k+694*g*j*k+5583*h*j*k-6706*i*j*k+12747*j^2*k-11412*a*k^2+8291*b*k^2+5104*c*k^2-3935*d*k^2+11263*e*k^2+12223*f*k^2+9884*g*k^2+7941*h*k^2-12747*i*k^2,
    f*i^2-6786*f*h*j+8913*g*h*j-7618*h^2*j+1661*a*i*j+6169*b*i*j+12458*c*i*j+4527*d*i*j+9240*e*i*j-9425*f*i*j+2902*g*i*j-9518*h*i*j-11444*i^2*j-2311*a*j^2-11661*b*j^2-11829*c*j^2-11259*d*j^2+7146*e*j^2-3220*f*j^2-5911*g*j^2-7139*h*j^2-6419*i*j^2-2077*j^3-10343*f*h*k-1278*g*h*k+1303*h^2*k-2262*a*i*k-18*b*i*k-720*c*i*k-6253*d*i*k+10714*e*i*k+9724*f*i*k+5510*g*i*k+6301*h*i*k+85*i^2*k+1152*a*j*k+6069*b*j*k-4995*c*j*k-31*d*j*k+7010*e*j*k-3679*f*j*k-12486*g*j*k-5408*h*j*k+4870*i*j*k+6321*j^2*k+1953*a*k^2-7849*b*k^2-6974*c*k^2-7520*d*k^2+12295*e*k^2-3133*f*k^2+11742*g*k^2-2793*h*k^2-6321*i*k^2,
    e*i^2+1959*f*h*j+11808*g*h*j-11471*h^2*j+1445*a*i*j-3674*b*i*j-12261*c*i*j+5948*d*i*j-4011*e*i*j-9167*f*i*j-13114*g*i*j+6472*h*i*j+1516*i^2*j-2766*a*j^2-1371*b*j^2-13880*c*j^2-6121*d*j^2+6112*e*j^2+10869*f*j^2+6406*g*j^2-10785*h*j^2+2038*i*j^2+9991*j^3-6278*f*h*k+11337*g*h*k+1446*h^2*k-9425*a*i*k+11264*b*i*k-7789*c*i*k+7790*d*i*k-3290*e*i*k+9542*f*i*k+7966*g*i*k+13784*h*i*k-6920*i^2*k-12817*a*j*k+12596*b*j*k+7282*c*j*k-8886*d*j*k+4834*e*j*k+13947*f*j*k+1043*g*j*k+5413*h*j*k-2042*i*j*k-12022*j^2*k+3124*a*k^2-6899*b*k^2+9701*c*k^2-11997*d*k^2-8336*e*k^2-5558*f*k^2-531*g*k^2-7949*h*k^2+12022*i*k^2,
    d*i^2-11434*f*h*j-683*g*h*j-3097*h^2*j+2568*a*i*j+6557*b*i*j+13941*c*i*j-1326*d*i*j+13933*e*i*j-11649*f*i*j-6607*g*i*j+7073*h*i*j+4510*i^2*j-8962*a*j^2+10254*b*j^2+4609*c*j^2-5569*d*j^2-5562*e*j^2-8342*f*j^2-6196*g*j^2+11845*h*j^2+6856*i*j^2-4982*j^3+6650*f*h*k+10887*g*h*k+3835*h^2*k+5728*a*i*k+7596*b*i*k+6851*c*i*k-12807*d*i*k+10090*e*i*k+3854*f*i*k+13141*g*i*k-3462*h*i*k-11932*i^2*k+3299*a*j*k-4017*b*j*k+10129*c*j*k-4118*d*j*k-3525*e*j*k-4466*f*j*k-264*g*j*k+2101*h*j*k+3068*i*j*k+7334*j^2*k-6117*a*k^2+5748*b*k^2+5272*c*k^2+6830*d*k^2-13387*e*k^2-12629*f*k^2+2975*g*k^2+1914*h*k^2-7334*i*k^2,
    c*i^2+306*f*h*j-13932*g*h*j-10600*h^2*j+12909*a*i*j+10783*b*i*j+10823*c*i*j+6315*d*i*j-10105*e*i*j+240*f*i*j+1374*g*i*j+6094*h*i*j+4604*i^2*j-7005*a*j^2+2902*b*j^2+1760*c*j^2-11553*d*j^2-2321*e*j^2-8467*f*j^2+13751*g*j^2+11173*h*j^2-6108*i*j^2-12664*j^3-1937*f*h*k-7592*g*h*k+1500*h^2*k+6892*a*i*k-11415*b*i*k-8729*c*i*k-4152*d*i*k-6306*e*i*k-13333*f*i*k+6078*g*i*k+2855*h*i*k+13641*i^2*k+6184*a*j*k+7465*b*j*k-2610*c*j*k-13035*d*j*k-13819*e*j*k+10303*f*j*k+2120*g*j*k-5422*h*j*k-12279*i*j*k+9419*j^2*k-6811*a*k^2+213*b*k^2+9390*c*k^2-3365*d*k^2-9825*e*k^2+7149*f*k^2-2111*g*k^2-2958*h*k^2-9419*i*k^2,
    b*i^2-9381*f*h*j-12981*g*h*j+1277*h^2*j-1412*a*i*j+1872*b*i*j-2552*c*i*j+12170*d*i*j-7980*e*i*j-8617*f*i*j-10095*g*i*j-10164*h*i*j-11011*i^2*j-2980*a*j^2+7138*b*j^2-5860*c*j^2+9331*d*j^2+7317*e*j^2+5474*f*j^2+8666*g*j^2-6026*h*j^2-25*i*j^2-9999*j^3-2811*f*h*k-8053*g*h*k+5432*h^2*k-13905*a*i*k-3796*b*i*k-794*c*i*k-12784*d*i*k-2489*e*i*k+9557*f*i*k+7008*g*i*k+8597*h*i*k+5756*i^2*k-13031*a*j*k-4414*b*j*k-6204*c*j*k-2148*d*j*k+8361*e*j*k-11322*f*j*k-10408*g*j*k-13779*h*j*k-4863*i*j*k+7184*j^2*k-1102*a*k^2-883*b*k^2-6172*c*k^2-6521*d*k^2+380*e*k^2-9053*f*k^2+8048*g*k^2-13039*h*k^2-7184*i*k^2,
    a*i^2+10444*f*h*j-5349*g*h*j+10973*h^2*j-570*a*i*j+6850*b*i*j+4572*c*i*j-8306*d*i*j+11780*e*i*j+9284*f*i*j+4735*g*i*j+1488*h*i*j-8967*i^2*j+13563*a*j^2+3116*b*j^2+2878*c*j^2-12287*d*j^2-7835*e*j^2-3964*f*j^2-7750*g*j^2-3865*h*j^2-6993*i*j^2+5207*j^3-6469*f*h*k+8021*g*h*k+5098*h^2*k-6824*a*i*k-9264*b*i*k+4579*c*i*k+8610*d*i*k-3234*e*i*k+12496*f*i*k+13937*g*i*k-3012*h*i*k+6157*i^2*k-3628*a*j*k+11744*b*j*k-10291*c*j*k-1913*d*j*k-8947*e*j*k-4671*f*j*k-5246*g*j*k+4961*h*j*k-12460*i*j*k+13763*j^2*k-10895*a*k^2-8256*b*k^2-12385*c*k^2+4587*d*k^2-8102*e*k^2-6811*f*k^2-4125*g*k^2+7253*h*k^2-13763*i*k^2,
    h^2*i+6525*f*h*j-12614*g*h*j-13902*h^2*j-4965*a*i*j-9250*b*i*j-7792*c*i*j+5639*d*i*j+4859*e*i*j+2757*f*i*j+11520*g*i*j+9953*h*i*j-2627*i^2*j-6131*a*j^2-7237*b*j^2+12604*c*j^2+2747*d*j^2+7080*e*j^2+140*f*j^2-12385*g*j^2-11719*h*j^2-839*i*j^2+4563*j^3+11882*f*h*k-12960*g*h*k-10829*h^2*k-4869*a*i*k+7263*b*i*k+5510*c*i*k+13539*d*i*k+5788*e*i*k-7422*f*i*k+3255*g*i*k-13165*h*i*k-11077*i^2*k-8095*a*j*k+12623*b*j*k+4514*c*j*k+7920*d*j*k-2470*e*j*k-4758*f*j*k-8282*g*j*k-1846*h*j*k+6465*i*j*k-9821*j^2*k-8475*a*k^2-4283*b*k^2+5087*c*k^2-2807*d*k^2-13137*e*k^2+7892*f*k^2+13762*g*k^2-11028*h*k^2+9821*i*k^2,
    g*h*i-3398*f*h*j-340*g*h*j+2522*h^2*j-162*a*i*j+2418*b*i*j-4729*c*i*j-5800*d*i*j-13619*e*i*j+9162*f*i*j+11331*g*i*j+8589*h*i*j-8101*i^2*j+8629*a*j^2+6015*b*j^2-4878*c*j^2+13679*d*j^2+12410*e*j^2+10655*f*j^2-13355*g*j^2-6083*h*j^2+3665*i*j^2+1765*j^3-2109*f*h*k+11693*g*h*k-5786*h^2*k+688*a*i*k-8874*b*i*k-5429*c*i*k+10232*d*i*k+7387*e*i*k-982*f*i*k-10816*g*i*k-2016*h*i*k+4817*i^2*k-9775*a*j*k+9187*b*j*k-7892*c*j*k+2425*d*j*k-3593*e*j*k-4979*f*j*k+2045*g*j*k+3030*h*j*k-4613*i*j*k+10993*j^2*k+6920*a*k^2+997*b*k^2-1034*c*k^2-3725*d*k^2-1554*e*k^2-13746*f*k^2-11512*g*k^2+2848*h*k^2-10993*i*k^2,
    f*h*i+12293*f*h*j+3328*g*h*j-11357*h^2*j-8177*a*i*j-7209*b*i*j+3219*c*i*j-2945*d*i*j+7786*e*i*j+5840*f*i*j-2238*g*i*j-8212*h*i*j-10841*i^2*j-13266*a*j^2+9578*b*j^2-3743*c*j^2+9897*d*j^2-2829*e*j^2+12287*f*j^2+5907*g*j^2-2637*h*j^2-5367*i*j^2+13842*j^3+7605*f*h*k+2005*g*h*k+4086*h^2*k+6604*a*i*k+7922*b*i*k+12581*c*i*k+11650*d*i*k+1717*e*i*k-7842*f*i*k+1251*g*i*k+6548*h*i*k-11626*i^2*k-9556*a*j*k-8929*b*j*k+94*c*j*k-10750*d*j*k-10434*e*j*k-7950*f*j*k-7472*g*j*k-4832*h*j*k+8928*i*j*k+10747*j^2*k+3035*a*k^2-13819*b*k^2-4911*c*k^2-10322*d*k^2+4918*e*k^2-13499*f*k^2-6076*g*k^2+5131*h*k^2-10747*i*k^2,
    h^3+6852*f*h*j+12078*g*h*j-3214*h^2*j-9013*a*i*j+6269*b*i*j+4611*c*i*j+3804*d*i*j+13344*e*i*j-3164*f*i*j+9351*g*i*j+6326*h*i*j+8222*i^2*j+4410*a*j^2+9795*b*j^2-9249*c*j^2-11407*d*j^2-85*e*j^2-3790*f*j^2-13667*g*j^2+7405*h*j^2+9765*i*j^2+8538*j^3+7247*f*h*k+2342*g*h*k-1038*h^2*k+7520*a*i*k+12872*b*i*k+7140*c*i*k+1463*d*i*k-9984*e*i*k+2969*f*i*k-12825*g*i*k-5183*h*i*k+1129*i^2*k+3829*a*j*k+3072*b*j*k-4053*c*j*k-11993*d*j*k-10598*e*j*k-2794*f*j*k+2679*g*j*k+12134*h*j*k+6933*i*j*k+3678*j^2*k-4767*a*k^2-6200*b*k^2+5901*c*k^2+2940*d*k^2-3903*e*k^2-13123*f*k^2+4873*g*k^2+12430*h*k^2-3678*i*k^2,
    g*h^2+7385*f*h*j+227*g*h*j-4000*h^2*j+13868*a*i*j+5702*b*i*j-4846*c*i*j-6478*d*i*j-1224*e*i*j+7899*f*i*j-1048*g*i*j+4878*h*i*j+4384*i^2*j-7228*a*j^2+3928*b*j^2-9634*c*j^2+2157*d*j^2-8058*e*j^2+7230*f*j^2+13158*g*j^2-4197*h*j^2-11355*i*j^2+8023*j^3+906*f*h*k+2699*g*h*k-10257*h^2*k-12311*a*i*k-4730*b*i*k+3116*c*i*k-4241*d*i*k-2251*e*i*k+8434*f*i*k+2763*g*i*k+13671*h*i*k+5898*i^2*k+8923*a*j*k+9586*b*j*k-9073*c*j*k-10317*d*j*k+9399*e*j*k-11119*f*j*k-11728*g*j*k+5368*h*j*k-6817*i*j*k+7819*j^2*k+3409*a*k^2+4996*b*k^2+11594*c*k^2+5187*d*k^2+577*e*k^2-2130*f*k^2+89*g*k^2-1206*h*k^2-7819*i*k^2,
    f*h^2-13069*f*h*j+8885*g*h*j+4681*h^2*j+3056*a*i*j+10800*b*i*j-6139*c*i*j-12227*d*i*j-8097*e*i*j+13107*f*i*j+9996*g*i*j-5023*h*i*j-8596*i^2*j+3593*a*j^2+4861*b*j^2-3157*c*j^2-5605*d*j^2+3867*e*j^2-623*f*j^2+13896*g*j^2+11348*h*j^2-605*i*j^2+4932*j^3-10973*f*h*k+1149*g*h*k+12797*h^2*k+12021*a*i*k+12146*b*i*k-5724*c*i*k-11471*d*i*k-12564*e*i*k-1385*f*i*k-11123*g*i*k-13327*h*i*k-2358*i^2*k-7484*a*j*k-9199*b*j*k+12087*c*j*k-9338*d*j*k-5559*e*j*k-742*f*j*k-12469*g*j*k+13647*h*j*k+12242*i*j*k-6532*j^2*k+2406*a*k^2-1746*b*k^2+7016*c*k^2-8259*d*k^2-9805*e*k^2-4857*f*k^2-10684*g*k^2+10727*h*k^2+6532*i*k^2
    );
