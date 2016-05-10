restart 
loadPackage "PHCpack"
R = CC[x,y];
f =  {x^3*y^5 + y^2 + x^2*y, x*y + x^2 - 1};
I = ideal f
(mv,q,qsols) = mixedVolume(f,StartSystem=>true)
fsols = trackPaths(f,q,qsols, interactive=>true, saveSettingsPath=>"settings.phc")
fsols2 = trackPaths(f,q,qsols, loadSettingsPath=>"settings.phc")
fsols == fsols2
fIntermediateSols = trackPaths(f,q,qsols, loadSettingsPath=>"settings.phc", intermediateSolutions=>true);
netList fIntermediateSols
fsols3 = flatten for l in fIntermediateSols list ( last l )
fsols3 == fsols