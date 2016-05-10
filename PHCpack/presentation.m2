restart 
loadPackage "PHCpack"
R = CC[x,y];
f =  {x^3*y^5 + y^2 + x^2*y, x*y + x^2 - 1};
I = ideal f
--(mv,smv,q,qsols) = mixedVolume(f,StableMixedVolume=>true,StartSystem=>true)
(mv,q,qsols) = mixedVolume(f,StartSystem=>true)
fsols = trackPaths(f,q,qsols, interactive=>true, saveSettingsPath=>"settings.phc")
fsols2 = trackPaths(f,q,qsols, loadSettingsPath=>"settings.phc")
fsols == fsols2
fIntermediateSols = trackPaths(f,q,qsols, loadSettingsPath=>"settings.phc", intermediateSolutions=>true)
fsols3 = flatten for l in fIntermediateSols list ( last l )
fsols3 == fsols
--m = mixedVolume(f)
--(mv,sv) = mixedVolume(f,StableMixedVolume => true)
--mv = mixedVolume(f,interactive=>true)
--(mv,smv,q,qsols) = mixedVolume(f,interactive=>true)
--mixedVolume(f,interactive=>true)
--fsols = trackPaths(f,q,qsols)