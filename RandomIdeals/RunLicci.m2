needsPackage "RandomIdeal"
needsPackage "ResidualIntersections"
runLicci := (R,d,k) -> (
    L = apply(100,i->randomShelling(#R,d,k));
    LI = apply(L,P -> idealFromShelling(R,P) );
    Licci = apply(LI,isLicci);
    (L,LI,Licci)
    )

getLicciStats := (L,LI,Licci) -> (
    #select(Licci, i -> i)
    )

end--

R=QQ[x_1..x_5]
(L,LI,Licci) = runLicci(R,1,5)
getLicciStats(L,LI,Licci)
