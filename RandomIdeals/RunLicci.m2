needsPackage "RandomIdeal"
needsPackage "ResidualIntersections"
runLicci := (R,d,k) -> (
    L = apply(100,i->randomShelling(numgens R,d,k));
    LI = apply(L,P -> idealFromShelling(R,P) );
    Licci = apply(LI,isLicci);
    (L,LI,Licci)
    )

getLicciStats := (L,LI,Licci) -> (
    #select(Licci, i -> i)
    )

end--

R=QQ[x_1..x_6]
{*
scan((1..10),i-> ((L,LI,Licci) = runLicci(R,1,i);
        print getLicciStats(L,LI,Licci)))
*}

needsPackage "Visualize"
R=QQ[x_1..x_6];
(L,LI,Licci) = runLicci(R,1,6);
print getLicciStats(L,LI,Licci)
scan(100,i -> if Licci_i then visualize graph L_i)
