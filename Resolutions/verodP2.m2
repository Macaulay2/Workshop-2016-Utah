restart
deg = 5
P2 = ZZ/101[a,b,c]
B = basis(deg, P2)
Pn = ZZ/101[x_1..x_(numColumns B)]
phi = map(P2, Pn, B)
I = kernel phi;
betti I
gbTrace=2
elapsedTime C = res(I, FastNonminimal=>true)
betti(C, Minimize=>true)
