
-- old change of rings --
--pushFwdChainComplex = method()
--pushFwdChainComplex(ChainComplex,RingMap) := (C,f) -> (
--    D := new ChainComplex;
--    D.ring = source f;
--    for i from min C to max C do
--    D.dd _i = pushFwd(C.dd_i,f);    
--    D
--    )

-- changeOfRingsTor = method()
-- changeOfRingsTor(Module,Module,RingMap) := (M,N,f) -> (
--    -- f : R --> S, N an S module, M an R module
--    F := complete res N;
--    FR := pushFwdChainComplex(F,f);
--    G := complete res M;
--    spectralSequence((G) ** (filteredComplex FR) )
--    )




-- change of rings --

--Compute change of rings for Tor

changeOfRingsTor = method()
changeOfRingsTor(Module,Module,RingMap) := (M,N,f) -> (
    -- f : R --> S a finite ring map, N an S module, M an R module
    F := complete res N;
    pushFwdF := pushFwd(f,F);
    G := complete res M;
    E := spectralSequence(filteredComplex(G) ** pushFwdF);
    EE := spectralSequence(G ** (filteredComplex pushFwdF));
    (E,EE) 
)
doc ///
     Key
       changeOfRingsTor
     Headline
         compute the change of rings spectral sequence
     Usage
          E = changeOfRingsTor(M,N,f)
     Inputs
     	 M:Module
	 N:Module
	 f:RingMap
     Outputs
     	  E:Sequence
     Description
          Text
	       This method computes the change of ring spectral sequence for cetain kinds of ring maps.
	  Example
	       k=QQ;
	       R=k[a,b,c];
	       S=k[s,t];
	       f = map(S,R,{s^2,s*t,t^2});
--	       kappa = coker vars S;
--	       kkappa = coker vars R;
--	       (E,EE) = changeOfRingsTor(kkappa,kappa,f);
--	       e = prune E
--	       ee = prune EE
--	       e^0
--	       e^1
--	       e^2
--	       e^infinity
--	       ee^0
--	       ee^1
--	       ee^2
--	       (ee^2).dd
--	       ee^3
--	       ee^infinity   	      
///