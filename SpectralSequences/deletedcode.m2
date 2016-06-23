
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
