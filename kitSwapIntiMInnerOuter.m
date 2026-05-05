function newIntiM = kitSwapIntiMInnerOuter(oldIntiM)
% KITSWAPINTIMINNEROUTER rearranges an intiM so that the 'inner'
% measurements become the 'outer' measurements, and vice versa
%
% C.C.Conway 2023

newIntiM = oldIntiM;
newIntiM.intensity.mean.inner = oldIntiM.intensity.mean.outer;
newIntiM.intensity.mean.outer = oldIntiM.intensity.mean.inner;

newIntiM.intensity.max.inner = oldIntiM.intensity.max.outer;
newIntiM.intensity.max.outer = oldIntiM.intensity.max.inner;

newIntiM.intensity.bg.inner = oldIntiM.intensity.bg.outer;
newIntiM.intensity.bg.outer = oldIntiM.intensity.bg.inner;

end