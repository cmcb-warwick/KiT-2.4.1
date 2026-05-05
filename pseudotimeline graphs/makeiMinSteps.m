function makeiMinSteps(nBins, IntiM, iM, normalised, markerdist, intID, expnum, expcond)
% Before you start: navigate to BEDCA_v1/BEDCA/ExptData in your MATLAB
% folder list (generally on left of screen). The structure created here
% will be saved in this folder.
% 
% Function to split an iM into 'nBins' bins based on Mad2 intensity.
% Assumes that the marker you are binning by is 'outer' marker in intiM.
% This will create and save a parent structure with 'nBins' substructures, 
% each named based on which bin they are. The main structure will be named
% as follows, with values in curly brackets below (user-defined unless
% otherwise stated):
% iM_exp{expnum}_dist{markerdist}_{normalised}int{intID}_ss{stepsize - not user-defined}_{expcond}.mat
%
% Suggested values:
% nBins: 25 or other integer
% IntiM: the intiM corresponding to your iM
% iM: the distance matrix between your two markers
% normalised: string. Put 'norm' or 'normalised' to normalise 'outer' intensity to 
%   'inner' intensity of your intiM. Write '' if you have combined experiments
%   using kitNormCombineiMsIntiMs (as this will already be normalised). If
%   you write 'norm' or 'normalised', the resulting variable will have 'N'
%   appended before the intID.
% markerdist: string. The distance between two markers. Shorten if possible
%   e.g. for CenpC to Ndc80-N, you could write 'CCNN'.
% intID: string. What intensity you've binned by. e.g. 'Mad2'. If you've
%   combined and normalised using kitNormCombineiMsIntiMs, you could write
%   'NMad2' here.
% expnum: string. Either single number or combination of multiple
%   experiment numbers (e.g. '010512' for a combination of expts 1, 5 and
%   12). You may also include a date here if you are doing multiple repeats
%   of the same experiment, e.g. '010512260130' for a run of expts 1, 5, and
%   12 performed on 30th Jan 2026. (only if desired!)
% expcond: string. '' if no treatment has been used, 'DMSO', 'ZM', 'Mon' etc
%   as befits your experiment. As you will have lots of bins, this portion
%   is not used in BEDCA, but is good for your own organisation of data.
%   Could also put date here if you like.
%   
% Copyright (c) 2021 T. E. Germanova, edited 2022 C. C. Conway

Total = size(find(~isnan(IntiM.intensity.mean.inner(:))),1);
stepsize = floor(Total/nBins);
bin = 1;
n = stepsize;
k = n;
if strcmp(normalised,'norm') || strcmp(normalised,'normalised')
    nStr = 'N';
    nv = IntiM.intensity.mean.outer./IntiM.intensity.mean.inner;
else
    nStr = '';
    nv = IntiM.intensity.mean.outer;
end

while n < Total+1
    iM_copy = iM;
    Subset = mink(maxk(nv(:), n), k); %20250430
    
    positionLow = nv<min(Subset);
    positionLowSisSep = max(positionLow, [], 2);
    positionHigh = nv>max(Subset);
    positionHighSisSep = max(positionHigh, [], 2);
    
    iM_copy.microscope.raw.delta.x.all(positionLow) = NaN;
    iM_copy.microscope.raw.delta.y.all(positionLow) = NaN;
    iM_copy.microscope.raw.delta.z.all(positionLow) = NaN;
    if isfield (iM_copy.microscope, 'sisSep')
        iM_copy.microscope.sisSep.threeD(positionLowSisSep) = NaN;
        iM_copy.microscope.raw.swivel.threeD.all(positionLow) = NaN;
    end
    
    iM_copy.microscope.raw.delta.x.all(positionHigh) = NaN;
    iM_copy.microscope.raw.delta.y.all(positionHigh) = NaN;
    iM_copy.microscope.raw.delta.z.all(positionHigh) = NaN;
    if isfield (iM_copy.microscope, 'sisSep')
        iM_copy.microscope.sisSep.threeD(positionHighSisSep) = NaN;
        iM_copy.microscope.raw.swivel.threeD.all(positionHigh) = NaN;
    end
    
    fieldname = sprintf('iM_exp%s_dist%s_%sint%s_ss%d_%sBin%02d', expnum, markerdist, nStr, intID, stepsize, expcond, bin);
    workingStruct = struct(fieldname, iM_copy);
    f = fieldnames(workingStruct);
    for i = 1:length(f)
        masterStruct.(f{i}) = workingStruct.(f{i});
    end
    n = n+k;
    bin = bin+1;

end
fieldname = sprintf('iM_exp%s_dist%s_%sint%s_ss%d_%s.mat', expnum, markerdist, nStr, intID, stepsize, expcond);
save(fieldname, '-struct', 'masterStruct');

end
%%
