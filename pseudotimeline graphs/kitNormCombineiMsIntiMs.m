function [iMcomb, intiMcomb] = kitNormCombineiMsIntiMs(iMs, intiMs, ctrlintiMs, varargin)
% KITNORMCOMBINEIMSINTIMS is a function to combine iMs and intiMs from two
% or more experimental repeats, where the resultant iMs and intiMs will be
% used to rerun BEDCA. This function will order the intensity values for
% each experiment from highest to lowest, split into bins, and use the
% median of the highest intensity bin for normalisation, where the
% intensity values will be divided by this value per experiment. In the
% case that you have a perturbation experiment, provided control intiMs
% will provide the normalisation values per experiment, and your intiMs
% from the perturbation will be normalised to the control (i.e. DMSO).
% The resultant iM will be a simple concatenation, whilst the resultant
% intiM will be a concatenation with values from each experiment normalised
% to 1 in the intensity.mean.outer portion of the structure. Do not use
% 'normalised' as an option in makeiMinSteps function (even though the data
% may be normalised to CenpC if you have used 'outerDivInner' as
% normalisation marker).
%
%
%   KITNORMCOMBINEIMSINTIMS(iMs, intiMs, ctrlintiMs)
%       - iMs should be organised as {iM1, iM2, ...}
%       - intiMs should be organised as {intiM1, intiM2, ...}
%       - ctrlintiMs are optional, these should only be used if you have an
%         experimental treatment (e.g. monastrol) whose iMs you wish to
%         concatenate. In this case, the DMSO/control treatment intiMs
%         should be organised as {ctrlintiM1, ctrlintiM2, ...}. If you do
%         not need to use this feature, please use [] in its place.
%
%   Options are available, defaults in {}:
%
%   nBins: {25} or other integer. The number of bins to split each intiM
%       structure into (the highest bin will be used for normalisation).
%
%   normalisationMarker: 'outer', {'outerDivInner'}, 'inner', or 'multi'.
%       Which position in the intiM the intensity marker is located in. if
%       'multi', please also provide option 'normalisationOrder'.
%
%   normalisationOrder: {nan} or cell of strings in form {'outer', 'inner',
%       'outerDivInner', ...} in the same order that your intiMs have been
%       provided. Ensure that you have the same number of values as iMs.
%       Only use if normalisationMarker value is 'multi'.
%
%
% Copyright (c) 2024 C. C. Conway

opts.nBins = 25;
opts.normalisationMarker = 'outerDivInner';
opts.normalisationOrder = nan;
opts = processOptions(opts,varargin{:});

nBins = opts.nBins;
nExpts = length(iMs);

if ~strcmp('multi', opts.normalisationMarker)
    opts.normalisationOrder = repmat({opts.normalisationMarker}, 1, nExpts);
end

%indicate if need to normalise to non-perturbed intiM
if isempty(ctrlintiMs)
    controlNorm = 0;
else
    controlNorm = 1;
end

iMcomb = struct;
intiMcomb = struct;

for iExpt = 1:nExpts
    wkiM = iMs{iExpt};
    wkintiM = intiMs{iExpt};
    wkintiM.intensity.mean.inner = rmLowCenpC(wkintiM); %20250430
    if controlNorm
        wkDivIntiM = ctrlintiMs{iExpt};
    else
        wkDivIntiM = intiMs{iExpt};
    end
    wkDivIntiM.intensity.mean.inner = rmLowCenpC(wkDivIntiM); %20250430
    %bin data on intiM to divide by (either untreated partner or original)
    binnedExpt = makeBinsintiMs(wkDivIntiM, nBins, opts.normalisationOrder{iExpt});
    maxMedianBin = max(median(binnedExpt));
    %divide relevant intensity by max median value of control intiM.
    if strcmp(opts.normalisationOrder{iExpt}, 'outer')
        wkIntensity = wkintiM.intensity.mean.outer/maxMedianBin;
    elseif strcmp(opts.normalisationOrder{iExpt}, 'outerDivInner')
        wkIntensity = wkintiM.intensity.mean.outer ./ wkintiM.intensity.mean.inner;
        wkIntensity = wkIntensity/maxMedianBin;
    elseif strcmp(opts.normalisationOrder{iExpt}, 'inner')
        wkIntensity = wkintiM.intensity.mean.inner/maxMedianBin;
    end
    
    %iMcomb structure
    increaseLabel = (iExpt-1)*100000;
    %if new iMcomb then create fields. If not, update.
    %kitVersion - intiMcomb is same
    if ~isfield(iMcomb, 'kitVersion')
        iMcomb.kitVersion = {wkiM.kitVersion};
    else
        iMcomb.kitVersion = vertcat(iMcomb.kitVersion, {wkiM.kitVersion});
    end
    intiMcomb.kitVersion = iMcomb.kitVersion; %replicating

    %label - intiMcomb is same
    if ~isfield(iMcomb, 'label')
        iMcomb.label = wkiM.label;
    else
        labNums = str2num(wkiM.label);
        labNums = labNums + increaseLabel;
        newLabels = num2str(labNums, '%07d');
        iMcomb.label = vertcat(iMcomb.label, newLabels);
    end
    intiMcomb.label = iMcomb.label; %replicating
    
    %direction
    if ~isfield(iMcomb, 'direction')
        iMcomb.direction = wkiM.direction;
    else
        iMcomb.direction.P  = vertcat(iMcomb.direction.P, wkiM.direction.P);
        iMcomb.direction.AP = vertcat(iMcomb.direction.AP, wkiM.direction.AP);
        iMcomb.direction.N  = vertcat(iMcomb.direction.N, wkiM.direction.N);
        iMcomb.direction.S  = vertcat(iMcomb.direction.S, wkiM.direction.S);
    end

    %microscope
    if ~isfield(iMcomb, 'microscope')
        iMcomb.microscope = wkiM.microscope;
    else
        %% coords x y z
        iMcomb.microscope.coords.x = vertcat(iMcomb.microscope.coords.x, wkiM.microscope.coords.x);
        iMcomb.microscope.coords.y = vertcat(iMcomb.microscope.coords.y, wkiM.microscope.coords.y);
        iMcomb.microscope.coords.z = vertcat(iMcomb.microscope.coords.z, wkiM.microscope.coords.z);
        %% sisSep x y z twoD threeD
        iMcomb.microscope.sisSep.x = vertcat(iMcomb.microscope.sisSep.x, wkiM.microscope.sisSep.x);
        iMcomb.microscope.sisSep.y = vertcat(iMcomb.microscope.sisSep.y, wkiM.microscope.sisSep.y);
        iMcomb.microscope.sisSep.z = vertcat(iMcomb.microscope.sisSep.z, wkiM.microscope.sisSep.z);
        iMcomb.microscope.sisSep.twoD = vertcat(iMcomb.microscope.sisSep.twoD, wkiM.microscope.sisSep.twoD);
        iMcomb.microscope.sisSep.threeD = vertcat(iMcomb.microscope.sisSep.threeD, wkiM.microscope.sisSep.threeD);
        %% raw.delta
        % raw.delta.x all P AP N S
        iMcomb.microscope.raw.delta.x.all = vertcat(iMcomb.microscope.raw.delta.x.all, wkiM.microscope.raw.delta.x.all);
        iMcomb.microscope.raw.delta.x.P = vertcat(iMcomb.microscope.raw.delta.x.P, wkiM.microscope.raw.delta.x.P);
        iMcomb.microscope.raw.delta.x.AP = vertcat(iMcomb.microscope.raw.delta.x.AP, wkiM.microscope.raw.delta.x.AP);
        iMcomb.microscope.raw.delta.x.N = vertcat(iMcomb.microscope.raw.delta.x.N, wkiM.microscope.raw.delta.x.N);
        iMcomb.microscope.raw.delta.x.S = vertcat(iMcomb.microscope.raw.delta.x.S, wkiM.microscope.raw.delta.x.S);
        % raw.delta.y all P AP N S
        iMcomb.microscope.raw.delta.y.all = vertcat(iMcomb.microscope.raw.delta.y.all, wkiM.microscope.raw.delta.y.all);
        iMcomb.microscope.raw.delta.y.P = vertcat(iMcomb.microscope.raw.delta.y.P, wkiM.microscope.raw.delta.y.P);
        iMcomb.microscope.raw.delta.y.AP = vertcat(iMcomb.microscope.raw.delta.y.AP, wkiM.microscope.raw.delta.y.AP);
        iMcomb.microscope.raw.delta.y.N = vertcat(iMcomb.microscope.raw.delta.y.N, wkiM.microscope.raw.delta.y.N);
        iMcomb.microscope.raw.delta.y.S = vertcat(iMcomb.microscope.raw.delta.y.S, wkiM.microscope.raw.delta.y.S);
        % raw.delta.z all P AP N S
        iMcomb.microscope.raw.delta.z.all = vertcat(iMcomb.microscope.raw.delta.z.all, wkiM.microscope.raw.delta.z.all);
        iMcomb.microscope.raw.delta.z.P = vertcat(iMcomb.microscope.raw.delta.z.P, wkiM.microscope.raw.delta.z.P);
        iMcomb.microscope.raw.delta.z.AP = vertcat(iMcomb.microscope.raw.delta.z.AP, wkiM.microscope.raw.delta.z.AP);
        iMcomb.microscope.raw.delta.z.N = vertcat(iMcomb.microscope.raw.delta.z.N, wkiM.microscope.raw.delta.z.N);
        iMcomb.microscope.raw.delta.z.S = vertcat(iMcomb.microscope.raw.delta.z.S, wkiM.microscope.raw.delta.z.S);
        % raw.delta.oneD
        iMcomb.microscope.raw.delta.oneD = vertcat(iMcomb.microscope.raw.delta.oneD, wkiM.microscope.raw.delta.oneD);
        % raw.delta.twoD all P AP N S
        iMcomb.microscope.raw.delta.twoD.all = vertcat(iMcomb.microscope.raw.delta.twoD.all, wkiM.microscope.raw.delta.twoD.all);
        iMcomb.microscope.raw.delta.twoD.P = vertcat(iMcomb.microscope.raw.delta.twoD.P, wkiM.microscope.raw.delta.twoD.P);
        iMcomb.microscope.raw.delta.twoD.AP = vertcat(iMcomb.microscope.raw.delta.twoD.AP, wkiM.microscope.raw.delta.twoD.AP);
        iMcomb.microscope.raw.delta.twoD.N = vertcat(iMcomb.microscope.raw.delta.twoD.N, wkiM.microscope.raw.delta.twoD.N);
        iMcomb.microscope.raw.delta.twoD.S = vertcat(iMcomb.microscope.raw.delta.twoD.S, wkiM.microscope.raw.delta.twoD.S);
        % raw.delta.threeD all P AP N S
        iMcomb.microscope.raw.delta.threeD.all = vertcat(iMcomb.microscope.raw.delta.threeD.all, wkiM.microscope.raw.delta.threeD.all);
        iMcomb.microscope.raw.delta.threeD.P = vertcat(iMcomb.microscope.raw.delta.threeD.P, wkiM.microscope.raw.delta.threeD.P);
        iMcomb.microscope.raw.delta.threeD.AP = vertcat(iMcomb.microscope.raw.delta.threeD.AP, wkiM.microscope.raw.delta.threeD.AP);
        iMcomb.microscope.raw.delta.threeD.N = vertcat(iMcomb.microscope.raw.delta.threeD.N, wkiM.microscope.raw.delta.threeD.N);
        iMcomb.microscope.raw.delta.threeD.S = vertcat(iMcomb.microscope.raw.delta.threeD.S, wkiM.microscope.raw.delta.threeD.S);
        %% raw.swivel
        % raw.swivel.y all P AP N S
        iMcomb.microscope.raw.swivel.y.all = vertcat(iMcomb.microscope.raw.swivel.y.all, wkiM.microscope.raw.swivel.y.all);
        iMcomb.microscope.raw.swivel.y.P = vertcat(iMcomb.microscope.raw.swivel.y.P, wkiM.microscope.raw.swivel.y.P);
        iMcomb.microscope.raw.swivel.y.AP = vertcat(iMcomb.microscope.raw.swivel.y.AP, wkiM.microscope.raw.swivel.y.AP);
        iMcomb.microscope.raw.swivel.y.N = vertcat(iMcomb.microscope.raw.swivel.y.N, wkiM.microscope.raw.swivel.y.N);
        iMcomb.microscope.raw.swivel.y.S = vertcat(iMcomb.microscope.raw.swivel.y.S, wkiM.microscope.raw.swivel.y.S);
        % raw.swivel.z all P AP N S
        iMcomb.microscope.raw.swivel.z.all = vertcat(iMcomb.microscope.raw.swivel.z.all, wkiM.microscope.raw.swivel.z.all);
        iMcomb.microscope.raw.swivel.z.P = vertcat(iMcomb.microscope.raw.swivel.z.P, wkiM.microscope.raw.swivel.z.P);
        iMcomb.microscope.raw.swivel.z.AP = vertcat(iMcomb.microscope.raw.swivel.z.AP, wkiM.microscope.raw.swivel.z.AP);
        iMcomb.microscope.raw.swivel.z.N = vertcat(iMcomb.microscope.raw.swivel.z.N, wkiM.microscope.raw.swivel.z.N);
        iMcomb.microscope.raw.swivel.z.S = vertcat(iMcomb.microscope.raw.swivel.z.S, wkiM.microscope.raw.swivel.z.S);
        % raw.swivel.threeD all P AP N S
        iMcomb.microscope.raw.swivel.threeD.all = vertcat(iMcomb.microscope.raw.swivel.threeD.all, wkiM.microscope.raw.swivel.threeD.all);
        iMcomb.microscope.raw.swivel.threeD.P = vertcat(iMcomb.microscope.raw.swivel.threeD.P, wkiM.microscope.raw.swivel.threeD.P);
        iMcomb.microscope.raw.swivel.threeD.AP = vertcat(iMcomb.microscope.raw.swivel.threeD.AP, wkiM.microscope.raw.swivel.threeD.AP);
        iMcomb.microscope.raw.swivel.threeD.N = vertcat(iMcomb.microscope.raw.swivel.threeD.N, wkiM.microscope.raw.swivel.threeD.N);
        iMcomb.microscope.raw.swivel.threeD.S = vertcat(iMcomb.microscope.raw.swivel.threeD.S, wkiM.microscope.raw.swivel.threeD.S);
        % raw.swivel.kMT
        iMcomb.microscope.raw.swivel.kMT = vertcat(iMcomb.microscope.raw.swivel.kMT, wkiM.microscope.raw.swivel.kMT);
        %% depthFilter.delta
        % depthFilter.delta.x all P AP N S
        iMcomb.microscope.depthFilter.delta.x.all = vertcat(iMcomb.microscope.depthFilter.delta.x.all, wkiM.microscope.depthFilter.delta.x.all);
        iMcomb.microscope.depthFilter.delta.x.P = vertcat(iMcomb.microscope.depthFilter.delta.x.P, wkiM.microscope.depthFilter.delta.x.P);
        iMcomb.microscope.depthFilter.delta.x.AP = vertcat(iMcomb.microscope.depthFilter.delta.x.AP, wkiM.microscope.depthFilter.delta.x.AP);
        iMcomb.microscope.depthFilter.delta.x.N = vertcat(iMcomb.microscope.depthFilter.delta.x.N, wkiM.microscope.depthFilter.delta.x.N);
        iMcomb.microscope.depthFilter.delta.x.S = vertcat(iMcomb.microscope.depthFilter.delta.x.S, wkiM.microscope.depthFilter.delta.x.S);
        % depthFilter.delta.y all P AP N S
        iMcomb.microscope.depthFilter.delta.y.all = vertcat(iMcomb.microscope.depthFilter.delta.y.all, wkiM.microscope.depthFilter.delta.y.all);
        iMcomb.microscope.depthFilter.delta.y.P = vertcat(iMcomb.microscope.depthFilter.delta.y.P, wkiM.microscope.depthFilter.delta.y.P);
        iMcomb.microscope.depthFilter.delta.y.AP = vertcat(iMcomb.microscope.depthFilter.delta.y.AP, wkiM.microscope.depthFilter.delta.y.AP);
        iMcomb.microscope.depthFilter.delta.y.N = vertcat(iMcomb.microscope.depthFilter.delta.y.N, wkiM.microscope.depthFilter.delta.y.N);
        iMcomb.microscope.depthFilter.delta.y.S = vertcat(iMcomb.microscope.depthFilter.delta.y.S, wkiM.microscope.depthFilter.delta.y.S);
        % depthFilter.delta.z all P AP N S
        iMcomb.microscope.depthFilter.delta.z.all = vertcat(iMcomb.microscope.depthFilter.delta.z.all, wkiM.microscope.depthFilter.delta.z.all);
        iMcomb.microscope.depthFilter.delta.z.P = vertcat(iMcomb.microscope.depthFilter.delta.z.P, wkiM.microscope.depthFilter.delta.z.P);
        iMcomb.microscope.depthFilter.delta.z.AP = vertcat(iMcomb.microscope.depthFilter.delta.z.AP, wkiM.microscope.depthFilter.delta.z.AP);
        iMcomb.microscope.depthFilter.delta.z.N = vertcat(iMcomb.microscope.depthFilter.delta.z.N, wkiM.microscope.depthFilter.delta.z.N);
        iMcomb.microscope.depthFilter.delta.z.S = vertcat(iMcomb.microscope.depthFilter.delta.z.S, wkiM.microscope.depthFilter.delta.z.S);
        % depthFilter.delta.oneD
        iMcomb.microscope.depthFilter.delta.oneD = vertcat(iMcomb.microscope.depthFilter.delta.oneD, wkiM.microscope.depthFilter.delta.oneD);
        % depthFilter.delta.twoD all P AP N S
        iMcomb.microscope.depthFilter.delta.twoD.all = vertcat(iMcomb.microscope.depthFilter.delta.twoD.all, wkiM.microscope.depthFilter.delta.twoD.all);
        iMcomb.microscope.depthFilter.delta.twoD.P = vertcat(iMcomb.microscope.depthFilter.delta.twoD.P, wkiM.microscope.depthFilter.delta.twoD.P);
        iMcomb.microscope.depthFilter.delta.twoD.AP = vertcat(iMcomb.microscope.depthFilter.delta.twoD.AP, wkiM.microscope.depthFilter.delta.twoD.AP);
        iMcomb.microscope.depthFilter.delta.twoD.N = vertcat(iMcomb.microscope.depthFilter.delta.twoD.N, wkiM.microscope.depthFilter.delta.twoD.N);
        iMcomb.microscope.depthFilter.delta.twoD.S = vertcat(iMcomb.microscope.depthFilter.delta.twoD.S, wkiM.microscope.depthFilter.delta.twoD.S);
        % depthFilter.delta.threeD all P AP N S
        iMcomb.microscope.depthFilter.delta.threeD.all = vertcat(iMcomb.microscope.depthFilter.delta.threeD.all, wkiM.microscope.depthFilter.delta.threeD.all);
        iMcomb.microscope.depthFilter.delta.threeD.P = vertcat(iMcomb.microscope.depthFilter.delta.threeD.P, wkiM.microscope.depthFilter.delta.threeD.P);
        iMcomb.microscope.depthFilter.delta.threeD.AP = vertcat(iMcomb.microscope.depthFilter.delta.threeD.AP, wkiM.microscope.depthFilter.delta.threeD.AP);
        iMcomb.microscope.depthFilter.delta.threeD.N = vertcat(iMcomb.microscope.depthFilter.delta.threeD.N, wkiM.microscope.depthFilter.delta.threeD.N);
        iMcomb.microscope.depthFilter.delta.threeD.S = vertcat(iMcomb.microscope.depthFilter.delta.threeD.S, wkiM.microscope.depthFilter.delta.threeD.S);
        %% depthFilter.swivel
        % depthFilter.swivel.y all P AP N S
        iMcomb.microscope.depthFilter.swivel.y.all = vertcat(iMcomb.microscope.depthFilter.swivel.y.all, wkiM.microscope.depthFilter.swivel.y.all);
        iMcomb.microscope.depthFilter.swivel.y.P = vertcat(iMcomb.microscope.depthFilter.swivel.y.P, wkiM.microscope.depthFilter.swivel.y.P);
        iMcomb.microscope.depthFilter.swivel.y.AP = vertcat(iMcomb.microscope.depthFilter.swivel.y.AP, wkiM.microscope.depthFilter.swivel.y.AP);
        iMcomb.microscope.depthFilter.swivel.y.N = vertcat(iMcomb.microscope.depthFilter.swivel.y.N, wkiM.microscope.depthFilter.swivel.y.N);
        iMcomb.microscope.depthFilter.swivel.y.S = vertcat(iMcomb.microscope.depthFilter.swivel.y.S, wkiM.microscope.depthFilter.swivel.y.S);
        % depthFilter.swivel.z all P AP N S
        iMcomb.microscope.depthFilter.swivel.z.all = vertcat(iMcomb.microscope.depthFilter.swivel.z.all, wkiM.microscope.depthFilter.swivel.z.all);
        iMcomb.microscope.depthFilter.swivel.z.P = vertcat(iMcomb.microscope.depthFilter.swivel.z.P, wkiM.microscope.depthFilter.swivel.z.P);
        iMcomb.microscope.depthFilter.swivel.z.AP = vertcat(iMcomb.microscope.depthFilter.swivel.z.AP, wkiM.microscope.depthFilter.swivel.z.AP);
        iMcomb.microscope.depthFilter.swivel.z.N = vertcat(iMcomb.microscope.depthFilter.swivel.z.N, wkiM.microscope.depthFilter.swivel.z.N);
        iMcomb.microscope.depthFilter.swivel.z.S = vertcat(iMcomb.microscope.depthFilter.swivel.z.S, wkiM.microscope.depthFilter.swivel.z.S);
        % depthFilter.swivel.threeD all P AP N S
        iMcomb.microscope.depthFilter.swivel.threeD.all = vertcat(iMcomb.microscope.depthFilter.swivel.threeD.all, wkiM.microscope.depthFilter.swivel.threeD.all);
        iMcomb.microscope.depthFilter.swivel.threeD.P = vertcat(iMcomb.microscope.depthFilter.swivel.threeD.P, wkiM.microscope.depthFilter.swivel.threeD.P);
        iMcomb.microscope.depthFilter.swivel.threeD.AP = vertcat(iMcomb.microscope.depthFilter.swivel.threeD.AP, wkiM.microscope.depthFilter.swivel.threeD.AP);
        iMcomb.microscope.depthFilter.swivel.threeD.N = vertcat(iMcomb.microscope.depthFilter.swivel.threeD.N, wkiM.microscope.depthFilter.swivel.threeD.N);
        iMcomb.microscope.depthFilter.swivel.threeD.S = vertcat(iMcomb.microscope.depthFilter.swivel.threeD.S, wkiM.microscope.depthFilter.swivel.threeD.S);
        % depthFilter.swivel.kMT
        iMcomb.microscope.depthFilter.swivel.kMT = vertcat(iMcomb.microscope.depthFilter.swivel.kMT, wkiM.microscope.depthFilter.swivel.kMT);
    end

    %plate
    if ~isfield(iMcomb, 'plate')
        iMcomb.plate = wkiM.plate;
    else
        %% coords x y z
        iMcomb.plate.coords.x = vertcat(iMcomb.plate.coords.x, wkiM.plate.coords.x);
        iMcomb.plate.coords.y = vertcat(iMcomb.plate.coords.y, wkiM.plate.coords.y);
        iMcomb.plate.coords.z = vertcat(iMcomb.plate.coords.z, wkiM.plate.coords.z);
        %% sisSep x y z twoD threeD
        iMcomb.plate.sisSep.x = vertcat(iMcomb.plate.sisSep.x, wkiM.plate.sisSep.x);
        iMcomb.plate.sisSep.y = vertcat(iMcomb.plate.sisSep.y, wkiM.plate.sisSep.y);
        iMcomb.plate.sisSep.z = vertcat(iMcomb.plate.sisSep.z, wkiM.plate.sisSep.z);
        iMcomb.plate.sisSep.twoD = vertcat(iMcomb.plate.sisSep.twoD, wkiM.plate.sisSep.twoD);
        iMcomb.plate.sisSep.threeD = vertcat(iMcomb.plate.sisSep.threeD, wkiM.plate.sisSep.threeD);
        %% raw.delta
        % raw.delta.x all P AP N S
        iMcomb.plate.raw.delta.x.all = vertcat(iMcomb.plate.raw.delta.x.all, wkiM.plate.raw.delta.x.all);
        iMcomb.plate.raw.delta.x.P = vertcat(iMcomb.plate.raw.delta.x.P, wkiM.plate.raw.delta.x.P);
        iMcomb.plate.raw.delta.x.AP = vertcat(iMcomb.plate.raw.delta.x.AP, wkiM.plate.raw.delta.x.AP);
        iMcomb.plate.raw.delta.x.N = vertcat(iMcomb.plate.raw.delta.x.N, wkiM.plate.raw.delta.x.N);
        iMcomb.plate.raw.delta.x.S = vertcat(iMcomb.plate.raw.delta.x.S, wkiM.plate.raw.delta.x.S);
        % raw.delta.y all P AP N S
        iMcomb.plate.raw.delta.y.all = vertcat(iMcomb.plate.raw.delta.y.all, wkiM.plate.raw.delta.y.all);
        iMcomb.plate.raw.delta.y.P = vertcat(iMcomb.plate.raw.delta.y.P, wkiM.plate.raw.delta.y.P);
        iMcomb.plate.raw.delta.y.AP = vertcat(iMcomb.plate.raw.delta.y.AP, wkiM.plate.raw.delta.y.AP);
        iMcomb.plate.raw.delta.y.N = vertcat(iMcomb.plate.raw.delta.y.N, wkiM.plate.raw.delta.y.N);
        iMcomb.plate.raw.delta.y.S = vertcat(iMcomb.plate.raw.delta.y.S, wkiM.plate.raw.delta.y.S);
        % raw.delta.z all P AP N S
        iMcomb.plate.raw.delta.z.all = vertcat(iMcomb.plate.raw.delta.z.all, wkiM.plate.raw.delta.z.all);
        iMcomb.plate.raw.delta.z.P = vertcat(iMcomb.plate.raw.delta.z.P, wkiM.plate.raw.delta.z.P);
        iMcomb.plate.raw.delta.z.AP = vertcat(iMcomb.plate.raw.delta.z.AP, wkiM.plate.raw.delta.z.AP);
        iMcomb.plate.raw.delta.z.N = vertcat(iMcomb.plate.raw.delta.z.N, wkiM.plate.raw.delta.z.N);
        iMcomb.plate.raw.delta.z.S = vertcat(iMcomb.plate.raw.delta.z.S, wkiM.plate.raw.delta.z.S);
        % raw.delta.oneD
        iMcomb.plate.raw.delta.oneD = vertcat(iMcomb.plate.raw.delta.oneD, wkiM.plate.raw.delta.oneD);
        % raw.delta.twoD all P AP N S
        iMcomb.plate.raw.delta.twoD.all = vertcat(iMcomb.plate.raw.delta.twoD.all, wkiM.plate.raw.delta.twoD.all);
        iMcomb.plate.raw.delta.twoD.P = vertcat(iMcomb.plate.raw.delta.twoD.P, wkiM.plate.raw.delta.twoD.P);
        iMcomb.plate.raw.delta.twoD.AP = vertcat(iMcomb.plate.raw.delta.twoD.AP, wkiM.plate.raw.delta.twoD.AP);
        iMcomb.plate.raw.delta.twoD.N = vertcat(iMcomb.plate.raw.delta.twoD.N, wkiM.plate.raw.delta.twoD.N);
        iMcomb.plate.raw.delta.twoD.S = vertcat(iMcomb.plate.raw.delta.twoD.S, wkiM.plate.raw.delta.twoD.S);
        % raw.delta.threeD all P AP N S
        iMcomb.plate.raw.delta.threeD.all = vertcat(iMcomb.plate.raw.delta.threeD.all, wkiM.plate.raw.delta.threeD.all);
        iMcomb.plate.raw.delta.threeD.P = vertcat(iMcomb.plate.raw.delta.threeD.P, wkiM.plate.raw.delta.threeD.P);
        iMcomb.plate.raw.delta.threeD.AP = vertcat(iMcomb.plate.raw.delta.threeD.AP, wkiM.plate.raw.delta.threeD.AP);
        iMcomb.plate.raw.delta.threeD.N = vertcat(iMcomb.plate.raw.delta.threeD.N, wkiM.plate.raw.delta.threeD.N);
        iMcomb.plate.raw.delta.threeD.S = vertcat(iMcomb.plate.raw.delta.threeD.S, wkiM.plate.raw.delta.threeD.S);
        %% raw.swivel
        % raw.swivel.y all P AP N S
        iMcomb.plate.raw.swivel.y.all = vertcat(iMcomb.plate.raw.swivel.y.all, wkiM.plate.raw.swivel.y.all);
        iMcomb.plate.raw.swivel.y.P = vertcat(iMcomb.plate.raw.swivel.y.P, wkiM.plate.raw.swivel.y.P);
        iMcomb.plate.raw.swivel.y.AP = vertcat(iMcomb.plate.raw.swivel.y.AP, wkiM.plate.raw.swivel.y.AP);
        iMcomb.plate.raw.swivel.y.N = vertcat(iMcomb.plate.raw.swivel.y.N, wkiM.plate.raw.swivel.y.N);
        iMcomb.plate.raw.swivel.y.S = vertcat(iMcomb.plate.raw.swivel.y.S, wkiM.plate.raw.swivel.y.S);
        % raw.swivel.z all P AP N S
        iMcomb.plate.raw.swivel.z.all = vertcat(iMcomb.plate.raw.swivel.z.all, wkiM.plate.raw.swivel.z.all);
        iMcomb.plate.raw.swivel.z.P = vertcat(iMcomb.plate.raw.swivel.z.P, wkiM.plate.raw.swivel.z.P);
        iMcomb.plate.raw.swivel.z.AP = vertcat(iMcomb.plate.raw.swivel.z.AP, wkiM.plate.raw.swivel.z.AP);
        iMcomb.plate.raw.swivel.z.N = vertcat(iMcomb.plate.raw.swivel.z.N, wkiM.plate.raw.swivel.z.N);
        iMcomb.plate.raw.swivel.z.S = vertcat(iMcomb.plate.raw.swivel.z.S, wkiM.plate.raw.swivel.z.S);
        % raw.swivel.threeD all P AP N S
        iMcomb.plate.raw.swivel.threeD.all = vertcat(iMcomb.plate.raw.swivel.threeD.all, wkiM.plate.raw.swivel.threeD.all);
        iMcomb.plate.raw.swivel.threeD.P = vertcat(iMcomb.plate.raw.swivel.threeD.P, wkiM.plate.raw.swivel.threeD.P);
        iMcomb.plate.raw.swivel.threeD.AP = vertcat(iMcomb.plate.raw.swivel.threeD.AP, wkiM.plate.raw.swivel.threeD.AP);
        iMcomb.plate.raw.swivel.threeD.N = vertcat(iMcomb.plate.raw.swivel.threeD.N, wkiM.plate.raw.swivel.threeD.N);
        iMcomb.plate.raw.swivel.threeD.S = vertcat(iMcomb.plate.raw.swivel.threeD.S, wkiM.plate.raw.swivel.threeD.S);
        % raw.swivel.kMT
        iMcomb.plate.raw.swivel.kMT = vertcat(iMcomb.plate.raw.swivel.kMT, wkiM.plate.raw.swivel.kMT);
        %% depthFilter.delta
        % depthFilter.delta.x all P AP N S
        iMcomb.plate.depthFilter.delta.x.all = vertcat(iMcomb.plate.depthFilter.delta.x.all, wkiM.plate.depthFilter.delta.x.all);
        iMcomb.plate.depthFilter.delta.x.P = vertcat(iMcomb.plate.depthFilter.delta.x.P, wkiM.plate.depthFilter.delta.x.P);
        iMcomb.plate.depthFilter.delta.x.AP = vertcat(iMcomb.plate.depthFilter.delta.x.AP, wkiM.plate.depthFilter.delta.x.AP);
        iMcomb.plate.depthFilter.delta.x.N = vertcat(iMcomb.plate.depthFilter.delta.x.N, wkiM.plate.depthFilter.delta.x.N);
        iMcomb.plate.depthFilter.delta.x.S = vertcat(iMcomb.plate.depthFilter.delta.x.S, wkiM.plate.depthFilter.delta.x.S);
        % depthFilter.delta.y all P AP N S
        iMcomb.plate.depthFilter.delta.y.all = vertcat(iMcomb.plate.depthFilter.delta.y.all, wkiM.plate.depthFilter.delta.y.all);
        iMcomb.plate.depthFilter.delta.y.P = vertcat(iMcomb.plate.depthFilter.delta.y.P, wkiM.plate.depthFilter.delta.y.P);
        iMcomb.plate.depthFilter.delta.y.AP = vertcat(iMcomb.plate.depthFilter.delta.y.AP, wkiM.plate.depthFilter.delta.y.AP);
        iMcomb.plate.depthFilter.delta.y.N = vertcat(iMcomb.plate.depthFilter.delta.y.N, wkiM.plate.depthFilter.delta.y.N);
        iMcomb.plate.depthFilter.delta.y.S = vertcat(iMcomb.plate.depthFilter.delta.y.S, wkiM.plate.depthFilter.delta.y.S);
        % depthFilter.delta.z all P AP N S
        iMcomb.plate.depthFilter.delta.z.all = vertcat(iMcomb.plate.depthFilter.delta.z.all, wkiM.plate.depthFilter.delta.z.all);
        iMcomb.plate.depthFilter.delta.z.P = vertcat(iMcomb.plate.depthFilter.delta.z.P, wkiM.plate.depthFilter.delta.z.P);
        iMcomb.plate.depthFilter.delta.z.AP = vertcat(iMcomb.plate.depthFilter.delta.z.AP, wkiM.plate.depthFilter.delta.z.AP);
        iMcomb.plate.depthFilter.delta.z.N = vertcat(iMcomb.plate.depthFilter.delta.z.N, wkiM.plate.depthFilter.delta.z.N);
        iMcomb.plate.depthFilter.delta.z.S = vertcat(iMcomb.plate.depthFilter.delta.z.S, wkiM.plate.depthFilter.delta.z.S);
        % depthFilter.delta.oneD
        iMcomb.plate.depthFilter.delta.oneD = vertcat(iMcomb.plate.depthFilter.delta.oneD, wkiM.plate.depthFilter.delta.oneD);
        % depthFilter.delta.twoD all P AP N S
        iMcomb.plate.depthFilter.delta.twoD.all = vertcat(iMcomb.plate.depthFilter.delta.twoD.all, wkiM.plate.depthFilter.delta.twoD.all);
        iMcomb.plate.depthFilter.delta.twoD.P = vertcat(iMcomb.plate.depthFilter.delta.twoD.P, wkiM.plate.depthFilter.delta.twoD.P);
        iMcomb.plate.depthFilter.delta.twoD.AP = vertcat(iMcomb.plate.depthFilter.delta.twoD.AP, wkiM.plate.depthFilter.delta.twoD.AP);
        iMcomb.plate.depthFilter.delta.twoD.N = vertcat(iMcomb.plate.depthFilter.delta.twoD.N, wkiM.plate.depthFilter.delta.twoD.N);
        iMcomb.plate.depthFilter.delta.twoD.S = vertcat(iMcomb.plate.depthFilter.delta.twoD.S, wkiM.plate.depthFilter.delta.twoD.S);
        % depthFilter.delta.threeD all P AP N S
        iMcomb.plate.depthFilter.delta.threeD.all = vertcat(iMcomb.plate.depthFilter.delta.threeD.all, wkiM.plate.depthFilter.delta.threeD.all);
        iMcomb.plate.depthFilter.delta.threeD.P = vertcat(iMcomb.plate.depthFilter.delta.threeD.P, wkiM.plate.depthFilter.delta.threeD.P);
        iMcomb.plate.depthFilter.delta.threeD.AP = vertcat(iMcomb.plate.depthFilter.delta.threeD.AP, wkiM.plate.depthFilter.delta.threeD.AP);
        iMcomb.plate.depthFilter.delta.threeD.N = vertcat(iMcomb.plate.depthFilter.delta.threeD.N, wkiM.plate.depthFilter.delta.threeD.N);
        iMcomb.plate.depthFilter.delta.threeD.S = vertcat(iMcomb.plate.depthFilter.delta.threeD.S, wkiM.plate.depthFilter.delta.threeD.S);
        %% depthFilter.swivel
        % depthFilter.swivel.y all P AP N S
        iMcomb.plate.depthFilter.swivel.y.all = vertcat(iMcomb.plate.depthFilter.swivel.y.all, wkiM.plate.depthFilter.swivel.y.all);
        iMcomb.plate.depthFilter.swivel.y.P = vertcat(iMcomb.plate.depthFilter.swivel.y.P, wkiM.plate.depthFilter.swivel.y.P);
        iMcomb.plate.depthFilter.swivel.y.AP = vertcat(iMcomb.plate.depthFilter.swivel.y.AP, wkiM.plate.depthFilter.swivel.y.AP);
        iMcomb.plate.depthFilter.swivel.y.N = vertcat(iMcomb.plate.depthFilter.swivel.y.N, wkiM.plate.depthFilter.swivel.y.N);
        iMcomb.plate.depthFilter.swivel.y.S = vertcat(iMcomb.plate.depthFilter.swivel.y.S, wkiM.plate.depthFilter.swivel.y.S);
        % depthFilter.swivel.z all P AP N S
        iMcomb.plate.depthFilter.swivel.z.all = vertcat(iMcomb.plate.depthFilter.swivel.z.all, wkiM.plate.depthFilter.swivel.z.all);
        iMcomb.plate.depthFilter.swivel.z.P = vertcat(iMcomb.plate.depthFilter.swivel.z.P, wkiM.plate.depthFilter.swivel.z.P);
        iMcomb.plate.depthFilter.swivel.z.AP = vertcat(iMcomb.plate.depthFilter.swivel.z.AP, wkiM.plate.depthFilter.swivel.z.AP);
        iMcomb.plate.depthFilter.swivel.z.N = vertcat(iMcomb.plate.depthFilter.swivel.z.N, wkiM.plate.depthFilter.swivel.z.N);
        iMcomb.plate.depthFilter.swivel.z.S = vertcat(iMcomb.plate.depthFilter.swivel.z.S, wkiM.plate.depthFilter.swivel.z.S);
        % depthFilter.swivel.threeD all P AP N S
        iMcomb.plate.depthFilter.swivel.threeD.all = vertcat(iMcomb.plate.depthFilter.swivel.threeD.all, wkiM.plate.depthFilter.swivel.threeD.all);
        iMcomb.plate.depthFilter.swivel.threeD.P = vertcat(iMcomb.plate.depthFilter.swivel.threeD.P, wkiM.plate.depthFilter.swivel.threeD.P);
        iMcomb.plate.depthFilter.swivel.threeD.AP = vertcat(iMcomb.plate.depthFilter.swivel.threeD.AP, wkiM.plate.depthFilter.swivel.threeD.AP);
        iMcomb.plate.depthFilter.swivel.threeD.N = vertcat(iMcomb.plate.depthFilter.swivel.threeD.N, wkiM.plate.depthFilter.swivel.threeD.N);
        iMcomb.plate.depthFilter.swivel.threeD.S = vertcat(iMcomb.plate.depthFilter.swivel.threeD.S, wkiM.plate.depthFilter.swivel.threeD.S);
        % depthFilter.swivel.kMT
        iMcomb.plate.depthFilter.swivel.kMT = vertcat(iMcomb.plate.depthFilter.swivel.kMT, wkiM.plate.depthFilter.swivel.kMT);
        %% twist y z threeD
        iMcomb.plate.twist.y = vertcat(iMcomb.plate.twist.y, wkiM.plate.twist.y);
        iMcomb.plate.twist.z = vertcat(iMcomb.plate.twist.z, wkiM.plate.twist.z);
        iMcomb.plate.twist.threeD = vertcat(iMcomb.plate.twist.threeD, wkiM.plate.twist.threeD);
        %% sisterCentreSpeed
        iMcomb.plate.sisterCentreSpeed = vertcat(iMcomb.plate.sisterCentreSpeed, wkiM.plate.sisterCentreSpeed);
        %% plateThickness
        iMcomb.plate.plateThickness = vertcat(iMcomb.plate.plateThickness, wkiM.plate.plateThickness);
    end

    %intensity
    if ~isfield(iMcomb, 'intensity')
        iMcomb.intensity = wkiM.intensity;
    else
        %% mean inner outer
        iMcomb.intensity.mean.inner = vertcat(iMcomb.intensity.mean.inner, wkiM.intensity.mean.inner);
        iMcomb.intensity.mean.outer = vertcat(iMcomb.intensity.mean.outer, wkiM.intensity.mean.outer);
        %% max inner outer
        iMcomb.intensity.max.inner = vertcat(iMcomb.intensity.max.inner, wkiM.intensity.max.inner);
        iMcomb.intensity.max.outer = vertcat(iMcomb.intensity.max.outer, wkiM.intensity.max.outer);
        %% bg inner outer
        iMcomb.intensity.bg.inner = vertcat(iMcomb.intensity.bg.inner, wkiM.intensity.bg.inner);
        iMcomb.intensity.bg.outer = vertcat(iMcomb.intensity.bg.outer, wkiM.intensity.bg.outer);
    end

    %intiMcomb structure
    if ~isfield(intiMcomb, 'intensity')
        intiMcomb.intensity = wkintiM.intensity;
        %remember to change the mean outer intensity!
        intiMcomb.intensity.mean.outer = wkIntensity;
    else
        %% mean inner outer
        intiMcomb.intensity.mean.inner = vertcat(intiMcomb.intensity.mean.inner, wkintiM.intensity.mean.inner);
        %remember to change the mean outer intensity!
        intiMcomb.intensity.mean.outer = vertcat(intiMcomb.intensity.mean.outer, wkIntensity);
        %% max inner outer
        intiMcomb.intensity.max.inner = vertcat(intiMcomb.intensity.max.inner, wkintiM.intensity.max.inner);
        intiMcomb.intensity.max.outer = vertcat(intiMcomb.intensity.max.outer, wkintiM.intensity.max.outer);
        %% bg inner outer
        intiMcomb.intensity.bg.inner = vertcat(intiMcomb.intensity.bg.inner, wkintiM.intensity.bg.inner);
        intiMcomb.intensity.bg.outer = vertcat(intiMcomb.intensity.bg.outer, wkintiM.intensity.bg.outer);
    end

end


end

%%
function options = processOptions(defaults, varargin)
% PROCESSOPTIONS Process option pairs into struct
%
% Takes a struct containing default values and a list of string/value pairs and
% updates the defaults according to the options found in the pairs. Case
% insensitive.
%
% Copyright (c) 2013 Jonathan Armond

options = defaults;
fields = fieldnames(options);
i = 1;
while i <= length(varargin)
  optname = varargin{i};
  if ~ischar(optname)
    error(['Expected string for parameter ' num2str(i)]);
  end

  % Find corresponding field name.
  idx = find(strcmpi(fields,optname));
  if isempty(idx)
    error(['Unrecognized option ''' optname '''']);
  end
  field = fields{idx};
  
  if i+1 > length(varargin)
    error(['Expected value to follow ''' optname '''']);
  end
  optvalue = varargin{i+1};

  % Store option value in struct.
  options.(field) = optvalue;
  
  i = i+2;
end

end

%%
function binnedMat = makeBinsintiMs(intiMToBin, nBins, marker)
    if strcmp(marker, 'outer')
        toBin = intiMToBin.intensity.mean.outer;
    elseif strcmp(marker, 'outerDivInner')
        toBin = intiMToBin.intensity.mean.outer ./ intiMToBin.intensity.mean.inner;
    elseif strcmp(marker, 'inner')
        toBin = intiMToBin.intensity.mean.inner;
    end
    toBin = toBin(:);
    toBin = rmmissing(toBin);
    toBin = sort(toBin, 'descend'); %edit 20250430
    total = length(toBin);
    stepsize = floor(total/nBins);
    first = 1;
    binnedMat = nan(stepsize, nBins);

    for a = 1:nBins
        last = first + stepsize - 1;
        
        binnedMat(:, a) = toBin(first:last);

        first = last + 1;

    end

end
%%
function newInner = rmLowCenpC(intiM)
nKTs = length(intiM.intensity.mean.inner);
if size(intiM.intensity.mean.inner, 2) < 1
    paired = 1;
else
    paired = 0;
end

for iKT = 1:nKTs
    bgTimesTwo = intiM.intensity.bg.inner(iKT)*2;
    if ~isnan(intiM.intensity.mean.inner(iKT,1))
        if intiM.intensity.mean.inner(iKT,1)<bgTimesTwo
            intiM.intensity.mean.inner(iKT,1) = nan;
        end
    end
    if paired
        if ~isnan(intiM.intensity.mean.inner(iKT,2))
            if intiM.intensity.mean.inner(iKT,2)<bgTimesTwo
                intiM.intensity.mean.inner(iKT,2) = nan;
            end
        end
    end
end
newInner = intiM.intensity.mean.inner;
end
