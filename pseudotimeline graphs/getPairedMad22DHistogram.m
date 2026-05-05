function data = getPairedMad22DHistogram(intiMs, varargin)
% GETPAIREDMAD22DHISTOGRAM gets 2D histogram showing which bin each sister
% in kinetochore pair occupies
%
%   GETPAIREDMAD22DHISTOGRAM(intiMs).
%   intiMs must be organised as {intiM1, intiM2, intiM3, ...}.
%
% Copyright (c) 2026 C. C. Conway
%
% options below:

opts.nBins = 25; %any integer
opts.normIndInt = 1; %normalise outer to inner intensity (i.e. Mad2 to CenpC)
opts.normalise = 1; %normalise to 1 or not. Keep as 1 when assessing multiple experiments.
opts.swapBin = 0; %if you want bin 25 to be highest intensity bin.
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = length(intiMs);
allIntiMScaled = [];
allcMax = [];
allPairs = [];

for iExpt = 1:nExpts
    %step 1: remove CenpC intensities that are too low
    [newInner, newOuter] = rmLowCenpC(intiMs{iExpt});
    intiMs{iExpt}.intensity.mean.inner = newInner;
    intiMs{iExpt}.intensity.mean.outer = newOuter;
    cIntiM = intiMs{iExpt}.intensity.mean.outer;
    if opts.normIndInt
        cIntiM = cIntiM ./ intiMs{iExpt}.intensity.mean.inner;
    end
    %then bin (max values retained)
    ctrlIntiM = makeBins(nBins, cIntiM);
    cMed = median(ctrlIntiM);
    
    wkcIntiM = rmmissing(cIntiM);
    
    if isNorm
        %get max median of ctrl
        cMax = max(cMed);
        %put in a matrix so we know what to divide by later on
        allcMax = cat(2, allcMax, cMax);
        %normalise intiM
        ctrlIntiM = ctrlIntiM/cMax;
        wkcIntiM = wkcIntiM/cMax;

    end
    allIntiMScaled = cat(1, allIntiMScaled, ctrlIntiM);
    allPairs = cat(1, allPairs, wkcIntiM);
end

ctrlIntiMScaledBinned = makeBins(nBins, allIntiMScaled);
binSz = size(ctrlIntiMScaledBinned, 1);
pairBins = [];
%toFindUnique = [];
for iPair = 1:length(allPairs)
    wkPair = [];
    bothFine = 1; %assume both in pair have been accepted into KT matrix

    for iSis = 1:2
        [row, col] = find(ctrlIntiMScaledBinned == allPairs(iPair, iSis));
        if ~isempty(col)
            if length(col) == 1
                wkPos = (col-0.5)+(row/binSz);
                wkPair = cat(2, wkPair, wkPos);
            else
                wkPos = (col(1)-0.5)+(row(1)/binSz);
                wkPair = cat(2, wkPair, wkPos);
            end
        else
            bothFine = 0;
        end
    end

    if bothFine
        wkPair = sort(wkPair);
        pairBins = cat(1, pairBins, wkPair);
        %wkUnique = (wkPair(1)^1.5)*((wkPair(2)^2)+1); %horrendous calculation, basically making sure that row1 column2 isn't equivalent to row2 column1 in the plotting. for 1:25, none of these values overlap
        %toFindUnique = cat(1, toFindUnique, wkUnique);
    end

end

data.nCells = getnCells(intiMs);
data.nKTs_all = length(ctrlIntiMScaledBinned(:));
data.nPairs = length(pairBins);
data.nExpts = nExpts;

%%
figure
hold on
colormap('jet')
numPerCoord = histcounts2(pairBins(:,1), pairBins(:,2), 0.5:1:(nBins+0.5), 0.5:1:(nBins+0.5));
data.histCounts = numPerCoord;
%coordDiag = diag(diag(numPerCoord));
%dataFor2DHist = numPerCoord + transpose(numPerCoord);
%dataFor2DHist = dataFor2DHist - coordDiag;
h = histogram2('XBinEdges', [0.5:1:(nBins+0.5)], 'YBinEdges', [0.5:1:(nBins+0.5)], 'BinCounts', numPerCoord, 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
plot([0 26], [0 26])
xlim([0.5 nBins+0.5])
ylim([0.5 nBins+0.5])
axis square
xticks([1:nBins])
yticks([1:nBins])
set(gca, 'TickDir', 'out', 'XTickLabelRotation', 0)
colorbar
end




%% processOptions
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
function [newInner, newOuter] = rmLowCenpC(intiM)
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
            intiM.intensity.mean.outer(iKT,1) = nan;
        end
    end
    if paired
        if ~isnan(intiM.intensity.mean.inner(iKT,2))
            if intiM.intensity.mean.inner(iKT,2)<bgTimesTwo
                intiM.intensity.mean.inner(iKT,2) = nan;
                intiM.intensity.mean.outer(iKT,2) = nan;
            end
        end
    end
end
newInner = intiM.intensity.mean.inner;
newOuter = intiM.intensity.mean.outer;
end

%%
function binnedMat = makeBins(nBins, toBin)
    toBin = toBin(:);
    toBin = rmmissing(toBin);
    toBin = sort(toBin,'descend'); %CCC 2025.04.28 change from 'ascend', this will fit better with later binning
    total = length(toBin);
    stepsize = floor(total/nBins);
    first = 1;
    binnedMat = nan(stepsize, nBins);

    for a = 1:nBins
        last = first + stepsize - 1;
        
        binnedMat(:, a) = toBin(first:last);

        first = last + 1;

    end

    %the matrix will go from max to min intensity (left to right)

end

%%
function nCells = getnCells(intiMs)
nExpts = length(intiMs);
nCells = 0;
for iExpt = 1:nExpts
    CellIDs = [];
    cellLabel = intiMs{iExpt}.label;
    nKTs = size(cellLabel,1);
    for iKT = 1:nKTs
         CellIDs = cat(2, CellIDs, str2num(cellLabel(iKT,3:4)));
    end
    wkCells = unique(CellIDs);
    nCells = nCells + length(wkCells);
end
end