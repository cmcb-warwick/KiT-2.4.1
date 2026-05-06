function intData = makeSingleIntPlots(ctrlIntiMs, varargin)
% MAKESINGLEINTPLOTS Produces a plot of kinetchore intensity
% measurements over multiple cells and experiments split into bins.
% 
%   MAKESINGLEINTPLOTS(ctrlIntiMs).
%   ctrlIntiMs must be organised as {Mad2intiM1; Mad2intiM2; ...}.
%
%   Options are available:
%
%    Options, defaults in {}:-
%
%   nBins: {25} or other integer. The number of bins to split each intiM
%       structure into.
%
%   normalise: 0 or {1}. Whether to normalise the median of the maximum bin for
%       each intensity class to one. 
%
%   normIndInt: 0 or {1}. Whether to normalise the independent intensity
%       marker (i.e. Mad2) to CenpC or not.
%
%   indMarker: {nan} or string of intensity marker that you binned on. Will
%       request this later if not filled in now.
%
%   minBinFirst: {0} or 1. Whether to have maximum intensity as maximum (0) or
%       minimum (1) bin number when plotting.
%
% Copyright (c) 2025 C. C. Conway
opts.nBins = 25;
opts.normalise = 1;
opts.normIndInt = 1;
opts.indMarker = 'Venus-Mad2';
opts.minBinFirst = 0;
opts.calcHalfChange = 0; %not yet implemented
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = length(ctrlIntiMs);

indCtrl = []; %idea: instead of using the binned data, use the binning to normalise all data, then get all Mad2 data out

for iExpt = 1:nExpts
    %test from here: removing CenpC data if it's not at least 2x background
    ctrlIntiMs{iExpt}.intensity.mean.inner = rmLowCenpC(ctrlIntiMs{iExpt});
    indIntiM = ctrlIntiMs{iExpt}.intensity.mean.outer;
    if opts.normIndInt
        indIntiM = indIntiM ./ ctrlIntiMs{iExpt}.intensity.mean.inner;
    end
    indIntiM = indIntiM(:); %20250428
    indIntiM = rmmissing(indIntiM); %20250428
    iCBinned = makeBins(nBins, indIntiM);
    iCMed = median(iCBinned);
    
    
    if isNorm
        iMax = max(iCMed);
        iCBinned = iCBinned/iMax;
        indIntiM = indIntiM/iMax; %20250428
    end
    %indCtrl = cat(1, indCtrl, iCBinned);
    indCtrl = cat(1, indCtrl, indIntiM); %20250428

end

iCIntsGrouped = makeBins(nBins, indCtrl);
cBinSz = size(iCIntsGrouped,1);
[indCIntsMed, indCIntsNeg, indCIntsPos] = makeQuartiles(iCIntsGrouped);
[~, indCIntsLCI, indCIntsUCI] = makeMedianConfInts(iCIntsGrouped);

if isNorm
    indCMaxMed = max(indCIntsMed);
    iCIntsGrouped = iCIntsGrouped/indCMaxMed;
    indCIntsMed = indCIntsMed/indCMaxMed;
    indCIntsNeg = indCIntsNeg/indCMaxMed;
    indCIntsPos = indCIntsPos/indCMaxMed;
    indCIntsLCI = indCIntsLCI/indCMaxMed;
    indCIntsUCI = indCIntsUCI/indCMaxMed;
end

intData.intensities = iCIntsGrouped;
CIstats.median = indCIntsMed;
CIstats.lowerCI = indCIntsLCI;
CIstats.upperCI = indCIntsUCI;
CIstats.lowerQuart = indCIntsNeg;
CIstats.upperQuart = indCIntsPos;
intData.stats = CIstats;


nCells_ctrl = getnCells(ctrlIntiMs);
intData.nCells = nCells_ctrl;
intData.binSize = cBinSz;



%%

ctrlXVals = 1:nBins;

indCols = [0 0 0;...
    0.4 0.4 0.4];
%black and grey

%% Make figure for control
if isnan(opts.indMarker)
    labelAnswers = inputdlg({'Enter independent intensity marker:'}, 'Marker proteins',[1 35], {'Mad2'});
else
    labelAnswers = {opts.indMarker};
end
ctrlFig = figure;
hold on
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [indCIntsNeg fliplr(indCIntsPos)], indCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [indCIntsLCI fliplr(indCIntsUCI)], indCols(1,:), ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
plot(ctrlXVals, indCIntsMed, '-', 'Color', indCols(1,:), ...
    'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);

hold off
ctrlax = gca;
IndCtrlTicks = yticks;
upperCtrlLimit = max(IndCtrlTicks(:));
lowerCtrlLimit = min(IndCtrlTicks(:));

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;

if opts.normIndInt
    indDivSignaller = '/CenpC';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'YLim', [lowerCtrlLimit upperCtrlLimit], 'TickDir', 'out', 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'YGrid', 'on', 'YMinorGrid', 'on');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'YGrid', 'on', 'YMinorGrid', 'on');
end
ylabel(labelInd);


if opts.minBinFirst
    set(gca, 'XDir', 'reverse');
end

xlabel(sprintf('%s%s pseudo-timeline bins', labelAnswers{1}, indDivSignaller));

   
%% Resizing windows and axes
if isNorm

        minValue = min(indCIntsNeg, [], 'all');
        possVec = [lowerCtrlLimit, -0.5, -0.2, -0.1, -0.05];
    
    possVec = sort(possVec);
    if minValue < possVec(2)
        lowerCtrlLimit = possVec(1);
    elseif minValue < possVec(3)
        lowerCtrlLimit = possVec(2);
    elseif minValue < possVec(4)
        lowerCtrlLimit = possVec(3);
    elseif minValue < possVec(5)
        lowerCtrlLimit = possVec(4);
    else
        lowerCtrlLimit = possVec(5);
    end
    figure(ctrlFig);
    set(ctrlax, 'YLim', [lowerCtrlLimit upperCtrlLimit]);
    ctrlRange = upperCtrlLimit - lowerCtrlLimit;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    %%set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, 275, ctrlRange*200])
    set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlRange*200])
    set(ctrlFig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (ctrlRange*200)+(ctrlPos(2)*2)])
    set(ctrlax, 'XTickLabelRotation', 0)
    
    
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
function [medianVal, lowerQuart, upperQuart] = makeQuartiles(binnedMatrix)
medianVal = median(binnedMatrix);
lowerQuart = quantile(binnedMatrix, 0.25);
upperQuart = quantile(binnedMatrix, 0.75);
end
%%
function [medianVal, lowerCI, upperCI] = makeMedianConfInts(binnedMatrix)
medianVal = median(binnedMatrix);
nVals = size(binnedMatrix, 1);
lPos = (nVals*0.5) - (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101 - first two 0.5 values indicate we are assessing median quantile, third 0.5 is from 1-0.5 (probability of not 0.5).
uPos = (nVals*0.5) + (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101
lPct = lPos/nVals;
uPct = uPos/nVals;
lowerCI = quantile(binnedMatrix, lPct);
upperCI = quantile(binnedMatrix, uPct);
end


%%
function binnedMat = makeBins(nBins, toBin)
    toBin = toBin(:);
    toBin = rmmissing(toBin);
    toBin = sort(toBin,'descend'); %CCC 2025.04.28 change from 'ascend', this will fit better with later binning
    total = length(toBin);
    stepsize = floor(total/nBins);
    %remainder = mod(total, nBins); %CCC 2025.04.28 to calculate number of non included data points
    %to decide: exclude more from the lower intensities or higher
    %intensities? Lower intensities are probably more consistent so keep
    %preferentially?
    %excludeHigh = ceil(remainder/2); %this way if there are uneven numbers of KTs to exclude, I favour excluding the higher intensity KT. To favour lower intensity exclusion, change ceil to floor.
    %first = excludeHigh + 1; %in case remainder is 0, then start at first position.
    %!! Excluded above on advice of NJB as very low Mad2 values are in Mad2
    %noise and are therefore just stochastic data !!
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