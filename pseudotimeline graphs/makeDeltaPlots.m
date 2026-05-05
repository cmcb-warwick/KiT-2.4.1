function deltaData = makeDeltaPlots(deltaIntiM, normIntiM, varargin)
% MAKEDELTAPLOTS Produces a single plot of multiple kinetchore delta
% measurements over multiple cells and experiments split into bins.
% 
%   MAKEDELTAPLOTS(deltaIntiM, normData).
%   deltaIntiM must be organised as {delta, stddev, poolMad2intiM}.
%   normIntiMs is optional so can be replaced with [].
%
%   Options are available:
%
%   Options, defaults in {}:-
%
%   normalise: 0 or {1}. Whether to normalise the median of the maximum bin for
%       each intensity class to one. Will normalise to normData or intiM
%       you have provided.
%
%   normIndInt: {0} or 1. Whether to normalise the independent intensity
%       marker (i.e. Mad2) to CenpC or not. Because the pooled data has been
%       binned on MAD2/CENPC intensity, this is generally 0.
%
%   normIndInt_ctrl: 0 or {1}. Whether to normalise the independent intensity
%       marker of the control data you've provided (i.e. Mad2) to CenpC or not.
%
%   intMarker: {nan} or string of intensity marker that you binned on. Will
%       request this later if not filled in now.
%
%   distPair: {nan} or string of proteins you are measuring delta between.
%       Will request this later if not filled in now.
%
%   swapBin: {0} or 1. Whether to have maximum intensity as minimum (0) or
%       maximum (1) bin number when plotting.
%
% Copyright (c) 2025 C. C. Conway

opts.normalise = 1;
opts.normIndInt = 0;
opts.normIndInt_ctrl = 1;
opts.intMarker = nan;
opts.distPair = nan;
opts.swapBin = 0;
opts = processOptions(opts, varargin{:});

nBins = length(deltaIntiM{1,1});
isNorm = opts.normalise;


intiM = deltaIntiM{3}.intensity.mean.outer;
if opts.normIndInt
    intiM = deltaIntiM{3}.intensity.mean.outer ./ deltaIntiM{3}.intensity.mean.inner;
end

ints = makeBins(nBins, intiM);

%get quartile data
[cMed, cLQ, cUQ] = makeQuartiles(ints);
[~, cLCI, cUCI] = makeMedianConfInts(ints);
cBinSz = size(ints, 1);


if isNorm
    allCtrlData = [];
    if ~isempty(normIntiM)
        for iExpt = 1:length(normIntiM)
            [innerCtrl, outerCtrl] = rmLowCenpC(normIntiM{iExpt});
            normIntiM{iExpt}.intensity.mean.inner = innerCtrl; %20250429
            normIntiM{iExpt}.intensity.mean.outer = outerCtrl;
            wkNormIntiM = outerCtrl;
            if opts.normIndInt_ctrl
                wkNormIntiM = wkNormIntiM ./ innerCtrl;
            end
            iCtrlBinned = makeBins(nBins,wkNormIntiM);
            iCtrlMax = max(median(iCtrlBinned));
            iCtrlBinned = iCtrlBinned/iCtrlMax;
            allCtrlData = cat(1, allCtrlData, iCtrlBinned);
        end
        allCtrlBinned = makeBins(nBins,allCtrlData);
        cMax = max(median(allCtrlBinned));
    else
        cMax = max(cMed);
    end

    
    ints = ints/cMax;
    cMed = cMed/cMax;
    cLQ = cLQ/cMax;
    cUQ = cUQ/cMax;
    cLCI = cLCI/cMax;
    cUCI = cUCI/cMax;
    
    
end

deltaData.ints.data = ints;
deltaData.ints.stats.median = cMed;
deltaData.ints.stats.lowerCI = cLCI;
deltaData.ints.stats.upperCI = cUCI;
deltaData.ints.stats.lowerQuart = cLQ;
deltaData.ints.stats.upperQuart = cUQ;



%%


ctrlXVals = 1:nBins;


indCols = [0 0 0;...
    0.4 0.4 0.4];
%black and grey

depCols = [0.8118 0.1333 0.5333;...
            0.4471 0.0667 0.4275;...
            0.9451 0.7176 0.8510;...
            0.3137 0.6667 0.1882;...
            0.1098 0.6745 0.3412;...
            0.7255 0.8784 0.5255];
%magenta, purple, light pink, green, dark green, light green

%control Delta with 1.96 x std dev (95 pct confidence interval)
cDeltaLower = deltaIntiM{1} - (1.96*deltaIntiM{2});
cDeltaUpper = deltaIntiM{1} + (1.96*deltaIntiM{2});

deltaData.delta.mean = transpose(deltaIntiM{1});
deltaData.delta.std = transpose(deltaIntiM{2});
deltaData.delta.lCI = transpose(cDeltaLower);
deltaData.delta.uCI = transpose(cDeltaUpper);

deltaData.nCells = getnCells_pool(deltaIntiM{3});
binSz = size(ints,1);
deltaData.binSize = binSz;
%now I want to get resampled Delta data for doing midpoint analysis later.
%I think I want to do 4 data runs (like I do with Delta) then randomly
%sample integers from this.
DeltaSample = [];
for iBin = 1:nBins
    randBin = [];
    for iRun = 1:4
        wkData = [];
        for iSample = 1:binSz
            wkSample = random('Normal',deltaIntiM{1}(iBin), deltaIntiM{2}(iBin));
            wkData = cat(1, wkData, wkSample);
        end
        randBin = cat(1, randBin, wkData);
    end
    selectData = randperm(binSz*4, binSz);
    selectData = sort(selectData);
    binData = randBin(selectData);
    DeltaSample = cat(2, DeltaSample, binData);
end

deltaData.delta.samples = DeltaSample;
halfChanges = getHalfChanges(deltaData, 'dataType', 'delta');
deltaData.ints.HCint = halfChanges.Mad2.HCintensity;
deltaData.ints.HCbin = halfChanges.Mad2.HCbin;
deltaData.delta.HCdelta = halfChanges.Other.HCintensityKKDelta;
deltaData.delta.HCbin = halfChanges.Other.HCbin;
deltaData.delta.HCMad2Int = halfChanges.Other.HCMad2Int;
%% Make figure 
if any([isnan(opts.intMarker), isnan(opts.distPair)])
    labelAnswers = inputdlg({'Enter independent intensity marker:', 'Enter distance between:'}, 'Marker proteins',[1 35], {'Mad2', 'CenpC to Ndc80N'});
else
    labelAnswers = {opts.intMarker, opts.distPair};
end
ctrlFig = figure;
hold on
yyaxis left
if opts.swapBin
    wkxVals = ctrlXVals - (nBins+1); %shift to negative 
    finxVals = abs(wkxVals); %make positive
    ctrlXVals = finxVals; %this is a better option than fliplr because that would mess up sub-bin localisation - this is less like mirroring the x-axis and more like just replacing the axis tick numbers
end
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [cLQ fliplr(cUQ)], indCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [cLCI fliplr(cUCI)], indCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence interval
plot(ctrlXVals, cMed, '-', 'Color', indCols(1,:), ...
    'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
line('XData',[deltaData.ints.HCbin, deltaData.ints.HCbin], 'YData',[0 1], 'Color', indCols(1,:), 'LineStyle', '--');
yyaxis right
% Delta plotted here
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [transpose(cDeltaLower) fliplr(transpose(cDeltaUpper))], depCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none')
plot(ctrlXVals, transpose(deltaIntiM{1}), '-', 'Color', depCols(1,:), ...
    'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
line('XData',[deltaData.delta.HCbin, deltaData.delta.HCbin], 'YData',[0 80], 'Color', depCols(1,:), 'LineStyle', '--');
hold off
ctrlax = gca;
yyaxis left
intCtrlTicks = yticks;
minIntCtrlTick = min(intCtrlTicks);
maxIntCtrlTick = max(intCtrlTicks);
yyaxis right

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on');
else
    labelInd = sprintf('%s%s intensity (AU)', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on');
end
ylabel(labelInd);

yyaxis right


    ctrlXTicks = 1:nBins;
    lowerX = min(ctrlXVals)-1;
    upperX = max(ctrlXVals)+1;
    

labelDep = sprintf('%s distance (nm)', labelAnswers{2});
set(gca, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on');
if opts.swapBin
    set(gca, 'XDir', 'reverse');
end

ylabel(labelDep);
xlabel(sprintf('%s%s pseudo-timeline bins', labelAnswers{1}, indDivSignaller));


%% Resizing windows and axes
%if you have treatment data this should rescale nicely. Otherwise shouldn't
%rescale. 20250715 test. Verified working with monastrol data 20250715
dividedMax = max(cMed)/5;
nearestRound = round(dividedMax,1);
multiplierTreat = nearestRound*5; %this scales graph to nearest 0.5
    
%int first
minIntValue = min(cLQ);
possIntVec = [minIntCtrlTick, -0.5, -0.2, -0.1, -0.05]*multiplierTreat;
%now Delta
minDeltaValue = min(cDeltaLower);

    %int first
    lowerIntLimit = [];
    possIntVec = sort(possIntVec);
    if minIntValue < possIntVec(2) || minDeltaValue < possIntVec(2)*75/multiplierTreat %changed from 90, added division
        lowerIntLimit = possIntVec(1);
    elseif minIntValue < possIntVec(3) || minDeltaValue < possIntVec(3)*75/multiplierTreat %changed from 90, added division
        lowerIntLimit = possIntVec(2);
    elseif minIntValue < possIntVec(4) || minDeltaValue < possIntVec(4)*75/multiplierTreat %changed from 90, added division
        lowerIntLimit = possIntVec(3);
    elseif minIntValue < possIntVec(5) || minDeltaValue < possIntVec(5)*75/multiplierTreat %changed from 90, added division
        lowerIntLimit = possIntVec(4);
    else
        lowerIntLimit = possIntVec(5);
    end
    
    lowerDeltaLimit = (lowerIntLimit*75/multiplierTreat); %changed from 90, added division
    figure(ctrlFig);
    yyaxis left
    set(ctrlax, 'YLim', [lowerIntLimit maxIntCtrlTick]);
    
    yyaxis right
    set(ctrlax, 'YLim', [lowerDeltaLimit (maxIntCtrlTick*75/multiplierTreat)]); %changed from 90, added division
    ctrlRange = maxIntCtrlTick - lowerIntLimit;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlRange*200/multiplierTreat]) %20250715, added division
    set(ctrlFig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (ctrlRange*200/multiplierTreat)+(ctrlPos(2)*2)]) %20250715, added division
    set(ctrlax, 'XTickLabelRotation', 0)
     
%end
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
function binnedMat = makeBins(nBins, toBin)
    toBin = toBin(:);
    toBin = rmmissing(toBin);
    toBin = sort(toBin, 'descend'); %20250502
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
function nCells = getnCells_pool(intiM)
    CellIDs = [];
    cellLabel = intiM.label;
    nKTs = size(cellLabel,1);
    for iKT = 1:nKTs
         CellIDs = cat(2, CellIDs, str2num(cellLabel(iKT,1:4)));
    end
    nCells = length(unique(CellIDs));

end
