function allStats = makeMultipleMad2CurveGraphs(intiMs, varargin)
%
% Works somewhat similarly to makeIntPlots, but only Mad2 intensities are
% plotted, and multiple Mad2 plots are overlaid on each other.
%
% intiMs should be organised as {Mad2intiM1, Mad2intiM2, ...}
%
% Copyright (c) 2025 C. C. Conway
%
opts.nBins = 25;
opts.normalise = 1; %normalise to inner or not
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = length(intiMs);
allStats.nExpts = nExpts;
allExptData = cell(1, nExpts);
nKTs = 0;
combinedData = [];
for iExpt = 1:nExpts
    intiM = intiMs{iExpt};
    intiM.intensity.mean.inner = rmLowCenpC(intiM);
    wkData = intiM.intensity.mean.outer;
    if isNorm
        wkData = intiM.intensity.mean.outer ./ intiM.intensity.mean.inner;
    end
    binnedData = makeBins(nBins, wkData);
    maxInt = max(median(binnedData));
    binnedData = binnedData/maxInt;
    combinedData = cat(1, combinedData, binnedData);
    allExptData{iExpt} = binnedData;
    nKTs = nKTs + numel(binnedData);
end
nCellsIntiMs = intiMs;
allStats.nCells = getnCells(nCellsIntiMs);
allStats.nKTs = nKTs;
newCombData = makeBins(nBins, combinedData);
newCombData = newCombData/max(median(newCombData));

%following colours are from matlab online orderedcolors
figCols = [0.0000    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    1.0000    0.8390    0.0390;...
    0.3960    0.5090    0.9920;...
    1.0000    0.2700    0.2270;...
    0.0000    0.6390    0.6390;...
    0.7960    0.5170    0.3640;...
    0.1490    0.5490    0.8660;...
    0.9600    0.4660    0.1600;...
    1.0000    0.9090    0.3920;...
    0.7520    0.3600    0.9840;...
    0.2860    0.8580    0.2500;...
    0.4230    0.9560    1.0000;...
    0.9490    0.4030    0.7720;...
    0.9960    0.7520    0.2980;...
    0.4900    0.6620    1.0000;...
    1.0000    0.4780    0.4540;...
    0.1210    0.8110    0.7450;...
    0.8620    0.6000    0.4230;...
    0.7170    0.1920    0.1720;...
    0.2310    0.6660    0.1960;...
    0.3680    0.1330    0.5880;...
    0.0660    0.4430    0.7450;...
    0.8660    0.3290    0.0000;...
    0.0070    0.4700    0.5010;...
    0.9130    0.3170    0.7210;...
    0.0620    0.2580    0.5010;...
    0.7170    0.1920    0.1720;...
    0.6110    0.4660    0.1250;...
    0.0070    0.3450    0.0540;...
    0.8620    0.6000    0.4230;...
    0.3720    0.1050    0.0310;...
    1.0000    0.8190    0.6190;...
    0.0070    0.3450    0.0540;...
    0.2270    0.7840    0.1920;...
    1.0000    0.8390    0.0390;...
    0.9600    0.4660    0.1600;...
    0.7520    0.2980    0.0430;...
    0.9800    0.5410    0.8310;...
    0.4900    0.6620    1.0000;...
    0.8660    0.3290    0.0000;...
    0.3290    0.7130    1.0000;...
    0.0660    0.4430    0.7450;...
    0.9960    0.5640    0.2620;...
    0.4540    0.9210    0.8540;...
    0.0000    0.6390    0.6390;...
    0.0620    0.2580    0.5010;...
    0.3290    0.7130    1.0000;...
    1.0000    0.2700    0.2270;...
    0.5640    0.1490    0.1330;...
    0.0660    0.4430    0.7450];

xVals = 1:nBins;
for iFig = 1:2

    figure
    hold on
    %this will plot everything first
    for iExpt = 1:nExpts
        [medInts, LCIints, UCIints, LQints, UQints] = makeMedianConfIntsQuartiles(allExptData{iExpt});
        patch([xVals fliplr(xVals)], ...
            [LQints fliplr(UQints)], figCols(iExpt,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], ...
            [LCIints fliplr(UCIints)], figCols(iExpt,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, medInts, '-', 'Color', figCols(iExpt,:), ...
            'MarkerFaceColor', figCols(1,:), 'LineWidth', 1.5); %median
    end
    [medInts, LCIints, UCIints, LQints, UQints] = makeMedianConfIntsQuartiles(newCombData);
    patch([xVals fliplr(xVals)], ...
        [LQints fliplr(UQints)], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([xVals fliplr(xVals)], ...
        [LCIints fliplr(UCIints)], [0 0 0], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(xVals, medInts, '-', 'Color', [0 0 0], ...
        'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5); %median
    ctrlXTicks = 1:nBins;
    lowerX = min(xVals)-0.5;
    upperX = max(xVals)+0.5;

    ctrlFig = gcf;
    ctrlax = gca;

    ylabel('Venus-MAD2 intensity, normalised (AU)');
    lowerCtrlLimit = -0.05;
    if iFig == 1
        upperCtrlLimit = max(yticks);
    else
        upperCtrlLimit = 1.2;
    end

    set(gca,'YLim', [lowerCtrlLimit upperCtrlLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);

    xlabel('MAD2 pseudo-timeline bins');
    
    ctrlRange = upperCtrlLimit - lowerCtrlLimit;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    if iFig == 1
        set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlRange*200])
    else
        set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, 275, ctrlRange*200])
    end
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
function [medianVal, lowerCI, upperCI, lowerQuart, upperQuart] = makeMedianConfIntsQuartiles(binnedMatrix)
%made to work with data with nan values
medianVal = median(binnedMatrix, 1, 'omitnan');
lowerCI = [];
upperCI = [];
lowerQuart = [];
upperQuart = [];
nCols = size(binnedMatrix,2);
for iCol = 1:nCols
    nVals = length(find(~isnan(binnedMatrix(:,iCol))));
    lPos = (nVals*0.5) - (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101 - first two 0.5 values indicate we are assessing median quantile, third 0.5 is from 1-0.5 (probability of not 0.5).
    uPos = (nVals*0.5) + (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101
    lPct = lPos/nVals;
    uPct = uPos/nVals;
    wklowerCI = quantile(binnedMatrix(:,iCol), lPct);
    wkupperCI = quantile(binnedMatrix(:,iCol), uPct);
    wklowerQuart = quantile(binnedMatrix(:,iCol), 0.25);
    wkupperQuart = quantile(binnedMatrix(:,iCol), 0.75);
    lowerCI = cat(2, lowerCI, wklowerCI);
    upperCI = cat(2, upperCI, wkupperCI);
    lowerQuart = cat(2, lowerQuart, wklowerQuart);
    upperQuart = cat(2, upperQuart, wkupperQuart);
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