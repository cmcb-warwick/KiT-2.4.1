function Cells = getPctCellStages(intiMCellStages, varargin)
% GETPCTCELLSTAGES Produces plot of KT stages in bins. Now works on paired
% data!
% 
%   GETPCTCELLSTAGES(intiMCellStages).
%   intiMCellStages must be organised as {Mad2intiM1, cellStageTable1;
%   Mad2intiM2, cellStageTable2; ...}.
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
%   swapBin: {0} or 1. Whether to have maximum intensity as maximum (0) or
%       minimum (1) bin number when plotting.
%
%   makeGraph: 0 or {1}. Whether to plot graph or not.
%
% Copyright (c) 2024 C. C. Conway

opts.nBins = 25;
opts.normalise = 1;
opts.normIndInt = 1;
opts.indMarker = 'Venus-Mad2';
opts.swapBin = 0;
opts.makeGraph = 1;
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
normInt = opts.normIndInt;
nExpts = size(intiMCellStages, 1);
allMad2 = [];
allLabels = [];
allCells = {};

for iExpt = 1:nExpts
    cellStages = intiMCellStages{iExpt, 2}; %gets table of cell stages
    intData = intiMCellStages{iExpt, 1}.intensity.mean; %gets intensity
    wkCells.Rosette = [];
    wkCells.Congressing = [];
    wkCells.LatePM = [];
    wkCells.Meta = [];
    for iCell = 1:height(cellStages)
        %if cell class is xyz, will put that cell number in xyz substructure
        cellPhase = cellStages.Phase(iCell);
        if strcmp(cellPhase, 'Rosette')
            wkCells.Rosette = cat(1, wkCells.Rosette, iCell); 
        elseif strcmp(cellPhase, 'Congressing')
            wkCells.Congressing = cat(1, wkCells.Congressing, iCell);
        elseif strcmp(cellPhase, 'Late prometa')
            wkCells.LatePM = cat(1, wkCells.LatePM, iCell);
        elseif strcmp(cellPhase, 'Metaphase')
            wkCells.Meta = cat(1, wkCells.Meta, iCell);
        end
    end
    allCells{iExpt} = wkCells;
    intData.inner = rmLowCenpC(intiMCellStages{iExpt, 1}); %20250506
    Mad2Int = intData.outer; %assumes Mad2 data is in outer
    if normInt
        Mad2Int = Mad2Int ./intData.inner; %normalises to inner
    end
    
    cellLabels = intiMCellStages{iExpt, 1}.label; %get labels for intiM
    [Mad2WkMatrix, labelWkMatrix] = makeBinnedIntLabels(Mad2Int, cellLabels, nBins, iExpt, opts.swapBin); %make binned intiM with labels also binned
    if isNorm
        %normalise to max median
        Mad2median = median(Mad2WkMatrix);
        maxMad2 = max(Mad2median);
        Mad2WkMatrix = Mad2WkMatrix/maxMad2;
    end
    %concatenate all Mad2 and all labels
    allMad2 = cat(1, allMad2, Mad2WkMatrix);
    allLabels = cat(1, allLabels, labelWkMatrix);
end

[finalMad2, finalLabels] = reBinMad2Label(allMad2, allLabels); %rebin everything including labels
if isNorm
    maxNormMad2 = max(median(finalMad2));
    finalMad2 = finalMad2/maxNormMad2;
end

%%
allCellRosette = 0;
allCellCongressing = 0;
allCellLatePM = 0;
allCellMeta = 0;
allKTRosette = zeros(1, nBins);
allKTCongressing = zeros(1, nBins);
allKTLatePM = zeros(1, nBins);
allKTMeta = zeros(1, nBins);
Cells.Stats.Bins.Rosette = nan(nExpts, nBins);
Cells.Stats.Bins.Congressing = nan(nExpts, nBins);
Cells.Stats.Bins.LatePM = nan(nExpts, nBins);
Cells.Stats.Bins.Meta = nan(nExpts, nBins);
Cells.Ints.Bins.PerExpt = cell(nExpts, nBins);


allIntRosette = [];
allIntCongressing = [];
allIntLatePM = [];
allIntMeta = [];
allNCells = 0;
for iExpt = 1:nExpts
    exptIntRosette = [];
    exptIntCongressing = [];
    exptIntLatePM = [];
    exptIntMeta = [];

    exptCells = allCells{iExpt};
    nCell.Rosette = length(exptCells.Rosette);
    allCellRosette = allCellRosette + length(exptCells.Rosette);
    nCell.Congressing = length(exptCells.Congressing);
    allCellCongressing = allCellCongressing + length(exptCells.Congressing);
    nCell.LatePM = length(exptCells.LatePM);
    allCellLatePM = allCellLatePM + length(exptCells.LatePM);
    nCell.Meta = length(exptCells.Meta);
    allCellMeta = allCellMeta + length(exptCells.Meta);
    Cells.Stats.IndExpts(iExpt).nCells = nCell;
    numCells = length(exptCells.Rosette) + length(exptCells.Congressing) + length(exptCells.LatePM) + length(exptCells.Meta);
    allNCells = allNCells + numCells;
    pctCell.Rosette = length(exptCells.Rosette)/numCells;
    pctCell.Congressing = length(exptCells.Congressing)/numCells;
    pctCell.LatePM = length(exptCells.LatePM)/numCells;
    pctCell.Meta = length(exptCells.Meta)/numCells;
    Cells.Stats.IndExpts(iExpt).pctCells = pctCell;
    
    nKRosette = 0;
    nKCongressing = 0;
    nKLatePM = 0;
    nKMeta = 0;
    

    for iBin = 1:nBins
        %moved from above loop 2026.01.26
        iKRosette = [];
        iKCongressing = [];
        iKLatePM = [];
        iKMeta = [];
        %end
        nBRosette = 0;
        nBCongressing = 0;
        nBLatePM = 0;
        nBMeta = 0;
        for iKT = 1:length(finalMad2(:,iBin))
            exptID = str2num(finalLabels{iKT, iBin}(1:2));
            if exptID == iExpt
                cellID = str2num(finalLabels{iKT, iBin}(3:4));
                if ismember(cellID, exptCells.Rosette)
                    nKRosette = nKRosette + 1;
                    nBRosette = nBRosette + 1;
                    allIntRosette = vertcat(allIntRosette, finalMad2(iKT, iBin));
                    exptIntRosette = vertcat(exptIntRosette, finalMad2(iKT, iBin));
                    iKRosette = vertcat(iKRosette, finalMad2(iKT, iBin));
                elseif ismember(cellID, exptCells.Congressing)
                    nKCongressing = nKCongressing + 1;
                    nBCongressing = nBCongressing + 1;
                    allIntCongressing = vertcat(allIntCongressing, finalMad2(iKT, iBin));
                    exptIntCongressing = vertcat(exptIntCongressing, finalMad2(iKT, iBin));
                    iKCongressing = vertcat(iKCongressing, finalMad2(iKT, iBin));
                elseif ismember(cellID, exptCells.LatePM)
                    nKLatePM = nKLatePM + 1;
                    nBLatePM = nBLatePM + 1;
                    allIntLatePM = vertcat(allIntLatePM, finalMad2(iKT, iBin));
                    exptIntLatePM = vertcat(exptIntLatePM, finalMad2(iKT, iBin));
                    iKLatePM = vertcat(iKLatePM, finalMad2(iKT, iBin));
                elseif ismember(cellID, exptCells.Meta)
                    nKMeta = nKMeta + 1;
                    nBMeta = nBMeta + 1;
                    allIntMeta = vertcat(allIntMeta, finalMad2(iKT, iBin));
                    exptIntMeta = vertcat(exptIntMeta, finalMad2(iKT, iBin));
                    iKMeta = vertcat(iKMeta, finalMad2(iKT, iBin));
                end
            end
        end
        allKTRosette(iBin) = allKTRosette(iBin) + nBRosette;
        allKTCongressing(iBin) = allKTCongressing(iBin)+nBCongressing;
        allKTLatePM(iBin) = allKTLatePM(iBin)+nBLatePM;
        allKTMeta(iBin) = allKTMeta(iBin)+nBMeta;
        Cells.Stats.Bins.Rosette(iExpt, iBin) = nBRosette;
        Cells.Stats.Bins.Congressing(iExpt, iBin) = nBCongressing;
        Cells.Stats.Bins.LatePM(iExpt, iBin) = nBLatePM;
        Cells.Stats.Bins.Meta(iExpt, iBin) = nBMeta;
        
        intStruct.Rosette = iKRosette;
        intStruct.Congressing = iKCongressing;
        intStruct.LatePM = iKLatePM;
        intStruct.Meta = iKMeta;
        Cells.Ints.Bins.PerExpt{iExpt, iBin} = intStruct;

    end
    
    nKT.Rosette = nKRosette;
    nKT.Congressing = nKCongressing;
    nKT.LatePM = nKLatePM;
    nKT.Meta = nKMeta;
    Cells.Stats.IndExpts(iExpt).nKTs = nKT;

    sumKTs = nKRosette + nKCongressing + nKLatePM + nKMeta;
    pctKT.Rosette = nKRosette/sumKTs;
    pctKT.Congressing = nKCongressing/sumKTs;
    pctKT.LatePM = nKLatePM/sumKTs;
    pctKT.Meta = nKMeta/sumKTs;
    Cells.Stats.IndExpts(iExpt).pctKTs = pctKT;

    exptMedRosette = median(exptIntRosette);
    exptMedCongressing = median(exptIntCongressing);
    exptMedLatePM = median(exptIntLatePM);
    exptMedMeta = median(exptIntMeta);
    
    Cells.Ints.Expts(iExpt).medRosette = exptMedRosette;
    Cells.Ints.Expts(iExpt).medCongressing = exptMedCongressing;
    Cells.Ints.Expts(iExpt).medLatePM = exptMedLatePM;
    Cells.Ints.Expts(iExpt).medMeta = exptMedMeta;

    Cells.Ints.Expts(iExpt).allRosette = exptIntRosette;
    Cells.Ints.Expts(iExpt).allCongressing = exptIntCongressing;
    Cells.Ints.Expts(iExpt).allLatePM = exptIntLatePM;
    Cells.Ints.Expts(iExpt).allMeta = exptIntMeta;
end
Cells.Ints.Stages.Rosette = allIntRosette;
Cells.Ints.Stages.Congressing = allIntCongressing;
Cells.Ints.Stages.LatePM = allIntLatePM;
Cells.Ints.Stages.Meta = allIntMeta;
Cells.Ints.Bins.AllInts = finalMad2;
nKTsArray = [allKTRosette; allKTCongressing; allKTLatePM; allKTMeta];

Cells.Stats.Bins.Summary.nKTs.Table = array2table(nKTsArray, 'RowNames',{'Rosette', 'Congressing', 'LatePM', 'Meta'});
Cells.Stats.Bins.Summary.nKTs.Array = nKTsArray;
Cells.Stats.nKTs = numel(finalMad2);
Cells.Stats.nCells = allNCells;
Cells.Stats.nExpts = nExpts;
binSz = size(finalMad2, 1);
pctKTsArray = nKTsArray/binSz;

Cells.Stats.Bins.Summary.pctKTs.Table = array2table(pctKTsArray, 'RowNames',{'Rosette', 'Congressing', 'LatePM', 'Meta'});
Cells.Stats.Bins.Summary.pctKTs.Array = pctKTsArray;

%%
if opts.makeGraph
    stackedBarChart = [];
    if opts.swapBin
        for iBin = nBins:-1:1
            workingdata = [pctKTsArray(1, iBin) pctKTsArray(2, iBin) pctKTsArray(3, iBin) pctKTsArray(4, iBin)];
            stackedBarChart = cat(1, stackedBarChart, workingdata);
        end
    else
        for iBin = 1:nBins
            workingdata = [pctKTsArray(1, iBin) pctKTsArray(2, iBin) pctKTsArray(3, iBin) pctKTsArray(4, iBin)];
            stackedBarChart = cat(1, stackedBarChart, workingdata);
        end
    
    end
    
    
    %%
    fig = figure; %('Position',[300 300 900 420])
    ctrlXVals = 1:nBins;
    yyaxis right
    bar(ctrlXVals, stackedBarChart, 'stacked')
    if opts.normIndInt
        indDivSignaller = '/CenpC';
    else
        indDivSignaller = '';
    end
    ylabel('Proportion of KTs in bin from each mitotic phase')
    
    if opts.swapBin
        xlabel('MAD2 pseudo-timeline');
        set(gca, 'XDir', 'reverse', 'XLim', [0.5 nBins+0.5], 'XTick', 1:nBins, 'XTickLabel', 1:nBins, 'TickDir', 'out', 'FontSize', 9)
    else
        xlabel('MAD2 pseudo-timeline');
        set(gca, 'XLim', [0.5 nBins+0.5], 'XTick', 1:nBins, 'XTickLabel', 1:nBins, 'TickDir', 'out', 'FontSize', 9)
    end
    hold on
    [medianVal, lowerCI, upperCI, lowerQuart, upperQuart] = makeMedianConfIntsQuartiles(finalMad2);
    yyaxis left
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [lowerQuart fliplr(upperQuart)], 'k', ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [lowerCI fliplr(upperCI)], 'k', ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, medianVal, '-', 'Color', 'k', ...
        'MarkerFaceColor', 'k', 'LineWidth', 1.5);
    ylabel(sprintf('Venus-Mad2/%d intensity, normalised (AU)', indDivSignaller))
    
    
    ctrlAx = gca;
    set(ctrlAx, 'Units', 'pixels');
    ctrlPos = get(ctrlAx, 'Position');
    yyaxis left
    set(ctrlAx, 'ylim', [-0.05 1.2]);
    yyaxis right
    set(ctrlAx, 'ylim', [0 1])
    set(ctrlAx, 'Position', [ctrlPos(1), ctrlPos(2), 275, 1.25*200]) %this is to make all x ax size same and to scale up so all are 1.2 high
    set(fig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (1.25*200)+(ctrlPos(2)*2)])
    set(ctrlAx, 'XTickLabelRotation', 0)
end

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


%% makeBinnedIntLabels EDITED 20251031
function [IndIntsBinned, LabelsBinned] = makeBinnedIntLabels(BinByIntiM, intiMlabels, numBins, iExpt, swapBin)
%makes a binned intensity matrix with an equivalent label matrix
nLabels = length(intiMlabels);
for iLabel = 1:nLabels
    %for each KT (/KT pair), updates that label to be cell 1 etc
    intiMlabels(iLabel, 1:2) = num2str(iExpt, '%02d');
end

%despite my fear, the below duplicates the labels as ROWS, not columns.
nCols = size(BinByIntiM, 2);
if nCols > 1
    intiMlabels = repmat(intiMlabels, nCols, 1);
end
indIntsCol = BinByIntiM(:);

labelsCol = cellstr(intiMlabels);

Total = size(find(~isnan(indIntsCol)),1); %all numbers
n = floor(Total/numBins);
k = n;
bin = 0;
IndIntsBinned = nan(k, numBins);
LabelsBinned = {};


while n < (Total + 1)
    bin = bin + 1;
    
    [maxIndSubset, maxIdx] = maxk(indIntsCol, n); %gets highest n values and their positions (ordered highest to lowest)
    [minIndSubset, minIdx] = mink(maxIndSubset, k); %gets lowest k values from highest n values and their positions (ordered lowest to highest)
    
    maxLabSubset = labelsCol(maxIdx); %gets corresponding label values from highest n indInt values
    minLabSubset = maxLabSubset(minIdx); %gets corresponding label values from lowest k of highest n indInt values


    IndIntsBinned(:, bin) = minIndSubset;
    LabelsBinned = cat(2, LabelsBinned, minLabSubset);
    
    n = n + k;
end

if swapBin
    IndIntsBinned = fliplr(IndIntsBinned); %if first column is min, then top left should be min value
    LabelsBinned = fliplr(LabelsBinned); %flips correspondingly

else
    IndIntsBinned = flipud(IndIntsBinned); %if last column is min, flip matrix accordingly
    LabelsBinned = flipud(LabelsBinned); %flips correspondingly
end


% the following lines were lifted from binIntiMs_labels.
%[IndIntsWithNAN, indSortidxWithNAN] = sort(indIntsCol, 'descend');
%LabelsWithNAN = labelsCol(indSortidxWithNAN);

%rmPos = find(isnan(IndIntsWithNAN));
%keepAfter = max(rmPos)+1;
%IndIntsAll = IndIntsWithNAN(keepAfter:end);
%LabelsAll = LabelsWithNAN(keepAfter:end);

end


%%
function [Mad2NewBin, labelNewBin] = reBinMad2Label(Mad2Matrix, labelMatrix)
numBins = size(Mad2Matrix, 2);
binSize = size(Mad2Matrix, 1);
colInd = Mad2Matrix(:);
colDep = labelMatrix(:);
[sortInd, sortOrder] = sort(colInd, 'descend'); %20250506
sortDep = colDep(sortOrder); %sort labels by resorting order of Mad2
firstn = 1;
lastn = binSize;
Mad2NewBin = [];
labelNewBin = {}; %important to make this a cell!! otherwise labels will all merge together
for iBin = 1:numBins
    iCol = sortInd(firstn:lastn); %put sorted data in correct bin order
    dCol = sortDep(firstn:lastn); %put sorted data in correct bin order
    Mad2NewBin = cat(2, Mad2NewBin, iCol);
    labelNewBin = cat(2, labelNewBin, dCol);
    firstn = firstn + binSize;
    lastn = lastn + binSize;
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