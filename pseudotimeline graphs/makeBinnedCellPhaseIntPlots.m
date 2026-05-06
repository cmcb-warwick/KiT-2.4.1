function intData = makeBinnedCellPhaseIntPlots(intiMsCellStages, varargin)
%lots of the code is lifted from makeIntPlots
% MAKECELLPHASEINTPLOTS Produces plots of multiple kinetchore intensity
% measurements over multiple cells and experiments split into bins based on
% cell phase.
% 
%   MAKECELLPHASEINTPLOTS(intiMsCellStages).
%   intiMsCellStages must be organised as {Mad2intiM1, otherintiM1, cellPhasesTable1;
%   Mad2intiM2, otherintiM2, cellPhasesTable2; ...}.
%   No option of treatIntiMs or extraIntiMs in this function.
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
%   normDepInt: 0 or {1}. Whether to normalise the dependent intensity
%       marker to CenpC or not.
%
%   indMarker: {nan} or string of intensity marker that you binned on. Will
%       request this later if not filled in now.
%
%   depMarker: {nan} or string of second intensity marker. Will request
%       this later if not filled in now.
%
%   depInner: {0} or 1. Whether to use the inner marker of the second
%       (dependent) intiM for plotting.
%
%   minBinFirst: {0} or 1. Whether to have maximum intensity as maximum (0) or
%       minimum (1) bin number when plotting.
%
% Copyright (c) 2025 C. C. Conway
opts.nBins = 25;
opts.normalise = 1;
opts.normIndInt = 1;
opts.normDepInt = 1;
opts.indMarker = 'Venus-MAD2';
opts.depMarker = nan;
opts.depInner = 0;
opts.minBinFirst = 0;
opts.calcHalfChange = 0;
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = size(intiMsCellStages, 1);

indCtrl = [];
depCtrl = [];
allLabels = {};
allCells = {};

for iExpt = 1:nExpts
    cellStages = intiMsCellStages{iExpt, 3}; %gets table of cell stages
    wkCells.Rosette = [];
    wkCells.Congressing = [];
    wkCells.LatePM = [];
    wkCells.Metaphase = [];
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
            wkCells.Metaphase = cat(1, wkCells.Metaphase, iCell);
        end
    end
    allCells{iExpt} = wkCells;
    [innerCI, outerCI] = rmLowCenpC(intiMsCellStages{iExpt,1}); %20250429
    [innerCD, outerCD] = rmLowCenpC(intiMsCellStages{iExpt,2});
    intiMsCellStages{iExpt,1}.intensity.mean.inner = innerCI; %20250429
    intiMsCellStages{iExpt,1}.intensity.mean.outer = outerCI;
    intiMsCellStages{iExpt,2}.intensity.mean.inner = innerCD; %20250429
    intiMsCellStages{iExpt,2}.intensity.mean.outer = outerCD;

    [iCBinned, dCBinned, labelsBinned, iCAll, dCAll, labelsAll] = binIntiMs_labels(intiMsCellStages{iExpt,1}, intiMsCellStages{iExpt,2}, iExpt, 'nBins', nBins, 'IndNorm', opts.normIndInt, 'DepNorm', opts.normDepInt, 'InnerDep', opts.depInner, 'minIntFirstCol', opts.minBinFirst);
    iCMed = median(iCBinned);
    dCMed = median(dCBinned);
    
    if isNorm
        iMax = max(iCMed);
        dMax = max(dCMed);
        iCBinned = iCBinned/iMax;
        dCBinned = dCBinned/dMax;
        iCAll = iCAll/iMax;
        dCAll = dCAll/dMax;
    end
    %the following is so the whole data set is considered rather than just
    %the data points that had already been binned, potentially increasing
    %the number of data points in the final part
    %indCtrl_old = cat(1, indCtrl_old, iCBinned);
    indCtrl = cat(1, indCtrl, iCAll); %20250429
    %depCtrl_old = cat(1, depCtrl, dCBinned);
    depCtrl = cat(1, depCtrl, dCAll); %20250429

    allLabels = cat(1, allLabels, labelsAll);

end
   
[iCIntsGrouped, dCIntsGrouped, finalLabels] = rebinIntensitiesLabels(indCtrl, depCtrl, allLabels, 'nBins', nBins, 'minIntFirstCol', opts.minBinFirst);
cBinSz = size(iCIntsGrouped,1);

[indCIntsMed, indCIntsLCI, indCIntsUCI, indCIntsNeg, indCIntsPos] = makeMedianConfIntsQuartiles(iCIntsGrouped);
[depCIntsMed, depCIntsLCI, depCIntsUCI, depCIntsNeg, depCIntsPos] = makeMedianConfIntsQuartiles(dCIntsGrouped);

if isNorm
    indCMaxMed = max(indCIntsMed);
    iCIntsGrouped = iCIntsGrouped/indCMaxMed;
    indCIntsMed = indCIntsMed/indCMaxMed;
    indCIntsNeg = indCIntsNeg/indCMaxMed;
    indCIntsPos = indCIntsPos/indCMaxMed;
    indCIntsLCI = indCIntsLCI/indCMaxMed;
    indCIntsUCI = indCIntsUCI/indCMaxMed;

    depCMaxMed = max(depCIntsMed);
    dCIntsGrouped = dCIntsGrouped/depCMaxMed;
    depCIntsMed = depCIntsMed/depCMaxMed;
    depCIntsNeg = depCIntsNeg/depCMaxMed;
    depCIntsPos = depCIntsPos/depCMaxMed;
    depCIntsLCI = depCIntsLCI/depCMaxMed;
    depCIntsUCI = depCIntsUCI/depCMaxMed;
end
    
intData.all.independent.intensities = iCIntsGrouped;
CIstats.median = indCIntsMed;
CIstats.lowerCI = indCIntsLCI;
CIstats.upperCI = indCIntsUCI;
CIstats.lowerQuart = indCIntsNeg;
CIstats.upperQuart = indCIntsPos;
intData.all.independent.stats = CIstats;

intData.all.dependent.intensities = dCIntsGrouped;
CDstats.median = depCIntsMed;
CDstats.lowerCI = depCIntsLCI;
CDstats.upperCI = depCIntsUCI;
CDstats.lowerQuart = depCIntsNeg;
CDstats.upperQuart = depCIntsPos;
intData.all.dependent.stats = CDstats;

nCells_ctrl = getnCells(intiMsCellStages(:,1));
intData.all.nCells = nCells_ctrl;
intData.all.binSize = cBinSz;

%above is basically what I've done for makeIntPlots_fixedRebin

%below is now obtaining data specifically for each cell phase
%what I will do is make a copy of each lot of intensity data, then remove
%entries where they are not the desired cell phase. This will let me have
%consistent matrix sizes without having to correct for different number of
%nan values in each column.

RosetteIndData = iCIntsGrouped;
RosetteDepData = dCIntsGrouped;
RosetteLabels = finalLabels;

CongressingIndData = iCIntsGrouped;
CongressingDepData = dCIntsGrouped;
CongressingLabels = finalLabels;

LatePMIndData = iCIntsGrouped;
LatePMDepData = dCIntsGrouped;
LatePMLabels = finalLabels;

MetaphaseIndData = iCIntsGrouped;
MetaphaseDepData = dCIntsGrouped;
MetaphaseLabels = finalLabels;


for iBin = 1:nBins
    for iEntry = 1:cBinSz
        exptID = str2num(finalLabels{iEntry,iBin}(1:2));
        cellID = str2num(finalLabels{iEntry,iBin}(3:4));
        exptCells = allCells{exptID};
        if isempty(find(exptCells.Rosette == cellID))
            RosetteIndData(iEntry, iBin) = nan;
            RosetteDepData(iEntry, iBin) = nan;
            RosetteLabels{iEntry, iBin} = {''};
        end
        if isempty(find(exptCells.Congressing == cellID))
            CongressingIndData(iEntry, iBin) = nan;
            CongressingDepData(iEntry, iBin) = nan;
            CongressingLabels{iEntry, iBin} = {''};
        end
        if isempty(find(exptCells.LatePM == cellID))
            LatePMIndData(iEntry, iBin) = nan;
            LatePMDepData(iEntry, iBin) = nan;
            LatePMLabels{iEntry, iBin} = {''};
        end
        if isempty(find(exptCells.Metaphase == cellID))
            MetaphaseIndData(iEntry, iBin) = nan;
            MetaphaseDepData(iEntry, iBin) = nan;
            MetaphaseLabels{iEntry, iBin} = {''};
        end

    end
end

%% Rosette stats
[indRosetteIntsMed, indRosetteIntsLCI, indRosetteIntsUCI, indRosetteIntsNeg, indRosetteIntsPos] = makeMedianConfIntsQuartiles(RosetteIndData);
[depRosetteIntsMed, depRosetteIntsLCI, depRosetteIntsUCI, depRosetteIntsNeg, depRosetteIntsPos] = makeMedianConfIntsQuartiles(RosetteDepData);

intData.Rosette.independent.intensities = RosetteIndData;
RosetteIstats.median = indRosetteIntsMed;
RosetteIstats.lowerCI = indRosetteIntsLCI;
RosetteIstats.upperCI = indRosetteIntsUCI;
RosetteIstats.lowerQuart = indRosetteIntsNeg;
RosetteIstats.upperQuart = indRosetteIntsPos;
intData.Rosette.independent.stats = RosetteIstats;

intData.Rosette.dependent.intensities = RosetteDepData;
RosetteDstats.median = depRosetteIntsMed;
RosetteDstats.lowerCI = depRosetteIntsLCI;
RosetteDstats.upperCI = depRosetteIntsUCI;
RosetteDstats.lowerQuart = depRosetteIntsNeg;
RosetteDstats.upperQuart = depRosetteIntsPos;
intData.Rosette.dependent.stats = RosetteDstats;

%sprintf('Pause')
%pause(5)

intData.Rosette.nCells = getnCellsPhase(RosetteLabels);
intData.Rosette.binSize = getCellPhaseBinSize(RosetteIndData);

%% Congressing stats
[indCongressingIntsMed, indCongressingIntsLCI, indCongressingIntsUCI, indCongressingIntsNeg, indCongressingIntsPos] = makeMedianConfIntsQuartiles(CongressingIndData);
[depCongressingIntsMed, depCongressingIntsLCI, depCongressingIntsUCI, depCongressingIntsNeg, depCongressingIntsPos] = makeMedianConfIntsQuartiles(CongressingDepData);

intData.Congressing.independent.intensities = CongressingIndData;
CongressingIstats.median = indCongressingIntsMed;
CongressingIstats.lowerCI = indCongressingIntsLCI;
CongressingIstats.upperCI = indCongressingIntsUCI;
CongressingIstats.lowerQuart = indCongressingIntsNeg;
CongressingIstats.upperQuart = indCongressingIntsPos;
intData.Congressing.independent.stats = CongressingIstats;

intData.Congressing.dependent.intensities = CongressingDepData;
CongressingDstats.median = depCongressingIntsMed;
CongressingDstats.lowerCI = depCongressingIntsLCI;
CongressingDstats.upperCI = depCongressingIntsUCI;
CongressingDstats.lowerQuart = depCongressingIntsNeg;
CongressingDstats.upperQuart = depCongressingIntsPos;
intData.Congressing.dependent.stats = CongressingDstats;

intData.Congressing.nCells = getnCellsPhase(CongressingLabels);
intData.Congressing.binSize = getCellPhaseBinSize(CongressingIndData);

%% LatePM stats
[indLatePMIntsMed, indLatePMIntsLCI, indLatePMIntsUCI, indLatePMIntsNeg, indLatePMIntsPos] = makeMedianConfIntsQuartiles(LatePMIndData);
[depLatePMIntsMed, depLatePMIntsLCI, depLatePMIntsUCI, depLatePMIntsNeg, depLatePMIntsPos] = makeMedianConfIntsQuartiles(LatePMDepData);

intData.LatePM.independent.intensities = LatePMIndData;
LatePMIstats.median = indLatePMIntsMed;
LatePMIstats.lowerCI = indLatePMIntsLCI;
LatePMIstats.upperCI = indLatePMIntsUCI;
LatePMIstats.lowerQuart = indLatePMIntsNeg;
LatePMIstats.upperQuart = indLatePMIntsPos;
intData.LatePM.independent.stats = LatePMIstats;

intData.LatePM.dependent.intensities = LatePMDepData;
LatePMDstats.median = depLatePMIntsMed;
LatePMDstats.lowerCI = depLatePMIntsLCI;
LatePMDstats.upperCI = depLatePMIntsUCI;
LatePMDstats.lowerQuart = depLatePMIntsNeg;
LatePMDstats.upperQuart = depLatePMIntsPos;
intData.LatePM.dependent.stats = LatePMDstats;

intData.LatePM.nCells = getnCellsPhase(LatePMLabels);
intData.LatePM.binSize = getCellPhaseBinSize(LatePMIndData);

%% Metaphase stats
%sprintf('Pause')
%pause(5)
[indMetaphaseIntsMed, indMetaphaseIntsLCI, indMetaphaseIntsUCI, indMetaphaseIntsNeg, indMetaphaseIntsPos] = makeMedianConfIntsQuartiles(MetaphaseIndData);
[depMetaphaseIntsMed, depMetaphaseIntsLCI, depMetaphaseIntsUCI, depMetaphaseIntsNeg, depMetaphaseIntsPos] = makeMedianConfIntsQuartiles(MetaphaseDepData);

intData.Metaphase.independent.intensities = MetaphaseIndData;
MetaphaseIstats.median = indMetaphaseIntsMed;
MetaphaseIstats.lowerCI = indMetaphaseIntsLCI;
MetaphaseIstats.upperCI = indMetaphaseIntsUCI;
MetaphaseIstats.lowerQuart = indMetaphaseIntsNeg;
MetaphaseIstats.upperQuart = indMetaphaseIntsPos;
intData.Metaphase.independent.stats = MetaphaseIstats;

intData.Metaphase.dependent.intensities = MetaphaseDepData;
MetaphaseDstats.median = depMetaphaseIntsMed;
MetaphaseDstats.lowerCI = depMetaphaseIntsLCI;
MetaphaseDstats.upperCI = depMetaphaseIntsUCI;
MetaphaseDstats.lowerQuart = depMetaphaseIntsNeg;
MetaphaseDstats.upperQuart = depMetaphaseIntsPos;
intData.Metaphase.dependent.stats = MetaphaseDstats;

intData.Metaphase.nCells = getnCellsPhase(MetaphaseLabels);
intData.Metaphase.binSize = getCellPhaseBinSize(MetaphaseIndData);

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
%% Make figure for control
if anynan([opts.indMarker, opts.depMarker])
    labelAnswers = inputdlg({'Enter independent intensity marker:', 'Enter dependent intensity marker:'}, 'Marker proteins',[1 35], {'Mad2', 'BubR1'});
else
    labelAnswers = {opts.indMarker, opts.depMarker};
end
ctrlFig = figure;
hold on
yyaxis left
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [indCIntsNeg fliplr(indCIntsPos)], indCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [indCIntsLCI fliplr(indCIntsUCI)], indCols(1,:), ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
plot(ctrlXVals, indCIntsMed, '-', 'Color', indCols(1,:), ...
    'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);

yyaxis right
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [depCIntsNeg fliplr(depCIntsPos)], depCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [depCIntsLCI fliplr(depCIntsUCI)], depCols(1,:), ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
plot(ctrlXVals, depCIntsMed, '-', 'Color', depCols(1,:), ...
    'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);

hold off
ctrlax = gca;
yyaxis left
IndCtrlTicks = yticks;
yyaxis right
DepCtrlTicks = yticks;
upperCtrlLimit = max([IndCtrlTicks(:); DepCtrlTicks(:)]);
lowerCtrlLimit = min([IndCtrlTicks(:); DepCtrlTicks(:)]);

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(ctrlax, 'YColor', indCols(1,:), 'YLim', [lowerCtrlLimit upperCtrlLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(ctrlax, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = '/CENP-C';
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;

if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(ctrlax, 'YColor', depCols(1,:), 'YLim', [lowerCtrlLimit upperCtrlLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(ctrlax, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(ctrlax, 'XDir', 'reverse');
end
ylabel(labelDep);

xlabel(sprintf('%s%s pseudo-timeline bins', labelAnswers{1}, indDivSignaller));

%% figure for Rosette
RosetteFig = figure;
hold on
isNanPos = isnan(indRosetteIntsMed);
if any(isNanPos)
    if isNanPos(1)
        nNonNanGroups = 0; %start counting nGroups when there is the first lot of non nan data
        nNanGroups = 1;
    else
        nNonNanGroups = 1;
        nNanGroups = 0;
    end
    groupedNanXVals = []; %this is to make a dotted line on graph linking two data points where the data point between is missing
    groupedNonNanXVals = [];
    for iPos = 1:length(isNanPos)
        if ~isNanPos(iPos)
            wkNonNanData = [ctrlXVals(iPos); nNonNanGroups];
            groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData); %adds current x val to included data
            if iPos ~= length(isNanPos)
                if isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nan, then get x val halfway between
                    %nonNan and nan points as end point for line.
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                    
                    %also starts dotted line for next nan group
                    nNanGroups = nNanGroups + 1;
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                end
            end
        else
            if iPos == 1
                %if start of graph, we really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            elseif iPos == length(isNanPos)
                %if end of graph, we also really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            end
            if iPos ~= length(isNanPos)
                if ~isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nonNan, then get x val halfway between
                    %nan and nonNan points as end point for dotted line.
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                    
                    %also starts dotted line for next nan group
                    nNonNanGroups = nNonNanGroups + 1;
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                end
            end
        end
    end
    medInd_nonNan = [];
    LQInd_nonNan = [];
    UQInd_nonNan = [];
    LCIInd_nonNan = [];
    UCIInd_nonNan = [];
    
    medDep_nonNan = [];
    LQDep_nonNan = [];
    UQDep_nonNan = [];
    LCIDep_nonNan = [];
    UCIDep_nonNan = [];
    
    %get data for plotting nonNan vals
    for iVal = 1:size(groupedNonNanXVals,2)
        wkXVal = groupedNonNanXVals(1,iVal); %gives x value to plot
        wholeNumberTest = mod(wkXVal,1);
        if wholeNumberTest == 0
            %testing if there is a remainder or not after dividing
            %by 1, remainder of 0 is basically saying it is a whole
            medInd_nonNan = cat(2, medInd_nonNan, [indRosetteIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan = cat(2, LQInd_nonNan, [indRosetteIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan = cat(2, UQInd_nonNan, [indRosetteIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [indRosetteIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [indRosetteIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [depRosetteIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan = cat(2, LQDep_nonNan, [depRosetteIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan = cat(2, UQDep_nonNan, [depRosetteIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [depRosetteIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [depRosetteIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        else
            if iVal == 1
                %if you have decimal value but it's the first value
                %in the dataset, use y=mx+c, calculating m and c
                %from the next two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos+1;
                    checkVal = groupedNonNanXVals(1,checkPos); %check the next position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
                
            elseif iVal == size(groupedNonNanXVals,2)
                %if you have decimal value but it's the last value
                %in the dataset, use y=mx+c, calculating m and c
                %from the previous two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos-1; %minus rather than plus!!!
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
            else
                %if in middle of vector, fair assumption that there
                %is whole number either immediately next to it or 2
                %away
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos-1; %going back first
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x1 = checkVal;
                    end
                end
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos+1; %going forward next
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x2 = checkVal;
                    end
                end
            end

            wkMedInd = getSlopeIntersectPoint(x1, indRosetteIntsMed(x1), x2, indRosetteIntsMed(x2), wkXVal);
            wkLQInd  = getSlopeIntersectPoint(x1, indRosetteIntsNeg(x1), x2, indRosetteIntsNeg(x2), wkXVal);
            wkUQInd  = getSlopeIntersectPoint(x1, indRosetteIntsPos(x1), x2, indRosetteIntsPos(x2), wkXVal);
            wkLCIInd = getSlopeIntersectPoint(x1, indRosetteIntsLCI(x1), x2, indRosetteIntsLCI(x2), wkXVal);
            wkUCIInd = getSlopeIntersectPoint(x1, indRosetteIntsUCI(x1), x2, indRosetteIntsUCI(x2), wkXVal);

            wkMedDep = getSlopeIntersectPoint(x1, depRosetteIntsMed(x1), x2, depRosetteIntsMed(x2), wkXVal);
            wkLQDep  = getSlopeIntersectPoint(x1, depRosetteIntsNeg(x1), x2, depRosetteIntsNeg(x2), wkXVal);
            wkUQDep  = getSlopeIntersectPoint(x1, depRosetteIntsPos(x1), x2, depRosetteIntsPos(x2), wkXVal);
            wkLCIDep = getSlopeIntersectPoint(x1, depRosetteIntsLCI(x1), x2, depRosetteIntsLCI(x2), wkXVal);
            wkUCIDep = getSlopeIntersectPoint(x1, depRosetteIntsUCI(x1), x2, depRosetteIntsUCI(x2), wkXVal);

            medInd_nonNan = cat(2, medInd_nonNan, [wkMedInd; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan  = cat(2, LQInd_nonNan,  [wkLQInd;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan  = cat(2, UQInd_nonNan,  [wkUQInd;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [wkLCIInd; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [wkUCIInd; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [wkMedDep; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan  = cat(2, LQDep_nonNan,  [wkLQDep;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan  = cat(2, UQDep_nonNan,  [wkUQDep;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [wkLCIDep; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [wkUCIDep; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        end
    end

    % now do the same for the nan values
    % first set up nan median intensity matrices
    medInd_nan = [];
    medDep_nan = [];
    
    for iVal = 1:size(groupedNanXVals,2)
        wkXVal = groupedNanXVals(1,iVal); %gives x value to plot
        if wkXVal == 1
            %if you have a nan value at the start of your dataset then the
            %first 2 nonNan vals should work perfectly well to get the
            %datapoint you need (you've already calculated slope etc)
            x1 = groupedNonNanXVals(1,1);
            x2 = groupedNonNanXVals(1,2);
            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,1), x2, medInd_nonNan(1,2), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,1), x2, medDep_nonNan(1,2), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);

        elseif wkXVal == max(ctrlXVals)
            %if you have a nan value at the end of your dataset then the
            %last 2 nonNan vals should be sufficient to get datapoint you
            %need as you've already calculated slope
            nNonNanPos = size(groupedNonNanXVals,2);
            x1 = groupedNonNanXVals(1,nNonNanPos);
            x2 = groupedNonNanXVals(1,nNonNanPos-1);

            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,nNonNanPos), x2, medInd_nonNan(1,nNonNanPos-1), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,nNonNanPos), x2, medDep_nonNan(1,nNonNanPos-1), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        else
            %in theory all other numbers should be decimals... and also
            %share identity with the decimal values of the nonnan data. So
            %I can utilise my work from above to sort me out!!

            %check all vals, if the xVal matches then just take that value
            for iPos = 1:size(groupedNonNanXVals,2)
                wkNonNanVal = groupedNonNanXVals(1,iPos);
                if wkNonNanVal == wkXVal
                    wkMedInd = medInd_nonNan(1,iPos);
                    wkMedDep = medDep_nonNan(1,iPos);
                end
            end
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        end

    end
    
    %now to duplicate everything so that I can plot somewhat easily

    repNonNanXVals = repmat(groupedNonNanXVals, nNonNanGroups, 1);
    
    medInd_nonNan = repmat(medInd_nonNan, nNonNanGroups, 1);
    LQInd_nonNan  = repmat(LQInd_nonNan, nNonNanGroups, 1);
    UQInd_nonNan  = repmat(UQInd_nonNan, nNonNanGroups, 1);
    LCIInd_nonNan = repmat(LCIInd_nonNan, nNonNanGroups, 1);
    UCIInd_nonNan = repmat(UCIInd_nonNan, nNonNanGroups, 1);
    
    medDep_nonNan = repmat(medDep_nonNan, nNonNanGroups, 1);
    LQDep_nonNan  = repmat(LQDep_nonNan, nNonNanGroups, 1);
    UQDep_nonNan  = repmat(UQDep_nonNan, nNonNanGroups, 1);
    LCIDep_nonNan = repmat(LCIDep_nonNan, nNonNanGroups, 1);
    UCIDep_nonNan = repmat(UCIDep_nonNan, nNonNanGroups, 1);
        
    for iGroup = 1:nNonNanGroups
        for iVal = 1:size(groupedNonNanXVals,2)
            if repNonNanXVals(iGroup*2,iVal) ~= iGroup
                repNonNanXVals((iGroup*2)-1,iVal) = nan;

                medInd_nonNan((iGroup*2)-1,iVal) = nan;
                LQInd_nonNan((iGroup*2)-1,iVal) = nan;
                UQInd_nonNan((iGroup*2)-1,iVal) = nan;
                LCIInd_nonNan((iGroup*2)-1,iVal) = nan;
                UCIInd_nonNan((iGroup*2)-1,iVal) = nan;

                medDep_nonNan((iGroup*2)-1,iVal) = nan;
                LQDep_nonNan((iGroup*2)-1,iVal) = nan;
                UQDep_nonNan((iGroup*2)-1,iVal) = nan;
                LCIDep_nonNan((iGroup*2)-1,iVal) = nan;
                UCIDep_nonNan((iGroup*2)-1,iVal) = nan;
            end
        end
    end

    repXValsNoData = repmat(groupedNanXVals, nNanGroups, 1);
    medInd_nan = repmat(medInd_nan, nNanGroups, 1);
    medDep_nan = repmat(medDep_nan, nNanGroups, 1);

    for iGroup = 1:nNanGroups
        for iVal = 1:size(groupedNanXVals,2)
            if repXValsNoData(iGroup*2, iVal) ~= iGroup
                repXValsNoData((iGroup*2)-1,iVal) = nan;
                medInd_nan((iGroup*2)-1,iVal) = nan;
                medDep_nan((iGroup*2)-1,iVal) = nan;
            end
        end
    end


    for iGroup = 1:nNanGroups
    
        xVals = repXValsNoData((iGroup*2)-1,:);        xVals = rmmissing(xVals);
        indMed = medInd_nan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        depMed = medDep_nan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
    
        yyaxis left
        plot(xVals, indMed, ':', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
        yyaxis right
        plot(xVals, depMed, ':', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    end
    
    for iGroup = 1:nNonNanGroups
    
        xVals = repNonNanXVals((iGroup*2)-1,:);      xVals = rmmissing(xVals);
    
        indMed = medInd_nonNan((iGroup*2)-1,:);      indMed = rmmissing(indMed);
        indLQ = LQInd_nonNan((iGroup*2)-1,:);        indLQ = rmmissing(indLQ);
        indUQ = UQInd_nonNan((iGroup*2)-1,:);        indUQ = rmmissing(indUQ);
        indLCI = LCIInd_nonNan((iGroup*2)-1,:);      indLCI = rmmissing(indLCI);
        indUCI = UCIInd_nonNan((iGroup*2)-1,:);      indUCI = rmmissing(indUCI);

        depMed = medDep_nonNan((iGroup*2)-1,:);      depMed = rmmissing(depMed);
        depLQ = LQDep_nonNan((iGroup*2)-1,:);        depLQ = rmmissing(depLQ);
        depUQ = UQDep_nonNan((iGroup*2)-1,:);        depUQ = rmmissing(depUQ);
        depLCI = LCIDep_nonNan((iGroup*2)-1,:);      depLCI = rmmissing(depLCI);
        depUCI = UCIDep_nonNan((iGroup*2)-1,:);      depUCI = rmmissing(depUCI);
    
        yyaxis left
        patch([xVals fliplr(xVals)], [indLQ fliplr(indUQ)], indCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [indLCI fliplr(indUCI)], indCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, indMed, '-', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
        
        yyaxis right
        patch([xVals fliplr(xVals)], [depLQ fliplr(depUQ)], depCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [depLCI fliplr(depUCI)], depCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, depMed, '-', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    
    end
        
    
else
    yyaxis left
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indRosetteIntsNeg fliplr(indRosetteIntsPos)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indRosetteIntsLCI fliplr(indRosetteIntsUCI)], indCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, indRosetteIntsMed, '-', 'Color', indCols(1,:), ...
        'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
    yyaxis right
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depRosetteIntsNeg fliplr(depRosetteIntsPos)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depRosetteIntsLCI fliplr(depRosetteIntsUCI)], depCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, depRosetteIntsMed, '-', 'Color', depCols(1,:), ...
        'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
end


hold off
Rosetteax = gca;
yyaxis left
RosetteIndTicks = yticks;
yyaxis right
RosetteDepTicks = yticks;
RosetteUpperLimit = max([RosetteIndTicks(:); RosetteDepTicks(:)]);
RosetteLowerLimit = min([RosetteIndTicks(:); RosetteDepTicks(:)]);

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(Rosetteax, 'YColor', indCols(1,:), 'YLim', [RosetteLowerLimit RosetteUpperLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(Rosetteax, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = '/CENP-C';
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;


if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(Rosetteax, 'YColor', depCols(1,:), 'YLim', [RosetteLowerLimit RosetteUpperLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(Rosetteax, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(Rosetteax, 'XDir', 'reverse');
end
ylabel(labelDep);

xlabel(sprintf('%s%s PT bins, Rosette', labelAnswers{1}, indDivSignaller));

%% figure for Congressing
CongressingFig = figure;
hold on
isNanPos = isnan(indCongressingIntsMed);
if any(isNanPos)
    if isNanPos(1)
        nNonNanGroups = 0; %start counting nGroups when there is the first lot of non nan data
        nNanGroups = 1;
    else
        nNonNanGroups = 1;
        nNanGroups = 0;
    end
    groupedNanXVals = []; %this is to make a dotted line on graph linking two data points where the data point between is missing
    groupedNonNanXVals = [];
    for iPos = 1:length(isNanPos)
        if ~isNanPos(iPos)
            wkNonNanData = [ctrlXVals(iPos); nNonNanGroups];
            groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData); %adds current x val to included data
            if iPos ~= length(isNanPos)
                if isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nan, then get x val halfway between
                    %nonNan and nan points as end point for line.
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                    
                    %also starts dotted line for next nan group
                    nNanGroups = nNanGroups + 1;
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                end
            end
        else
            if iPos == 1
                %if start of graph, we really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            elseif iPos == length(isNanPos)
                %if end of graph, we also really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            end
            if iPos ~= length(isNanPos)
                if ~isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nonNan, then get x val halfway between
                    %nan and nonNan points as end point for dotted line.
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                    
                    %also starts dotted line for next nan group
                    nNonNanGroups = nNonNanGroups + 1;
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                end
            end
        end
    end
    medInd_nonNan = [];
    LQInd_nonNan = [];
    UQInd_nonNan = [];
    LCIInd_nonNan = [];
    UCIInd_nonNan = [];
    
    medDep_nonNan = [];
    LQDep_nonNan = [];
    UQDep_nonNan = [];
    LCIDep_nonNan = [];
    UCIDep_nonNan = [];
    
    %get data for plotting nonNan vals
    for iVal = 1:size(groupedNonNanXVals,2)
        wkXVal = groupedNonNanXVals(1,iVal); %gives x value to plot
        wholeNumberTest = mod(wkXVal,1);
        if wholeNumberTest == 0
            %testing if there is a remainder or not after dividing
            %by 1, remainder of 0 is basically saying it is a whole
            medInd_nonNan = cat(2, medInd_nonNan, [indCongressingIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan = cat(2, LQInd_nonNan, [indCongressingIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan = cat(2, UQInd_nonNan, [indCongressingIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [indCongressingIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [indCongressingIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [depCongressingIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan = cat(2, LQDep_nonNan, [depCongressingIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan = cat(2, UQDep_nonNan, [depCongressingIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [depCongressingIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [depCongressingIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        else
            if iVal == 1
                %if you have decimal value but it's the first value
                %in the dataset, use y=mx+c, calculating m and c
                %from the next two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos+1;
                    checkVal = groupedNonNanXVals(1,checkPos); %check the next position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
                
            elseif iVal == size(groupedNonNanXVals,2)
                %if you have decimal value but it's the last value
                %in the dataset, use y=mx+c, calculating m and c
                %from the previous two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos-1; %minus rather than plus!!!
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
            else
                %if in middle of vector, fair assumption that there
                %is whole number either immediately next to it or 2
                %away
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos-1; %going back first
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x1 = checkVal;
                    end
                end
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos+1; %going forward next
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x2 = checkVal;
                    end
                end
            end

            wkMedInd = getSlopeIntersectPoint(x1, indCongressingIntsMed(x1), x2, indCongressingIntsMed(x2), wkXVal);
            wkLQInd  = getSlopeIntersectPoint(x1, indCongressingIntsNeg(x1), x2, indCongressingIntsNeg(x2), wkXVal);
            wkUQInd  = getSlopeIntersectPoint(x1, indCongressingIntsPos(x1), x2, indCongressingIntsPos(x2), wkXVal);
            wkLCIInd = getSlopeIntersectPoint(x1, indCongressingIntsLCI(x1), x2, indCongressingIntsLCI(x2), wkXVal);
            wkUCIInd = getSlopeIntersectPoint(x1, indCongressingIntsUCI(x1), x2, indCongressingIntsUCI(x2), wkXVal);

            wkMedDep = getSlopeIntersectPoint(x1, depCongressingIntsMed(x1), x2, depCongressingIntsMed(x2), wkXVal);
            wkLQDep  = getSlopeIntersectPoint(x1, depCongressingIntsNeg(x1), x2, depCongressingIntsNeg(x2), wkXVal);
            wkUQDep  = getSlopeIntersectPoint(x1, depCongressingIntsPos(x1), x2, depCongressingIntsPos(x2), wkXVal);
            wkLCIDep = getSlopeIntersectPoint(x1, depCongressingIntsLCI(x1), x2, depCongressingIntsLCI(x2), wkXVal);
            wkUCIDep = getSlopeIntersectPoint(x1, depCongressingIntsUCI(x1), x2, depCongressingIntsUCI(x2), wkXVal);

            medInd_nonNan = cat(2, medInd_nonNan, [wkMedInd; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan  = cat(2, LQInd_nonNan,  [wkLQInd;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan  = cat(2, UQInd_nonNan,  [wkUQInd;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [wkLCIInd; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [wkUCIInd; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [wkMedDep; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan  = cat(2, LQDep_nonNan,  [wkLQDep;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan  = cat(2, UQDep_nonNan,  [wkUQDep;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [wkLCIDep; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [wkUCIDep; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        end
    end

    % now do the same for the nan values
    % first set up nan median intensity matrices
    medInd_nan = [];
    medDep_nan = [];
    
    for iVal = 1:size(groupedNanXVals,2)
        wkXVal = groupedNanXVals(1,iVal); %gives x value to plot
        if wkXVal == 1
            %if you have a nan value at the start of your dataset then the
            %first 2 nonNan vals should work perfectly well to get the
            %datapoint you need (you've already calculated slope etc)
            x1 = groupedNonNanXVals(1,1);
            x2 = groupedNonNanXVals(1,2);
            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,1), x2, medInd_nonNan(1,2), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,1), x2, medDep_nonNan(1,2), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);

        elseif wkXVal == max(ctrlXVals)
            %if you have a nan value at the end of your dataset then the
            %last 2 nonNan vals should be sufficient to get datapoint you
            %need as you've already calculated slope
            nNonNanPos = size(groupedNonNanXVals,2);
            x1 = groupedNonNanXVals(1,nNonNanPos);
            x2 = groupedNonNanXVals(1,nNonNanPos-1);

            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,nNonNanPos), x2, medInd_nonNan(1,nNonNanPos-1), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,nNonNanPos), x2, medDep_nonNan(1,nNonNanPos-1), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        else
            %in theory all other numbers should be decimals... and also
            %share identity with the decimal values of the nonnan data. So
            %I can utilise my work from above to sort me out!!

            %check all vals, if the xVal matches then just take that value
            for iPos = 1:size(groupedNonNanXVals,2)
                wkNonNanVal = groupedNonNanXVals(1,iPos);
                if wkNonNanVal == wkXVal
                    wkMedInd = medInd_nonNan(1,iPos);
                    wkMedDep = medDep_nonNan(1,iPos);
                end
            end
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        end

    end
    
    %now to duplicate everything so that I can plot somewhat easily

    repNonNanXVals = repmat(groupedNonNanXVals, nNonNanGroups, 1);
    
    medInd_nonNan = repmat(medInd_nonNan, nNonNanGroups, 1);
    LQInd_nonNan  = repmat(LQInd_nonNan, nNonNanGroups, 1);
    UQInd_nonNan  = repmat(UQInd_nonNan, nNonNanGroups, 1);
    LCIInd_nonNan = repmat(LCIInd_nonNan, nNonNanGroups, 1);
    UCIInd_nonNan = repmat(UCIInd_nonNan, nNonNanGroups, 1);
    
    medDep_nonNan = repmat(medDep_nonNan, nNonNanGroups, 1);
    LQDep_nonNan  = repmat(LQDep_nonNan, nNonNanGroups, 1);
    UQDep_nonNan  = repmat(UQDep_nonNan, nNonNanGroups, 1);
    LCIDep_nonNan = repmat(LCIDep_nonNan, nNonNanGroups, 1);
    UCIDep_nonNan = repmat(UCIDep_nonNan, nNonNanGroups, 1);
        
    for iGroup = 1:nNonNanGroups
        for iVal = 1:size(groupedNonNanXVals,2)
            if repNonNanXVals(iGroup*2,iVal) ~= iGroup
                repNonNanXVals((iGroup*2)-1,iVal) = nan;

                medInd_nonNan((iGroup*2)-1,iVal) = nan;
                LQInd_nonNan((iGroup*2)-1,iVal) = nan;
                UQInd_nonNan((iGroup*2)-1,iVal) = nan;
                LCIInd_nonNan((iGroup*2)-1,iVal) = nan;
                UCIInd_nonNan((iGroup*2)-1,iVal) = nan;

                medDep_nonNan((iGroup*2)-1,iVal) = nan;
                LQDep_nonNan((iGroup*2)-1,iVal) = nan;
                UQDep_nonNan((iGroup*2)-1,iVal) = nan;
                LCIDep_nonNan((iGroup*2)-1,iVal) = nan;
                UCIDep_nonNan((iGroup*2)-1,iVal) = nan;
            end
        end
    end

    repXValsNoData = repmat(groupedNanXVals, nNanGroups, 1);
    medInd_nan = repmat(medInd_nan, nNanGroups, 1);
    medDep_nan = repmat(medDep_nan, nNanGroups, 1);

    for iGroup = 1:nNanGroups
        for iVal = 1:size(groupedNanXVals,2)
            if repXValsNoData(iGroup*2, iVal) ~= iGroup
                repXValsNoData((iGroup*2)-1,iVal) = nan;
                medInd_nan((iGroup*2)-1,iVal) = nan;
                medDep_nan((iGroup*2)-1,iVal) = nan;
            end
        end
    end


    for iGroup = 1:nNanGroups
    
        xVals = repXValsNoData((iGroup*2)-1,:);        xVals = rmmissing(xVals);
        indMed = medInd_nan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        depMed = medDep_nan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
    
        yyaxis left
        plot(xVals, indMed, ':', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
        yyaxis right
        plot(xVals, depMed, ':', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    end
    
    for iGroup = 1:nNonNanGroups
    
        xVals = repNonNanXVals((iGroup*2)-1,:);        xVals = rmmissing(xVals);
    
        indMed = medInd_nonNan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        indLQ = LQInd_nonNan((iGroup*2)-1,:);        indLQ = rmmissing(indLQ);
        indUQ = UQInd_nonNan((iGroup*2)-1,:);        indUQ = rmmissing(indUQ);
        indLCI = LCIInd_nonNan((iGroup*2)-1,:);        indLCI = rmmissing(indLCI);
        indUCI = UCIInd_nonNan((iGroup*2)-1,:);        indUCI = rmmissing(indUCI);

        depMed = medDep_nonNan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
        depLQ = LQDep_nonNan((iGroup*2)-1,:);        depLQ = rmmissing(depLQ);
        depUQ = UQDep_nonNan((iGroup*2)-1,:);        depUQ = rmmissing(depUQ);
        depLCI = LCIDep_nonNan((iGroup*2)-1,:);        depLCI = rmmissing(depLCI);
        depUCI = UCIDep_nonNan((iGroup*2)-1,:);        depUCI = rmmissing(depUCI);
    
        yyaxis left
        patch([xVals fliplr(xVals)], [indLQ fliplr(indUQ)], indCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [indLCI fliplr(indUCI)], indCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, indMed, '-', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
        
        yyaxis right
        patch([xVals fliplr(xVals)], [depLQ fliplr(depUQ)], depCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [depLCI fliplr(depUCI)], depCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, depMed, '-', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    
    end
        
    
else
    yyaxis left
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indCongressingIntsNeg fliplr(indCongressingIntsPos)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indCongressingIntsLCI fliplr(indCongressingIntsUCI)], indCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, indCongressingIntsMed, '-', 'Color', indCols(1,:), ...
        'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
    yyaxis right
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depCongressingIntsNeg fliplr(depCongressingIntsPos)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depCongressingIntsLCI fliplr(depCongressingIntsUCI)], depCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, depCongressingIntsMed, '-', 'Color', depCols(1,:), ...
        'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
end


hold off
Congressingax = gca;
yyaxis left
CongressingIndTicks = yticks;
yyaxis right
CongressingDepTicks = yticks;
CongressingUpperLimit = max([CongressingIndTicks(:); CongressingDepTicks(:)]);
CongressingLowerLimit = min([CongressingIndTicks(:); CongressingDepTicks(:)]);

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(Congressingax, 'YColor', indCols(1,:), 'YLim', [CongressingLowerLimit CongressingUpperLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(Congressingax, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = '/CENP-C';
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;


if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(Congressingax, 'YColor', depCols(1,:), 'YLim', [CongressingLowerLimit CongressingUpperLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(Congressingax, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(Congressingax, 'XDir', 'reverse');
end
ylabel(labelDep);

xlabel(sprintf('%s%s PT bins, Congressing', labelAnswers{1}, indDivSignaller));

%% figure for LatePM
LatePMFig = figure;
hold on
isNanPos = isnan(indLatePMIntsMed);
if any(isNanPos)
    if isNanPos(1)
        nNonNanGroups = 0; %start counting nGroups when there is the first lot of non nan data
        nNanGroups = 1;
    else
        nNonNanGroups = 1;
        nNanGroups = 0;
    end
    groupedNanXVals = []; %this is to make a dotted line on graph linking two data points where the data point between is missing
    groupedNonNanXVals = [];
    for iPos = 1:length(isNanPos)
        if ~isNanPos(iPos)
            wkNonNanData = [ctrlXVals(iPos); nNonNanGroups];
            groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData); %adds current x val to included data
            if iPos ~= length(isNanPos)
                if isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nan, then get x val halfway between
                    %nonNan and nan points as end point for line.
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                    
                    %also starts dotted line for next nan group
                    nNanGroups = nNanGroups + 1;
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                end
            end
        else
            if iPos == 1
                %if start of graph, we really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            elseif iPos == length(isNanPos)
                %if end of graph, we also really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            end
            if iPos ~= length(isNanPos)
                if ~isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nonNan, then get x val halfway between
                    %nan and nonNan points as end point for dotted line.
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                    
                    %also starts dotted line for next nan group
                    nNonNanGroups = nNonNanGroups + 1;
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                end
            end
        end
    end
    medInd_nonNan = [];
    LQInd_nonNan = [];
    UQInd_nonNan = [];
    LCIInd_nonNan = [];
    UCIInd_nonNan = [];
    
    medDep_nonNan = [];
    LQDep_nonNan = [];
    UQDep_nonNan = [];
    LCIDep_nonNan = [];
    UCIDep_nonNan = [];
    
    %get data for plotting nonNan vals
    for iVal = 1:size(groupedNonNanXVals,2)
        wkXVal = groupedNonNanXVals(1,iVal); %gives x value to plot
        wholeNumberTest = mod(wkXVal,1);
        if wholeNumberTest == 0
            %testing if there is a remainder or not after dividing
            %by 1, remainder of 0 is basically saying it is a whole
            medInd_nonNan = cat(2, medInd_nonNan, [indLatePMIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan = cat(2, LQInd_nonNan, [indLatePMIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan = cat(2, UQInd_nonNan, [indLatePMIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [indLatePMIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [indLatePMIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [depLatePMIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan = cat(2, LQDep_nonNan, [depLatePMIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan = cat(2, UQDep_nonNan, [depLatePMIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [depLatePMIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [depLatePMIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        else
            if iVal == 1
                %if you have decimal value but it's the first value
                %in the dataset, use y=mx+c, calculating m and c
                %from the next two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos+1;
                    checkVal = groupedNonNanXVals(1,checkPos); %check the next position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
                
            elseif iVal == size(groupedNonNanXVals,2)
                %if you have decimal value but it's the last value
                %in the dataset, use y=mx+c, calculating m and c
                %from the previous two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos-1; %minus rather than plus!!!
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
            else
                %if in middle of vector, fair assumption that there
                %is whole number either immediately next to it or 2
                %away
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos-1; %going back first
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x1 = checkVal;
                    end
                end
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos+1; %going forward next
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x2 = checkVal;
                    end
                end
            end

            wkMedInd = getSlopeIntersectPoint(x1, indLatePMIntsMed(x1), x2, indLatePMIntsMed(x2), wkXVal);
            wkLQInd  = getSlopeIntersectPoint(x1, indLatePMIntsNeg(x1), x2, indLatePMIntsNeg(x2), wkXVal);
            wkUQInd  = getSlopeIntersectPoint(x1, indLatePMIntsPos(x1), x2, indLatePMIntsPos(x2), wkXVal);
            wkLCIInd = getSlopeIntersectPoint(x1, indLatePMIntsLCI(x1), x2, indLatePMIntsLCI(x2), wkXVal);
            wkUCIInd = getSlopeIntersectPoint(x1, indLatePMIntsUCI(x1), x2, indLatePMIntsUCI(x2), wkXVal);

            wkMedDep = getSlopeIntersectPoint(x1, depLatePMIntsMed(x1), x2, depLatePMIntsMed(x2), wkXVal);
            wkLQDep  = getSlopeIntersectPoint(x1, depLatePMIntsNeg(x1), x2, depLatePMIntsNeg(x2), wkXVal);
            wkUQDep  = getSlopeIntersectPoint(x1, depLatePMIntsPos(x1), x2, depLatePMIntsPos(x2), wkXVal);
            wkLCIDep = getSlopeIntersectPoint(x1, depLatePMIntsLCI(x1), x2, depLatePMIntsLCI(x2), wkXVal);
            wkUCIDep = getSlopeIntersectPoint(x1, depLatePMIntsUCI(x1), x2, depLatePMIntsUCI(x2), wkXVal);

            medInd_nonNan = cat(2, medInd_nonNan, [wkMedInd; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan  = cat(2, LQInd_nonNan,  [wkLQInd;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan  = cat(2, UQInd_nonNan,  [wkUQInd;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [wkLCIInd; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [wkUCIInd; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [wkMedDep; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan  = cat(2, LQDep_nonNan,  [wkLQDep;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan  = cat(2, UQDep_nonNan,  [wkUQDep;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [wkLCIDep; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [wkUCIDep; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        end
    end

    % now do the same for the nan values
    % first set up nan median intensity matrices
    medInd_nan = [];
    medDep_nan = [];
    
    for iVal = 1:size(groupedNanXVals,2)
        wkXVal = groupedNanXVals(1,iVal); %gives x value to plot
        if wkXVal == 1
            %if you have a nan value at the start of your dataset then the
            %first 2 nonNan vals should work perfectly well to get the
            %datapoint you need (you've already calculated slope etc)
            x1 = groupedNonNanXVals(1,1);
            x2 = groupedNonNanXVals(1,2);
            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,1), x2, medInd_nonNan(1,2), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,1), x2, medDep_nonNan(1,2), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);

        elseif wkXVal == max(ctrlXVals)
            %if you have a nan value at the end of your dataset then the
            %last 2 nonNan vals should be sufficient to get datapoint you
            %need as you've already calculated slope
            nNonNanPos = size(groupedNonNanXVals,2);
            x1 = groupedNonNanXVals(1,nNonNanPos);
            x2 = groupedNonNanXVals(1,nNonNanPos-1);

            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,nNonNanPos), x2, medInd_nonNan(1,nNonNanPos-1), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,nNonNanPos), x2, medDep_nonNan(1,nNonNanPos-1), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        else
            %in theory all other numbers should be decimals... and also
            %share identity with the decimal values of the nonnan data. So
            %I can utilise my work from above to sort me out!!

            %check all vals, if the xVal matches then just take that value
            for iPos = 1:size(groupedNonNanXVals,2)
                wkNonNanVal = groupedNonNanXVals(1,iPos);
                if wkNonNanVal == wkXVal
                    wkMedInd = medInd_nonNan(1,iPos);
                    wkMedDep = medDep_nonNan(1,iPos);
                end
            end
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        end

    end
    
    %now to duplicate everything so that I can plot somewhat easily

    repNonNanXVals = repmat(groupedNonNanXVals, nNonNanGroups, 1);
    
    medInd_nonNan = repmat(medInd_nonNan, nNonNanGroups, 1);
    LQInd_nonNan  = repmat(LQInd_nonNan, nNonNanGroups, 1);
    UQInd_nonNan  = repmat(UQInd_nonNan, nNonNanGroups, 1);
    LCIInd_nonNan = repmat(LCIInd_nonNan, nNonNanGroups, 1);
    UCIInd_nonNan = repmat(UCIInd_nonNan, nNonNanGroups, 1);
    
    medDep_nonNan = repmat(medDep_nonNan, nNonNanGroups, 1);
    LQDep_nonNan  = repmat(LQDep_nonNan, nNonNanGroups, 1);
    UQDep_nonNan  = repmat(UQDep_nonNan, nNonNanGroups, 1);
    LCIDep_nonNan = repmat(LCIDep_nonNan, nNonNanGroups, 1);
    UCIDep_nonNan = repmat(UCIDep_nonNan, nNonNanGroups, 1);
        
    for iGroup = 1:nNonNanGroups
        for iVal = 1:size(groupedNonNanXVals,2)
            if repNonNanXVals(iGroup*2,iVal) ~= iGroup
                repNonNanXVals((iGroup*2)-1,iVal) = nan;

                medInd_nonNan((iGroup*2)-1,iVal) = nan;
                LQInd_nonNan((iGroup*2)-1,iVal) = nan;
                UQInd_nonNan((iGroup*2)-1,iVal) = nan;
                LCIInd_nonNan((iGroup*2)-1,iVal) = nan;
                UCIInd_nonNan((iGroup*2)-1,iVal) = nan;

                medDep_nonNan((iGroup*2)-1,iVal) = nan;
                LQDep_nonNan((iGroup*2)-1,iVal) = nan;
                UQDep_nonNan((iGroup*2)-1,iVal) = nan;
                LCIDep_nonNan((iGroup*2)-1,iVal) = nan;
                UCIDep_nonNan((iGroup*2)-1,iVal) = nan;
            end
        end
    end

    repXValsNoData = repmat(groupedNanXVals, nNanGroups, 1);
    medInd_nan = repmat(medInd_nan, nNanGroups, 1);
    medDep_nan = repmat(medDep_nan, nNanGroups, 1);

    for iGroup = 1:nNanGroups
        for iVal = 1:size(groupedNanXVals,2)
            if repXValsNoData(iGroup*2, iVal) ~= iGroup
                repXValsNoData((iGroup*2)-1,iVal) = nan;
                medInd_nan((iGroup*2)-1,iVal) = nan;
                medDep_nan((iGroup*2)-1,iVal) = nan;
            end
        end
    end


    for iGroup = 1:nNanGroups
    
        xVals = repXValsNoData((iGroup*2)-1,:);        xVals = rmmissing(xVals);
        indMed = medInd_nan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        depMed = medDep_nan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
    
        yyaxis left
        plot(xVals, indMed, ':', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
        yyaxis right
        plot(xVals, depMed, ':', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    end
    
    for iGroup = 1:nNonNanGroups
    
        xVals = repNonNanXVals((iGroup*2)-1,:);        xVals = rmmissing(xVals);
    
        indMed = medInd_nonNan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        indLQ = LQInd_nonNan((iGroup*2)-1,:);        indLQ = rmmissing(indLQ);
        indUQ = UQInd_nonNan((iGroup*2)-1,:);        indUQ = rmmissing(indUQ);
        indLCI = LCIInd_nonNan((iGroup*2)-1,:);        indLCI = rmmissing(indLCI);
        indUCI = UCIInd_nonNan((iGroup*2)-1,:);        indUCI = rmmissing(indUCI);

        depMed = medDep_nonNan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
        depLQ = LQDep_nonNan((iGroup*2)-1,:);        depLQ = rmmissing(depLQ);
        depUQ = UQDep_nonNan((iGroup*2)-1,:);        depUQ = rmmissing(depUQ);
        depLCI = LCIDep_nonNan((iGroup*2)-1,:);        depLCI = rmmissing(depLCI);
        depUCI = UCIDep_nonNan((iGroup*2)-1,:);        depUCI = rmmissing(depUCI);
    
        yyaxis left
        patch([xVals fliplr(xVals)], [indLQ fliplr(indUQ)], indCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [indLCI fliplr(indUCI)], indCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, indMed, '-', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
        
        yyaxis right
        patch([xVals fliplr(xVals)], [depLQ fliplr(depUQ)], depCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [depLCI fliplr(depUCI)], depCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, depMed, '-', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    
    end
        
    
else
    yyaxis left
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indLatePMIntsNeg fliplr(indLatePMIntsPos)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indLatePMIntsLCI fliplr(indLatePMIntsUCI)], indCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, indLatePMIntsMed, '-', 'Color', indCols(1,:), ...
        'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
    yyaxis right
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depLatePMIntsNeg fliplr(depLatePMIntsPos)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depLatePMIntsLCI fliplr(depLatePMIntsUCI)], depCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, depLatePMIntsMed, '-', 'Color', depCols(1,:), ...
        'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
end


hold off
LatePMax = gca;
yyaxis left
LatePMIndTicks = yticks;
yyaxis right
LatePMDepTicks = yticks;
LatePMUpperLimit = max([LatePMIndTicks(:); LatePMDepTicks(:)]);
LatePMLowerLimit = min([LatePMIndTicks(:); LatePMDepTicks(:)]);

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(LatePMax, 'YColor', indCols(1,:), 'YLim', [LatePMLowerLimit LatePMUpperLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(LatePMax, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = '/CENP-C';
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;


if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(LatePMax, 'YColor', depCols(1,:), 'YLim', [LatePMLowerLimit LatePMUpperLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(LatePMax, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(LatePMax, 'XDir', 'reverse');
end
ylabel(labelDep);

xlabel(sprintf('%s%s PT bins, LatePM', labelAnswers{1}, indDivSignaller));

%% figure for Metaphase
MetaphaseFig = figure;
hold on
isNanPos = isnan(indMetaphaseIntsMed);
if any(isNanPos)
    if isNanPos(1)
        nNonNanGroups = 0; %start counting nGroups when there is the first lot of non nan data
        nNanGroups = 1;
    else
        nNonNanGroups = 1;
        nNanGroups = 0;
    end
    groupedNanXVals = []; %this is to make a dotted line on graph linking two data points where the data point between is missing
    groupedNonNanXVals = [];
    for iPos = 1:length(isNanPos)
        if ~isNanPos(iPos)
            wkNonNanData = [ctrlXVals(iPos); nNonNanGroups];
            groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData); %adds current x val to included data
            if iPos ~= length(isNanPos)
                if isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nan, then get x val halfway between
                    %nonNan and nan points as end point for line.
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                    
                    %also starts dotted line for next nan group
                    nNanGroups = nNanGroups + 1;
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                end
            end
        else
            if iPos == 1
                %if start of graph, we really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            elseif iPos == length(isNanPos)
                %if end of graph, we also really need an x value there
                wkNanData = [ctrlXVals(iPos); nNanGroups];
                groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
            end
            if iPos ~= length(isNanPos)
                if ~isNanPos(iPos+1)
                    %excludes if iPos = nBins because there is no iPos+1

                    %if next val is nonNan, then get x val halfway between
                    %nan and nonNan points as end point for dotted line.
                    wkNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNanGroups];
                    groupedNanXVals = cat(2, groupedNanXVals, wkNanData);
                    
                    %also starts dotted line for next nan group
                    nNonNanGroups = nNonNanGroups + 1;
                    wkNonNanData = [mean([ctrlXVals(iPos), ctrlXVals(iPos+1)]); nNonNanGroups];
                    groupedNonNanXVals = cat(2, groupedNonNanXVals, wkNonNanData);
                end
            end
        end
    end
    medInd_nonNan = [];
    LQInd_nonNan = [];
    UQInd_nonNan = [];
    LCIInd_nonNan = [];
    UCIInd_nonNan = [];
    
    medDep_nonNan = [];
    LQDep_nonNan = [];
    UQDep_nonNan = [];
    LCIDep_nonNan = [];
    UCIDep_nonNan = [];
    
    %get data for plotting nonNan vals
    for iVal = 1:size(groupedNonNanXVals,2)
        wkXVal = groupedNonNanXVals(1,iVal); %gives x value to plot
        wholeNumberTest = mod(wkXVal,1);
        if wholeNumberTest == 0
            %testing if there is a remainder or not after dividing
            %by 1, remainder of 0 is basically saying it is a whole
            medInd_nonNan = cat(2, medInd_nonNan, [indMetaphaseIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan = cat(2, LQInd_nonNan, [indMetaphaseIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan = cat(2, UQInd_nonNan, [indMetaphaseIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [indMetaphaseIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [indMetaphaseIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [depMetaphaseIntsMed(wkXVal); groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan = cat(2, LQDep_nonNan, [depMetaphaseIntsNeg(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan = cat(2, UQDep_nonNan, [depMetaphaseIntsPos(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [depMetaphaseIntsLCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [depMetaphaseIntsUCI(wkXVal); groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        else
            if iVal == 1
                %if you have decimal value but it's the first value
                %in the dataset, use y=mx+c, calculating m and c
                %from the next two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos+1;
                    checkVal = groupedNonNanXVals(1,checkPos); %check the next position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
                
            elseif iVal == size(groupedNonNanXVals,2)
                %if you have decimal value but it's the last value
                %in the dataset, use y=mx+c, calculating m and c
                %from the previous two whole values
                wholeNumberCounter = 0;
                checkPos = iVal; %which positions in the vector to check
                while wholeNumberCounter ~= 2
                    checkPos = checkPos-1; %minus rather than plus!!!
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumberCounter = wholeNumberCounter+1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        if wholeNumberCounter == 1
                            x1 = checkVal;
                        elseif wholeNumberCounter == 2
                            x2 = checkVal;
                        end
                    end
                end
            else
                %if in middle of vector, fair assumption that there
                %is whole number either immediately next to it or 2
                %away
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos-1; %going back first
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x1 = checkVal;
                    end
                end
                wholeNumEncountered = 0;
                checkPos = iVal;
                while wholeNumEncountered == 0
                    checkPos = checkPos+1; %going forward next
                    checkVal = groupedNonNanXVals(1,checkPos); %check the previous position
                    if mod(checkVal,1) == 0
                        wholeNumEncountered = 1; %if whole number, store this for calculating the slope. Only for first 2 whole numbers encountered
                        x2 = checkVal;
                    end
                end
            end

            wkMedInd = getSlopeIntersectPoint(x1, indMetaphaseIntsMed(x1), x2, indMetaphaseIntsMed(x2), wkXVal);
            wkLQInd  = getSlopeIntersectPoint(x1, indMetaphaseIntsNeg(x1), x2, indMetaphaseIntsNeg(x2), wkXVal);
            wkUQInd  = getSlopeIntersectPoint(x1, indMetaphaseIntsPos(x1), x2, indMetaphaseIntsPos(x2), wkXVal);
            wkLCIInd = getSlopeIntersectPoint(x1, indMetaphaseIntsLCI(x1), x2, indMetaphaseIntsLCI(x2), wkXVal);
            wkUCIInd = getSlopeIntersectPoint(x1, indMetaphaseIntsUCI(x1), x2, indMetaphaseIntsUCI(x2), wkXVal);

            wkMedDep = getSlopeIntersectPoint(x1, depMetaphaseIntsMed(x1), x2, depMetaphaseIntsMed(x2), wkXVal);
            wkLQDep  = getSlopeIntersectPoint(x1, depMetaphaseIntsNeg(x1), x2, depMetaphaseIntsNeg(x2), wkXVal);
            wkUQDep  = getSlopeIntersectPoint(x1, depMetaphaseIntsPos(x1), x2, depMetaphaseIntsPos(x2), wkXVal);
            wkLCIDep = getSlopeIntersectPoint(x1, depMetaphaseIntsLCI(x1), x2, depMetaphaseIntsLCI(x2), wkXVal);
            wkUCIDep = getSlopeIntersectPoint(x1, depMetaphaseIntsUCI(x1), x2, depMetaphaseIntsUCI(x2), wkXVal);

            medInd_nonNan = cat(2, medInd_nonNan, [wkMedInd; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQInd_nonNan  = cat(2, LQInd_nonNan,  [wkLQInd;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQInd_nonNan  = cat(2, UQInd_nonNan,  [wkUQInd;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIInd_nonNan = cat(2, LCIInd_nonNan, [wkLCIInd; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIInd_nonNan = cat(2, UCIInd_nonNan, [wkUCIInd; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
            
            medDep_nonNan = cat(2, medDep_nonNan, [wkMedDep; groupedNonNanXVals(2,iVal)]); %gives median data and the group
            LQDep_nonNan  = cat(2, LQDep_nonNan,  [wkLQDep;  groupedNonNanXVals(2,iVal)]); %gives LQ data and the group
            UQDep_nonNan  = cat(2, UQDep_nonNan,  [wkUQDep;  groupedNonNanXVals(2,iVal)]); %gives UQ data and the group
            LCIDep_nonNan = cat(2, LCIDep_nonNan, [wkLCIDep; groupedNonNanXVals(2,iVal)]); %gives LCI data and the group
            UCIDep_nonNan = cat(2, UCIDep_nonNan, [wkUCIDep; groupedNonNanXVals(2,iVal)]); %gives UCI data and the group
        end
    end

    % now do the same for the nan values
    % first set up nan median intensity matrices
    medInd_nan = [];
    medDep_nan = [];
    
    for iVal = 1:size(groupedNanXVals,2)
        wkXVal = groupedNanXVals(1,iVal); %gives x value to plot
        if wkXVal == 1
            %if you have a nan value at the start of your dataset then the
            %first 2 nonNan vals should work perfectly well to get the
            %datapoint you need (you've already calculated slope etc)
            x1 = groupedNonNanXVals(1,1);
            x2 = groupedNonNanXVals(1,2);
            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,1), x2, medInd_nonNan(1,2), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,1), x2, medDep_nonNan(1,2), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);

        elseif wkXVal == max(ctrlXVals)
            %if you have a nan value at the end of your dataset then the
            %last 2 nonNan vals should be sufficient to get datapoint you
            %need as you've already calculated slope
            nNonNanPos = size(groupedNonNanXVals,2);
            x1 = groupedNonNanXVals(1,nNonNanPos);
            x2 = groupedNonNanXVals(1,nNonNanPos-1);

            wkMedInd = getSlopeIntersectPoint(x1, medInd_nonNan(1,nNonNanPos), x2, medInd_nonNan(1,nNonNanPos-1), wkXVal);
            wkMedDep = getSlopeIntersectPoint(x1, meDep_nonNan(1,nNonNanPos), x2, medDep_nonNan(1,nNonNanPos-1), wkXVal);
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        else
            %in theory all other numbers should be decimals... and also
            %share identity with the decimal values of the nonnan data. So
            %I can utilise my work from above to sort me out!!

            %check all vals, if the xVal matches then just take that value
            for iPos = 1:size(groupedNonNanXVals,2)
                wkNonNanVal = groupedNonNanXVals(1,iPos);
                if wkNonNanVal == wkXVal
                    wkMedInd = medInd_nonNan(1,iPos);
                    wkMedDep = medDep_nonNan(1,iPos);
                end
            end
            medInd_nan = cat(2, medInd_nan, [wkMedInd; groupedNanXVals(2,iVal)]);
            medDep_nan = cat(2, medDep_nan, [wkMedDep; groupedNanXVals(2,iVal)]);
        end

    end
    
    %now to duplicate everything so that I can plot somewhat easily

    repNonNanXVals = repmat(groupedNonNanXVals, nNonNanGroups, 1);
    
    medInd_nonNan = repmat(medInd_nonNan, nNonNanGroups, 1);
    LQInd_nonNan  = repmat(LQInd_nonNan, nNonNanGroups, 1);
    UQInd_nonNan  = repmat(UQInd_nonNan, nNonNanGroups, 1);
    LCIInd_nonNan = repmat(LCIInd_nonNan, nNonNanGroups, 1);
    UCIInd_nonNan = repmat(UCIInd_nonNan, nNonNanGroups, 1);
    
    medDep_nonNan = repmat(medDep_nonNan, nNonNanGroups, 1);
    LQDep_nonNan  = repmat(LQDep_nonNan, nNonNanGroups, 1);
    UQDep_nonNan  = repmat(UQDep_nonNan, nNonNanGroups, 1);
    LCIDep_nonNan = repmat(LCIDep_nonNan, nNonNanGroups, 1);
    UCIDep_nonNan = repmat(UCIDep_nonNan, nNonNanGroups, 1);
        
    for iGroup = 1:nNonNanGroups
        for iVal = 1:size(groupedNonNanXVals,2)
            if repNonNanXVals(iGroup*2,iVal) ~= iGroup
                repNonNanXVals((iGroup*2)-1,iVal) = nan;

                medInd_nonNan((iGroup*2)-1,iVal) = nan;
                LQInd_nonNan((iGroup*2)-1,iVal) = nan;
                UQInd_nonNan((iGroup*2)-1,iVal) = nan;
                LCIInd_nonNan((iGroup*2)-1,iVal) = nan;
                UCIInd_nonNan((iGroup*2)-1,iVal) = nan;

                medDep_nonNan((iGroup*2)-1,iVal) = nan;
                LQDep_nonNan((iGroup*2)-1,iVal) = nan;
                UQDep_nonNan((iGroup*2)-1,iVal) = nan;
                LCIDep_nonNan((iGroup*2)-1,iVal) = nan;
                UCIDep_nonNan((iGroup*2)-1,iVal) = nan;
            end
        end
    end

    repXValsNoData = repmat(groupedNanXVals, nNanGroups, 1);
    medInd_nan = repmat(medInd_nan, nNanGroups, 1);
    medDep_nan = repmat(medDep_nan, nNanGroups, 1);

    for iGroup = 1:nNanGroups
        for iVal = 1:size(groupedNanXVals,2)
            if repXValsNoData(iGroup*2, iVal) ~= iGroup
                repXValsNoData((iGroup*2)-1,iVal) = nan;
                medInd_nan((iGroup*2)-1,iVal) = nan;
                medDep_nan((iGroup*2)-1,iVal) = nan;
            end
        end
    end


    for iGroup = 1:nNanGroups
    
        xVals = repXValsNoData((iGroup*2)-1,:);        xVals = rmmissing(xVals);
        indMed = medInd_nan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        depMed = medDep_nan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
    
        yyaxis left
        plot(xVals, indMed, ':', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
        yyaxis right
        plot(xVals, depMed, ':', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    end
    
    for iGroup = 1:nNonNanGroups
    
        xVals = repNonNanXVals((iGroup*2)-1,:);        xVals = rmmissing(xVals);
    
        indMed = medInd_nonNan((iGroup*2)-1,:);        indMed = rmmissing(indMed);
        indLQ = LQInd_nonNan((iGroup*2)-1,:);        indLQ = rmmissing(indLQ);
        indUQ = UQInd_nonNan((iGroup*2)-1,:);        indUQ = rmmissing(indUQ);
        indLCI = LCIInd_nonNan((iGroup*2)-1,:);        indLCI = rmmissing(indLCI);
        indUCI = UCIInd_nonNan((iGroup*2)-1,:);        indUCI = rmmissing(indUCI);

        depMed = medDep_nonNan((iGroup*2)-1,:);        depMed = rmmissing(depMed);
        depLQ = LQDep_nonNan((iGroup*2)-1,:);        depLQ = rmmissing(depLQ);
        depUQ = UQDep_nonNan((iGroup*2)-1,:);        depUQ = rmmissing(depUQ);
        depLCI = LCIDep_nonNan((iGroup*2)-1,:);        depLCI = rmmissing(depLCI);
        depUCI = UCIDep_nonNan((iGroup*2)-1,:);        depUCI = rmmissing(depUCI);
    
        yyaxis left
        patch([xVals fliplr(xVals)], [indLQ fliplr(indUQ)], indCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [indLCI fliplr(indUCI)], indCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, indMed, '-', 'Color', indCols(1,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
        
        yyaxis right
        patch([xVals fliplr(xVals)], [depLQ fliplr(depUQ)], depCols(1,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xVals fliplr(xVals)], [depLCI fliplr(depUCI)], depCols(1,:), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
        plot(xVals, depMed, '-', 'Color', depCols(1,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
    
    end
        
    
else
    yyaxis left
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indMetaphaseIntsNeg fliplr(indMetaphaseIntsPos)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [indMetaphaseIntsLCI fliplr(indMetaphaseIntsUCI)], indCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, indMetaphaseIntsMed, '-', 'Color', indCols(1,:), ...
        'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
    
    yyaxis right
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depMetaphaseIntsNeg fliplr(depMetaphaseIntsPos)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
    patch([ctrlXVals fliplr(ctrlXVals)], ...
        [depMetaphaseIntsLCI fliplr(depMetaphaseIntsUCI)], depCols(1,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval
    plot(ctrlXVals, depMetaphaseIntsMed, '-', 'Color', depCols(1,:), ...
        'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
end
hold off
Metaphaseax = gca;
yyaxis left
MetaphaseIndTicks = yticks;
yyaxis right
MetaphaseDepTicks = yticks;
MetaphaseUpperLimit = max([MetaphaseIndTicks(:); MetaphaseDepTicks(:)]);
MetaphaseLowerLimit = min([MetaphaseIndTicks(:); MetaphaseDepTicks(:)]);

yyaxis left
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(Metaphaseax, 'YColor', indCols(1,:), 'YLim', [MetaphaseLowerLimit MetaphaseUpperLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(Metaphaseax, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = '/CENP-C';
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;


if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(Metaphaseax, 'YColor', depCols(1,:), 'YLim', [MetaphaseLowerLimit MetaphaseUpperLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(Metaphaseax, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(Metaphaseax, 'XDir', 'reverse');
end
ylabel(labelDep);

xlabel(sprintf('%s%s PT bins, Metaphase', labelAnswers{1}, indDivSignaller));
%% Resizing windows and axes
if isNorm
    minValsVec = [indCIntsNeg, depCIntsNeg,...
        indRosetteIntsNeg, depRosetteIntsNeg,...
        indCongressingIntsNeg, depCongressingIntsNeg,...
        indLatePMIntsNeg, depLatePMIntsNeg,...
        indMetaphaseIntsNeg, depMetaphaseIntsNeg];
    minValue = min(minValsVec, [], 'all');

    lowerTick = min([lowerCtrlLimit, RosetteLowerLimit, CongressingLowerLimit, LatePMLowerLimit, MetaphaseLowerLimit]);
    upperTick = max([upperCtrlLimit, RosetteUpperLimit, CongressingUpperLimit, LatePMUpperLimit, MetaphaseUpperLimit]);
    possVec = [lowerTick, -0.5, -0.2, -0.1, -0.05];
    
    possVec = sort(possVec);
    if minValue < possVec(2)
        lowerTick = possVec(1);
    elseif minValue < possVec(3)
        lowerTick = possVec(2);
    elseif minValue < possVec(4)
        lowerTick = possVec(3);
    elseif minValue < possVec(5)
        lowerTick = possVec(4);
    else
        lowerTick = possVec(5);
    end
    figure(ctrlFig);
    yyaxis left
    set(ctrlax, 'YLim', [lowerTick upperTick]);
    yyaxis right
    set(ctrlax, 'YLim', [lowerTick upperTick]);
    ctrlRange = upperTick - lowerTick;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    %%set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, 275, ctrlRange*200])
    set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlRange*200])
    set(ctrlFig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (ctrlRange*200)+(ctrlPos(2)*2)])
    set(ctrlax, 'XTickLabelRotation', 0)

    figure(RosetteFig);
    yyaxis left
    set(Rosetteax, 'YLim', [lowerTick upperTick]);
    yyaxis right
    set(Rosetteax, 'YLim', [lowerTick upperTick]);
    RosetteRange = upperTick - lowerTick;
    set(Rosetteax, 'Units', 'pixels')
    RosettePos = get(Rosetteax, 'Position');
    %%set(Rosetteax, 'Position', [RosettePos(1), RosettePos(2)*1.25, 275, RosetteRange*200])
    set(Rosetteax, 'Position', [RosettePos(1), RosettePos(2)*1.25, RosettePos(3), RosetteRange*200])
    set(RosetteFig, 'Position', [100 60 RosettePos(3)+(RosettePos(1)*2) (RosetteRange*200)+(RosettePos(2)*2)])
    set(Rosetteax, 'XTickLabelRotation', 0)

    figure(CongressingFig);
    yyaxis left
    set(Congressingax, 'YLim', [lowerTick upperTick]);
    yyaxis right
    set(Congressingax, 'YLim', [lowerTick upperTick]);
    CongressingRange = upperTick - lowerTick;
    set(Congressingax, 'Units', 'pixels')
    CongressingPos = get(Congressingax, 'Position');
    %%set(Congressingax, 'Position', [CongressingPos(1), CongressingPos(2)*1.25, 275, CongressingRange*200])
    set(Congressingax, 'Position', [CongressingPos(1), CongressingPos(2)*1.25, CongressingPos(3), CongressingRange*200])
    set(CongressingFig, 'Position', [100 60 CongressingPos(3)+(CongressingPos(1)*2) (CongressingRange*200)+(CongressingPos(2)*2)])
    set(Congressingax, 'XTickLabelRotation', 0)

    figure(LatePMFig);
    yyaxis left
    set(LatePMax, 'YLim', [lowerTick upperTick]);
    yyaxis right
    set(LatePMax, 'YLim', [lowerTick upperTick]);
    LatePMRange = upperTick - lowerTick;
    set(LatePMax, 'Units', 'pixels')
    LatePMPos = get(LatePMax, 'Position');
    %%set(LatePMax, 'Position', [LatePMPos(1), LatePMPos(2)*1.25, 275, LatePMRange*200])
    set(LatePMax, 'Position', [LatePMPos(1), LatePMPos(2)*1.25, LatePMPos(3), LatePMRange*200])
    set(LatePMFig, 'Position', [100 60 LatePMPos(3)+(LatePMPos(1)*2) (LatePMRange*200)+(LatePMPos(2)*2)])
    set(LatePMax, 'XTickLabelRotation', 0)

    figure(MetaphaseFig);
    yyaxis left
    set(Metaphaseax, 'YLim', [lowerTick upperTick]);
    yyaxis right
    set(Metaphaseax, 'YLim', [lowerTick upperTick]);
    MetaphaseRange = upperTick - lowerTick;
    set(Metaphaseax, 'Units', 'pixels')
    MetaphasePos = get(Metaphaseax, 'Position');
    %%set(Metaphaseax, 'Position', [MetaphasePos(1), MetaphasePos(2)*1.25, 275, MetaphaseRange*200])
    set(Metaphaseax, 'Position', [MetaphasePos(1), MetaphasePos(2)*1.25, MetaphasePos(3), MetaphaseRange*200])
    set(MetaphaseFig, 'Position', [100 60 MetaphasePos(3)+(MetaphasePos(1)*2) (MetaphaseRange*200)+(MetaphasePos(2)*2)])
    set(Metaphaseax, 'XTickLabelRotation', 0)

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
%% rebinIntensitiesLabels
function [rebinnedIndInts, rebinnedDepInts, rebinnedLabels] = rebinIntensitiesLabels(IndIntsOrg, DepIntsOrg, LabelsOrg, varargin)

opts.nBins = 25; %number of bins to make
opts.minIntFirstCol = 0; %whether minimum intensities are in first column (1) or last column (0)
opts = processOptions(opts,varargin{:});

%making column vectors to work with mink/maxk
indIntsCol = IndIntsOrg(:);
depIntsCol = DepIntsOrg(:);
labelsCol = LabelsOrg(:);

%setting up loops
Total = length(indIntsCol);
n = floor(Total/opts.nBins); %bin starting position
k = n; %number of KTs per bin
bin = 0;

IndIntsBinned = nan(k, opts.nBins);
DepIntsBinned = nan(k, opts.nBins);
LabelsBinned = {};

while n < (Total + 1)
    bin = bin + 1;
    
    [maxIndSubset, maxIdx] = maxk(indIntsCol, n); %gets highest n values and their positions (ordered highest to lowest)
    [minIndSubset, minIdx] = mink(maxIndSubset, k); %gets lowest k values from highest n values and their positions (ordered lowest to highest)

    maxDepSubset = depIntsCol(maxIdx); %gets corresponding depInt values from highest n indInt values
    minDepSubset = maxDepSubset(minIdx); %gets corresponding depInt values from lowest k of highest n indInt values

    maxLabSubset = labelsCol(maxIdx); %gets corresponding label values from highest n indInt values
    minLabSubset = maxLabSubset(minIdx); %gets corresponding label values from lowest k of highest n indInt values

    IndIntsBinned(:, bin) = minIndSubset;
    DepIntsBinned(:, bin) = minDepSubset;
    LabelsBinned = cat(2, LabelsBinned, minLabSubset);
    
    n = n + k;
end

if opts.minIntFirstCol
    rebinnedIndInts = fliplr(IndIntsBinned); %if first column is min, then top left should be min value
    rebinnedDepInts = fliplr(DepIntsBinned); %flips correspondingly
    rebinnedLabels = fliplr(LabelsBinned); %flips correspondingly
else
    rebinnedIndInts = flipud(IndIntsBinned); %if last column is min, flip matrix accordingly
    rebinnedDepInts = flipud(DepIntsBinned); %flips correspondingly
    rebinnedLabels = flipud(LabelsBinned); %flips correspondingly
end


end


%% makeMedianConfIntsQuartiles
function [medianVal, lowerCI, upperCI, lowerQuart, upperQuart] = makeMedianConfIntsQuartiles(binnedMatrix)
medianVal = [];
lowerCI = [];
upperCI = [];
lowerQuart = [];
upperQuart = [];
for iCol = 1:size(binnedMatrix,2)
    wkData = binnedMatrix(:, iCol);
    wkData = rmmissing(wkData);
    if isempty(wkData)
        wkMedian = nan;
        wklowerCI = nan;
        wkupperCI = nan;
        wklowerQuart = nan;
        wkupperQuart = nan;
    else
        wkMedian = median(wkData);
        nVals = size(wkData, 1);
        lPos = (nVals*0.5) - (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101 - first two 0.5 values indicate we are assessing median quantile, third 0.5 is from 1-0.5 (probability of not 0.5).
        uPos = (nVals*0.5) + (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101
        lPct = lPos/nVals;
        uPct = uPos/nVals;
        if lPct < 0
            wklowerCI = wkMedian;
        else
            wklowerCI = quantile(wkData, lPct);
        end
        if uPct > 1
            wkupperCI = wkMedian;
        else
            wkupperCI = quantile(wkData, uPct);
        end
        wklowerQuart = quantile(wkData, 0.25);
        wkupperQuart = quantile(wkData, 0.75);
    end

    medianVal = cat(2, medianVal, wkMedian);
    lowerCI = cat(2, lowerCI, wklowerCI);
    upperCI = cat(2, upperCI, wkupperCI);
    lowerQuart = cat(2, lowerQuart, wklowerQuart);
    upperQuart = cat(2, upperQuart, wkupperQuart);
end

end

%% binIntiMs_labels
function [IndIntsBinned, DepIntsBinned, LabelsBinned, IndIntsAll, DepIntsAll, LabelsAll] = binIntiMs_labels(indIntiM, depIntiM, iExpt, varargin)
% Function to enable binning of two intiMs with cell phases. The order of intiM binning
% will be shared across vectors.
% (c) C. C. Conway, Feb 2025 Edit April 2025
% Options are available:

opts.IndNorm = 0; %whether to normalise indIntiM to KT marker
opts.DepNorm = 0; %whether to normalise depIntiM to KT marker
opts.InnerInd = 0; %whether inner intensity is to be assessed for indIntiM
opts.InnerDep = 0; %whether inner intensity is to be assessed for depIntiM
opts.nBins = 25; %number of bins to make
opts.minIntFirstCol = 0; %whether minimum intensities are in first column (1) or last column (0)
opts = processOptions(opts,varargin{:});


%getting independent marker intensities and making column vector to work
%with mink/maxk
if opts.InnerInd
    if opts.IndNorm
        indInts = indIntiM.intensity.mean.inner ./ indIntiM.intensity.mean.outer;
    else
        indInts = indIntiM.intensity.mean.inner;
    end
else
    if opts.IndNorm
        indInts = indIntiM.intensity.mean.outer ./ indIntiM.intensity.mean.inner;
    else
        indInts = indIntiM.intensity.mean.outer;
    end
end
indIntsCol = indInts(:);

%getting dependent marker intensities and making them column vector
if opts.InnerDep
    if opts.DepNorm
        depInts = depIntiM.intensity.mean.inner./depIntiM.intensity.mean.outer;
    else
        depInts = depIntiM.intensity.mean.inner;
    end
else
    if opts.DepNorm
        depInts = depIntiM.intensity.mean.outer./depIntiM.intensity.mean.inner;
    else
        depInts = depIntiM.intensity.mean.outer;
    end
end
depIntsCol = depInts(:);

intiMlabels = indIntiM.label;
nLabels = length(intiMlabels);
for iLabel = 1:nLabels
    intiMlabels(iLabel, 1:2) = num2str(iExpt, '%02d');
end

%fix to permit paired data
nCols = size(indIntiM.intensity.mean.outer, 2);
if nCols > 1
    intiMlabels = repmat(intiMlabels, nCols, 1);
end

labelsCol = cellstr(intiMlabels);

%setting up loops
Total = size(find(~isnan(indIntiM.intensity.mean.outer)),1);
n = floor(Total/opts.nBins); %bin starting position
k = n; %number of KTs per bin
bin = 0;
IndIntsBinned = nan(k, opts.nBins);
DepIntsBinned = nan(k, opts.nBins);
LabelsBinned = {};

while n < (Total + 1)
    bin = bin + 1;
    
    [maxIndSubset, maxIdx] = maxk(indIntsCol, n); %gets highest n values and their positions (ordered highest to lowest)
    [minIndSubset, minIdx] = mink(maxIndSubset, k); %gets lowest k values from highest n values and their positions (ordered lowest to highest)
    
    maxDepSubset = depIntsCol(maxIdx); %gets corresponding depInt values from highest n indInt values
    minDepSubset = maxDepSubset(minIdx); %gets corresponding depInt values from lowest k of highest n indInt values
    
    maxLabSubset = labelsCol(maxIdx); %gets corresponding label values from highest n indInt values
    minLabSubset = maxLabSubset(minIdx); %gets corresponding label values from lowest k of highest n indInt values


    IndIntsBinned(:, bin) = minIndSubset;
    DepIntsBinned(:, bin) = minDepSubset;
    LabelsBinned = cat(2, LabelsBinned, minLabSubset);
    
    n = n + k;
end

if opts.minIntFirstCol
    IndIntsBinned = fliplr(IndIntsBinned); %if first column is min, then top left should be min value
    DepIntsBinned = fliplr(DepIntsBinned); %flips correspondingly
    LabelsBinned = fliplr(LabelsBinned); %flips correspondingly

else
    IndIntsBinned = flipud(IndIntsBinned); %if last column is min, flip matrix accordingly
    DepIntsBinned = flipud(DepIntsBinned); %flips correspondingly
    LabelsBinned = flipud(LabelsBinned); %flips correspondingly
end

[IndIntsWithNAN, indSortidxWithNAN] = sort(indIntsCol, 'descend');
DepIntsWithNAN = depIntsCol(indSortidxWithNAN);
LabelsWithNAN = labelsCol(indSortidxWithNAN);

rmPos = find(isnan(IndIntsWithNAN));
keepAfter = max(rmPos)+1;
IndIntsAll = IndIntsWithNAN(keepAfter:end);
DepIntsAll = DepIntsWithNAN(keepAfter:end);
LabelsAll = LabelsWithNAN(keepAfter:end);

end

%% getnCells
function nCells = getnCells(intiMs)
nExpts = length(intiMs);
nCells = 0;
for iExpt = 1:nExpts
    CellIDs = [];
    cellLabel = intiMs{iExpt}.label;
    nKTs = size(cellLabel,1);
    for iKT = 1:nKTs
         CellIDs = cat(2, CellIDs, str2num(cellLabel(iKT,2:4)));
    end
    wkCells = unique(CellIDs);
    nCells = nCells + length(wkCells);
end
end
%% rmLowCenpC
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

%% getCellPhaseBinSize
function binSize = getCellPhaseBinSize(intMatrix)
binSize = [];
for iCol = 1:size(intMatrix,2)
    wkData = intMatrix(:, iCol);
    nonNanIdx = find(~isnan(wkData));
    if isempty(nonNanIdx)
        wkSize = 0;
    else
        wkSize = length(nonNanIdx);
    end
    binSize = cat(2, binSize, wkSize);
end

end

%% getnCellsPhase
function nCells = getnCellsPhase(labelData)
CellVector = [];
for iCol = 1:size(labelData,2)
    for iRow = 1:size(labelData,1)
        assessData = labelData{iRow, iCol};
        if ~iscell(assessData)
            wkCell = str2num(assessData(1:4)); %this is the first 4 numbers of the cell label, i.e. the experiment number and cell number
            CellVector = cat(1, CellVector, wkCell);
        end
    end
end

nCells = length(unique(CellVector));

end
%% getSlopeIntersectPoint
function desiredYVal = getSlopeIntersectPoint(xVal1, yVal1, xVal2, yVal2, desiredXVal)

slope = (yVal1-yVal2)/(xVal1-xVal2);
intercept = yVal1-(slope*xVal1); %y=mx+c to y-mx=c
desiredYVal = (slope*desiredXVal)+intercept;


end