function intData = makeIntPlots(ctrlIntiMs, treatIntiMs, varargin)
% MAKEINTPLOTS Produces two plots of multiple kinetchore intensity
% measurements over multiple cells and experiments split into bins (one control, one treatment).
% 
%   MAKEINTPLOTS(ctrlIntiMs, treatIntiMs).
%   ctrlIntiMs and treatIntiMs must be organised as {Mad2intiM1,
%   otherintiM1; Mad2intiM2, otherintiM2; ...}.
%   treatIntiMs are optional so can be replaced with []. treatIntiMs must
%   be in the same order as ctrlIntiMs (so that experiments are matched).
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
%   indMarker: {'Venus-MAD2'} or string of intensity marker that you binned
%       on. Will request this later if not filled in now.
%
%   depMarker: {nan} or string of second intensity marker. Will request
%       this later if not filled in now.
%
%   normMarker: {'CENP-C'} or string of reference kinetochore marker that
%       you (may) normalise to. Will not request this later.
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
opts.depMarker = NaN;
opts.normMarker = 'CENP-C';
opts.depInner = 0;
opts.minBinFirst = 0;
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = size(ctrlIntiMs, 1);
isTreat = ~isempty(treatIntiMs);

indCtrl_old = [];
indCtrl = [];
indTreat = [];
depCtrl = [];
depTreat = [];

for iExpt = 1:nExpts
    [innerCI, outerCI] = rmLowCenpC(ctrlIntiMs{iExpt,1}); %20250429
    [innerCD, outerCD] = rmLowCenpC(ctrlIntiMs{iExpt,2});
    ctrlIntiMs{iExpt,1}.intensity.mean.inner = innerCI; %20250429
    ctrlIntiMs{iExpt,1}.intensity.mean.outer = outerCI;
    ctrlIntiMs{iExpt,2}.intensity.mean.inner = innerCD; %20250429
    ctrlIntiMs{iExpt,2}.intensity.mean.outer = outerCD;
    [iCBinned, dCBinned, iCAll, dCAll] = binIntiMs(ctrlIntiMs{iExpt,1}, ctrlIntiMs{iExpt,2}, 'nBins', nBins, 'IndNorm', opts.normIndInt, 'DepNorm', opts.normDepInt, 'InnerDep', opts.depInner, 'minIntFirstCol', opts.minBinFirst);
    iCMed = median(iCBinned);
    dCMed = median(dCBinned);
    if isTreat
        [innerTI, outerTI] = rmLowCenpC(treatIntiMs{iExpt,1}); %20250429
        [innerTD, outerTD] = rmLowCenpC(treatIntiMs{iExpt,2});
        treatIntiMs{iExpt,1}.intensity.mean.inner = innerTI; %20250429
        treatIntiMs{iExpt,1}.intensity.mean.outer = outerTI;
        treatIntiMs{iExpt,2}.intensity.mean.inner = innerTD; %20250429
        treatIntiMs{iExpt,2}.intensity.mean.outer = outerTD;
        [iTBinned, dTBinned, iTAll, dTAll] = binIntiMs(treatIntiMs{iExpt,1}, treatIntiMs{iExpt,2}, 'nBins', nBins, 'IndNorm', opts.normIndInt, 'DepNorm', opts.normDepInt, 'InnerDep', opts.depInner, 'minIntFirstCol', opts.minBinFirst);
        iTMed = median(iTBinned);
        dTMed = median(dTBinned);
    end
    
    if isNorm
        iMax = max(iCMed);
        dMax = max(dCMed);
        iCBinned = iCBinned/iMax;
        dCBinned = dCBinned/dMax;
        iCAll = iCAll/iMax;
        dCAll = dCAll/dMax;
        if isTreat
            iTBinned = iTBinned/iMax;
            dTBinned = dTBinned/dMax;
            iTAll = iTAll/iMax;
            dTAll = dTAll/dMax;
        end
        
    end
    %the following is so the whole data set is considered rather than just
    %the data points that had already been binned, potentially increasing
    %the number of data points in the final part
    indCtrl_old = cat(1, indCtrl_old, iCBinned);
    indCtrl = cat(1, indCtrl, iCAll); %20250429
    %depCtrl = cat(1, depCtrl, dCBinned);
    depCtrl = cat(1, depCtrl, dCAll); %20250429
    if isTreat
        %indTreat = cat(1, indTreat, iTBinned);
        indTreat = cat(1, indTreat, iTAll); %20250429
        %depTreat = cat(1, depTreat, dTBinned);
        depTreat = cat(1, depTreat, dTAll); %20250429
    end
end
   
[iCIntsGrouped, dCIntsGrouped] = rebinIntensities(indCtrl, depCtrl, 'nBins', nBins, 'minIntFirstCol', opts.minBinFirst);
cBinSz = size(iCIntsGrouped,1);
[indCIntsMed, indCIntsNeg, indCIntsPos] = makeQuartiles(iCIntsGrouped);
[depCIntsMed, depCIntsNeg, depCIntsPos] = makeQuartiles(dCIntsGrouped);
[~, indCIntsLCI, indCIntsUCI] = makeMedianConfInts(iCIntsGrouped);
[~, depCIntsLCI, depCIntsUCI] = makeMedianConfInts(dCIntsGrouped);
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
intData.control.independent.intensities = iCIntsGrouped;
CIstats.median = indCIntsMed;
CIstats.lowerCI = indCIntsLCI;
CIstats.upperCI = indCIntsUCI;
CIstats.lowerQuart = indCIntsNeg;
CIstats.upperQuart = indCIntsPos;
intData.control.independent.stats = CIstats;

intData.control.dependent.intensities = dCIntsGrouped;
CDstats.median = depCIntsMed;
CDstats.lowerCI = depCIntsLCI;
CDstats.upperCI = depCIntsUCI;
CDstats.lowerQuart = depCIntsNeg;
CDstats.upperQuart = depCIntsPos;
intData.control.dependent.stats = CDstats;

nCells_ctrl = getnCells(ctrlIntiMs);
intData.control.nCells = nCells_ctrl;
intData.control.binSize = cBinSz;

if isTreat
    [iTIntsGrouped, dTIntsGrouped] = rebinIntensities(indTreat, depTreat, 'nBins', nBins, 'minIntFirstCol', opts.minBinFirst);
    tBinSz = size(iTIntsGrouped,1);
    [indTIntsMed, indTIntsNeg, indTIntsPos] = makeQuartiles(iTIntsGrouped);
    [depTIntsMed, depTIntsNeg, depTIntsPos] = makeQuartiles(dTIntsGrouped);
    [~, indTIntsLCI, indTIntsUCI] = makeMedianConfInts(iTIntsGrouped);
    [~, depTIntsLCI, depTIntsUCI] = makeMedianConfInts(dTIntsGrouped);
    if isNorm
        iTIntsGrouped = iTIntsGrouped/indCMaxMed;
        indTIntsMed = indTIntsMed/indCMaxMed;
        indTIntsNeg = indTIntsNeg/indCMaxMed;
        indTIntsPos = indTIntsPos/indCMaxMed;
        indTIntsLCI = indTIntsLCI/indCMaxMed;
        indTIntsUCI = indTIntsUCI/indCMaxMed;

        dTIntsGrouped = dTIntsGrouped/depCMaxMed;
        depTIntsMed = depTIntsMed/depCMaxMed;
        depTIntsNeg = depTIntsNeg/depCMaxMed;
        depTIntsPos = depTIntsPos/depCMaxMed;
        depTIntsLCI = depTIntsLCI/depCMaxMed;
        depTIntsUCI = depTIntsUCI/depCMaxMed;
    end
    intData.treatment.independent.intensities = iTIntsGrouped;
    TIstats.median = indTIntsMed;
    TIstats.lowerCI = indTIntsLCI;
    TIstats.upperCI = indTIntsUCI;
    TIstats.lowerQuart = indTIntsNeg;
    TIstats.upperQuart = indTIntsPos;
    intData.treatment.independent.stats = TIstats;
    
    intData.treatment.dependent.intensities = dTIntsGrouped;
    TDstats.median = depTIntsMed;
    TDstats.lowerCI = depTIntsLCI;
    TDstats.upperCI = depTIntsUCI;
    TDstats.lowerQuart = depTIntsNeg;
    TDstats.upperQuart = depTIntsPos;
    intData.treatment.dependent.stats = TDstats;
    
    nCells_treat = getnCells(treatIntiMs);
    intData.treatment.nCells = nCells_treat;
    intData.treatment.binSize = tBinSz;

end

halfChanges_ctrl = getHalfChanges(intData, 'dataType', 'int', 'isCtrl', 1);
intData.control.independent.stats.HCint = halfChanges_ctrl.Mad2.HCintensity;
intData.control.independent.stats.HCbin = halfChanges_ctrl.Mad2.HCbin;
intData.control.dependent.stats.HCint = halfChanges_ctrl.Other.HCintensityKKDelta;
intData.control.dependent.stats.HCbin = halfChanges_ctrl.Other.HCbin;
intData.control.dependent.stats.HCMad2Int = halfChanges_ctrl.Other.HCMad2Int;

if isTreat
    halfChanges_treat = getHalfChanges(intData, 'dataType', 'int', 'isCtrl', 0);
    intData.treatment.independent.stats.HCint = halfChanges_treat.Mad2.HCintensity;
    intData.treatment.independent.stats.HCbin = halfChanges_treat.Mad2.HCbin;
    intData.treatment.dependent.stats.HCint = halfChanges_treat.Other.HCintensityKKDelta;
    intData.treatment.dependent.stats.HCbin = halfChanges_treat.Other.HCbin;
    intData.treatment.dependent.stats.HCMad2Int = halfChanges_treat.Other.HCMad2Int;
end



%%
if isTreat
    treatXVals = 1:nBins;
end

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
if any([isnan(opts.indMarker), isnan(opts.depMarker), isnan(opts.normMarker)])
    labelAnswers = inputdlg({'Enter independent intensity marker:', 'Enter dependent intensity marker:', ...
        'Enter marker used for normalisation (if appropriate):'}, 'Marker proteins', [1 35], {'Venus-MAD2', 'SKAP', 'CENP-C'});
else
    labelAnswers = {opts.indMarker, opts.depMarker, opts.normMarker};
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
    'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5); %median

line('XData',[intData.control.independent.stats.HCbin, intData.control.independent.stats.HCbin], ...
    'YData',[0 1], 'Color', indCols(1,:), 'LineStyle', '--'); %half change

yyaxis right
patch([ctrlXVals fliplr(ctrlXVals)], ...
    [depCIntsNeg fliplr(depCIntsPos)], depCols(1,:), ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles

patch([ctrlXVals fliplr(ctrlXVals)], ...
    [depCIntsLCI fliplr(depCIntsUCI)], depCols(1,:), ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none') %95pct confidence interval

plot(ctrlXVals, depCIntsMed, '-', 'Color', depCols(1,:), ...
    'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5); %median

line('XData',[intData.control.dependent.stats.HCbin, intData.control.dependent.stats.HCbin], ...
    'YData',[0 1], 'Color', depCols(1,:), 'LineStyle', '--'); %half change


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
    indDivSignaller = strcat('/', labelAnswers{3});
else
    indDivSignaller = '';
end
if isNorm
    labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'YLim', [lowerCtrlLimit upperCtrlLimit], 'TickDir', 'out');
else
    labelInd = sprintf('%s%s intensity', labelAnswers{1}, indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out');
end
ylabel(labelInd);

yyaxis right
if opts.normDepInt
    depDivSignaller = strcat('/', labelAnswers{3});
else
    depDivSignaller = '';
end

ctrlXTicks = 1:nBins;
lowerX = min(ctrlXVals)-0.5;
upperX = max(ctrlXVals)+0.5;


if isNorm
    labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
    set(gca, 'YColor', depCols(1,:), 'YLim', [lowerCtrlLimit upperCtrlLimit], ...
        'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
    set(gca, 'YColor', depCols(1,:), 'XTick', ctrlXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
if opts.minBinFirst
    set(gca, 'XDir', 'reverse');
end
ylabel(labelDep);
xlabel(sprintf('%s%s pseudo-timeline bins', labelAnswers{1}, indDivSignaller));
%% treatment figure
if isTreat
    treatFig = figure;
    hold on
    yyaxis left
    % use this for quartiles + median
    patch([treatXVals fliplr(treatXVals)], ...
        [indTIntsNeg fliplr(indTIntsPos)], indCols(2,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles

    patch([treatXVals fliplr(treatXVals)], ...
        [indTIntsLCI fliplr(indTIntsUCI)], indCols(2,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence intervals

    plot(treatXVals, indTIntsMed, '-', 'Color', indCols(2,:), ...
        'MarkerFaceColor', indCols(2,:), 'LineWidth', 1.5); %median

    line('XData',[intData.treatment.independent.stats.HCbin, intData.treatment.independent.stats.HCbin], ...
        'YData',[0 1], 'Color', indCols(2,:), 'LineStyle', '--'); %half change

    yyaxis right
    patch([treatXVals fliplr(treatXVals)], ...
        [depTIntsNeg fliplr(depTIntsPos)], depCols(2,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles

    patch([treatXVals fliplr(treatXVals)], ...
        [depTIntsLCI fliplr(depTIntsUCI)], depCols(2,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence intervals

    plot(treatXVals, depTIntsMed, '-', 'Color', depCols(2,:), ...
        'MarkerFaceColor', depCols(2,:), 'LineWidth', 1.5); %median

    line('XData',[intData.treatment.dependent.stats.HCbin, intData.treatment.dependent.stats.HCbin], ...
        'YData',[0 1], 'Color', depCols(2,:), 'LineStyle', '--'); %half change
    
    

    hold off
    treatax = gca;
    yyaxis left
    IndTreatTicks = yticks;
    yyaxis right
    DepTreatTicks = yticks;
    upperTreatLimit = max([IndTreatTicks(:); DepTreatTicks(:)]);
    lowerTreatLimit = min([IndTreatTicks(:); DepTreatTicks(:)]);
    
    yyaxis left
    if opts.normIndInt
        indDivSignaller = strcat('/', labelAnswers{3});
    else
        indDivSignaller = '';
    end
    if isNorm
        labelInd = sprintf('%s%s intensity, normalised (AU)', labelAnswers{1}, indDivSignaller);
        set(gca, 'YColor', indCols(1,:), 'YLim', [lowerTreatLimit upperTreatLimit], 'TickDir', 'out');
    else
        labelInd = sprintf('%s%s intensity (AU)', labelAnswers{1}, indDivSignaller);
        set(gca, 'YColor', indCols(1,:), 'TickDir', 'out');
    end
    ylabel(labelInd);
    
    yyaxis right
    if opts.normDepInt
        depDivSignaller = strcat('/', labelAnswers{3});
    else
        depDivSignaller = '';
    end
    
    
    treatXTicks = 1:nBins;
    lowerX = min(treatXVals)-0.5;
    upperX = max(treatXVals)+0.5;
    
    
    if isNorm
        labelDep = sprintf('%s%s intensity, normalised (AU)', labelAnswers{2}, depDivSignaller);
        set(gca, 'YColor', depCols(1,:), 'YLim', [lowerTreatLimit upperTreatLimit], ...
            'XTick', treatXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
    else
        labelDep = sprintf('%s%s intensity (AU)', labelAnswers{2}, depDivSignaller);
        set(gca, 'YColor', depCols(1,:), 'XTick', treatXTicks, 'XLim', [lowerX upperX], 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
    end
    if opts.minBinFirst
        set(gca, 'XDir', 'normal');
    end
    ylabel(labelDep);

    xlabel(sprintf('%s%s pseudo-timeline bins', labelAnswers{1}, indDivSignaller));

end
   
%% Resizing windows and axes
if isNorm
    if isTreat
        minValue = min([indCIntsNeg, depCIntsNeg, indTIntsNeg, depTIntsNeg], [], 'all');
        maxMinAx = max([lowerCtrlLimit, lowerTreatLimit], [], 'all');
        minMinAx = min([lowerCtrlLimit, lowerTreatLimit], [], 'all');
        possVec = [maxMinAx, minMinAx, -0.2, -0.1, -0.05];
    else
        minValue = min([indCIntsNeg, depCIntsNeg], [], 'all');
        possVec = [lowerCtrlLimit, -0.5, -0.2, -0.1, -0.05];
    end
    
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
    yyaxis left
    set(ctrlax, 'YLim', [lowerCtrlLimit upperCtrlLimit]);
    yyaxis right
    set(ctrlax, 'YLim', [lowerCtrlLimit upperCtrlLimit]);
    ctrlRange = upperCtrlLimit - lowerCtrlLimit;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    %%set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, 275, ctrlRange*200])
    set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlRange*200])
    set(ctrlFig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (ctrlRange*200)+(ctrlPos(2)*2)])
    set(ctrlax, 'XTickLabelRotation', 0)
    
    if isTreat
        figure(treatFig);
        yyaxis left
        set(treatax, 'YLim', [lowerCtrlLimit upperTreatLimit]);
        yyaxis right
        set(treatax, 'YLim', [lowerCtrlLimit upperTreatLimit]);
        treatRange = upperTreatLimit - lowerCtrlLimit;
        set(treatax, 'Units', 'pixels')
        treatPos = get(treatax, 'Position');
        set(treatax, 'Position', [treatPos(1), treatPos(2)*1.25, 275, treatRange*200])
        set(treatFig, 'Position', [700 60 treatPos(3)+(treatPos(1)*2) (treatRange*200)+(treatPos(2)*2)])
        set(treatax, 'XTickLabelRotation', 0)
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
function [rebinnedIndInts, rebinnedDepInts] = rebinIntensities(IndIntsOrg, DepIntsOrg, varargin)

opts.nBins = 25; %number of bins to make
opts.minIntFirstCol = 0; %whether minimum intensities are in first column (1) or last column (0)
opts = processOptions(opts,varargin{:});

%making column vectors to work with mink/maxk
indIntsCol = IndIntsOrg(:);
depIntsCol = DepIntsOrg(:);

%setting up loops
Total = length(indIntsCol);
n = floor(Total/opts.nBins); %bin starting position
k = n; %number of KTs per bin
bin = 0;

IndIntsBinned = nan(k, opts.nBins);
DepIntsBinned = nan(k, opts.nBins);

while n < (Total + 1)
    bin = bin + 1;
    
    [maxIndSubset, maxIdx] = maxk(indIntsCol, n); %gets highest n values and their positions (ordered highest to lowest)
    [minIndSubset, minIdx] = mink(maxIndSubset, k); %gets lowest k values from highest n values and their positions (ordered lowest to highest)

    maxDepSubset = depIntsCol(maxIdx); %gets corresponding depInt values from highest n indInt values
    minDepSubset = maxDepSubset(minIdx); %gets corresponding depInt values from lowest k of highest n indInt values

    IndIntsBinned(:, bin) = minIndSubset;
    DepIntsBinned(:, bin) = minDepSubset;
    
    n = n + k;
end

if opts.minIntFirstCol
    rebinnedIndInts = fliplr(IndIntsBinned); %if first column is min, then top left should be min value
    rebinnedDepInts = fliplr(DepIntsBinned); %flips correspondingly
else
    rebinnedIndInts = flipud(IndIntsBinned); %if last column is min, flip matrix accordingly
    rebinnedDepInts = flipud(DepIntsBinned); %flips correspondingly
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
function [IndIntsBinned, DepIntsBinned, IndIntsAll, DepIntsAll] = binIntiMs(indIntiM, depIntiM, varargin)
% Function to enable binning of two intiMs. The order of intiM binning
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

%setting up loops
Total = size(find(~isnan(indIntiM.intensity.mean.outer)),1);
n = floor(Total/opts.nBins); %bin starting position
k = n; %number of KTs per bin
bin = 0;
IndIntsBinned = nan(k, opts.nBins);
DepIntsBinned = nan(k, opts.nBins);

while n < (Total + 1)
    bin = bin + 1;
    
    [maxIndSubset, maxIdx] = maxk(indIntsCol, n); %gets highest n values and their positions (ordered highest to lowest)
    [minIndSubset, minIdx] = mink(maxIndSubset, k); %gets lowest k values from highest n values and their positions (ordered lowest to highest)
    
    
    %minDepSubset = depIntsCol(minIdx); %gets corresponding depInt values from lowest n indInt values
    %maxDepSubset = minDepSubset(maxIdx); %gets corresponding depInt values from highest k of lowest n indInt values

    maxDepSubset = depIntsCol(maxIdx); %gets corresponding depInt values from highest n indInt values
    minDepSubset = maxDepSubset(minIdx); %gets corresponding depInt values from lowest k of highest n indInt values

    IndIntsBinned(:, bin) = minIndSubset;
    DepIntsBinned(:, bin) = minDepSubset;
    
    n = n + k;
end

if opts.minIntFirstCol
    IndIntsBinned = fliplr(IndIntsBinned); %if first column is min, then top left should be min value
    DepIntsBinned = fliplr(DepIntsBinned); %flips correspondingly
else
    IndIntsBinned = flipud(IndIntsBinned); %if last column is min, flip matrix accordingly
    DepIntsBinned = flipud(DepIntsBinned); %flips correspondingly
end

indIntsRM = rmmissing(indIntsCol);
[IndIntsAll, indSortidx] = sort(indIntsRM, 'descend');
depIntsRM = rmmissing(depIntsCol);
DepIntsAll = depIntsRM(indSortidx);

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
nExpts = size(intiMs,1);
nCells = 0;
for iExpt = 1:nExpts
    CellIDs = [];
    cellLabel = intiMs{iExpt,1}.label;
    nKTs = size(cellLabel,1);
    for iKT = 1:nKTs
         CellIDs = cat(2, CellIDs, str2num(cellLabel(iKT,2:4)));
    end
    wkCells = unique(CellIDs);
    nCells = nCells + length(wkCells);
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