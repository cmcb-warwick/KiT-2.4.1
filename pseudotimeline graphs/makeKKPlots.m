function KKIntData = makeKKPlots(ctrl_iMintiMPairs, treat_iMintiMPairs, varargin)
% MAKEKKPLOTS Produces a single plot of averaged KK distance
% measurements over multiple cells and experiments split into bins.
% Averaging is achieved by detecting whether the mean Mad2 bin of a pair of
% KTs with sequential bin values (i.e. bin 1 and 2) is closer to the lower
% or higher bin, based on position of intensity value within the bin.
% 
%    MAKEKKPLOTS(iMintiMPairs) intiMPairs must be organised as 
%   {iM intiM; iM2 intiM2} and so on. Maximum of five sets of iM/intiM
%   structures for now. Options are available:
%
%    Options, defaults in {}:-
%
%    nBins: {25} or other integer. The number of bins to split each intiM
%       structure into.
%
%   normIndInt: 0 or {1}. Whether to normalise the independent intensity
%       marker (i.e. Mad2) to CenpC or not.
%
%    normalise: 0 or {1}. Whether to normalise the median of the maximum bin for
%       each intensity class to one. 
%
%   swapBin: {0} or 1. Whether to have maximum intensity as minimum (0) or
%       maximum (1) bin number when plotting.
%
% Copyright (c) 2023 C. C. Conway

opts.nBins = 25;
opts.normIndInt = 1;
opts.normalise = 1;
opts.swapBin = 0;
opts = processOptions(opts, varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = size(ctrl_iMintiMPairs, 1);
isTreat = ~isempty(treat_iMintiMPairs);
nCCells = 0;
nTCells = 0;
ctrlIntiMScaled = [];
allcMax = [];
if isTreat
    treatIntiMScaled = [];
    alltMax = [];
end
for iExpt = 1:nExpts
    %step 1: remove CenpC intensities that are too low
    ctrl_iMintiMPairs{iExpt,2}.intensity.mean.inner = rmLowCenpC(ctrl_iMintiMPairs{iExpt,2});
    cIntiM = ctrl_iMintiMPairs{iExpt,2}.intensity.mean.outer;
    if opts.normIndInt
        cIntiM = cIntiM ./ ctrl_iMintiMPairs{iExpt,2}.intensity.mean.inner;
    end
    %then bin (max values retained)
    ctrlIntiM = makeBins(nBins, cIntiM);
    cMed = median(ctrlIntiM);


    if isTreat
        %remove low CenpC intensities
        treat_iMintiMPairs{iExpt,2}.intensity.mean.inner = rmLowCenpC(treat_iMintiMPairs{iExpt,2});
        tIntiM = treat_iMintiMPairs{iExpt,2}.intensity.mean.outer;
        if opts.normIndInt
            tIntiM = tIntiM ./ treat_iMintiMPairs{iExpt,2}.intensity.mean.inner;
        end
        %bin
        treatIntiM = makeBins(nBins, tIntiM);
    end
    
    if isNorm
        %get max median of ctrl
        cMax = max(cMed);
        %put in a matrix so we know what to divide by later on
        allcMax = cat(2, allcMax, cMax);
        %normalise intiM
        wkCtrlIntiM = ctrlIntiM/cMax;
        ctrlIntiMScaled = cat(1, ctrlIntiMScaled, wkCtrlIntiM);
        if isTreat
            wkTreatIntiM = treatIntiM/cMax;
            %normalise treatment intiM with control intensity data
            treatIntiMScaled = cat(1, treatIntiMScaled, wkTreatIntiM);
        end
    end
end

ctrlIntiMScaledBinned = makeBins(nBins, ctrlIntiMScaled);
ctrlKKnormCell = [];
if isTreat
    treatIntiMScaledBinned = makeBins(nBins, treatIntiMScaled);
    treatKKnormCell = [];
end
for iExpt = 1:nExpts
    [ctrlKKnormCell, ctrlnCells] = getKKBins_withNormInt(ctrl_iMintiMPairs{iExpt,1}, ctrl_iMintiMPairs{iExpt,2}, ctrlIntiMScaledBinned, allcMax(iExpt), nBins, ctrlKKnormCell, opts.normIndInt);
    nCCells = nCCells + ctrlnCells;
    if isTreat
        [treatKKnormCell, treatnCells] = getKKBins_withNormInt(treat_iMintiMPairs{iExpt,1}, treat_iMintiMPairs{iExpt,2}, treatIntiMScaledBinned, allcMax(iExpt), nBins, treatKKnormCell, opts.normIndInt);
        nTCells = nTCells + treatnCells;
    end
end

[MedianCtrlNormKK, LQCtrlNormKK, UQCtrlNormKK, avInBinCtrl, pairsBinsCtrl] = getKKquartiles(ctrlKKnormCell);
[~, LCICtrlNormKK, UCICtrlNormKK, ~, ~] = getKKConfInts(ctrlKKnormCell);
[iCtrlNormMed, iCtrlNormLQ, iCtrlNormUQ] = makeQuartiles(ctrlIntiMScaledBinned);
[~, iCtrlNormLCI, iCtrlNormUCI] = makeMedianConfInts(ctrlIntiMScaledBinned);
if isNorm
    maxInt = max(iCtrlNormMed);
    iCtrlNormMed = iCtrlNormMed/maxInt;
    iCtrlNormLQ = iCtrlNormLQ/maxInt;
    iCtrlNormUQ = iCtrlNormUQ/maxInt;
    iCtrlNormLCI = iCtrlNormLCI/maxInt;
    iCtrlNormUCI = iCtrlNormUCI/maxInt;
    ctrlIntiMScaledBinned = ctrlIntiMScaledBinned/maxInt;
end
KKIntData.control.KK = ctrlKKnormCell;
KKIntData.control.nKKCells = nCCells;
KKIntData.control.nIntCells = getnCells(ctrl_iMintiMPairs(:,2));
KKIntData.control.intensities = ctrlIntiMScaledBinned;
CIstats.median = iCtrlNormMed;
CIstats.lowerCI = iCtrlNormLCI;
CIstats.upperCI = iCtrlNormUCI;
CIstats.lowerQuart = iCtrlNormLQ;
CIstats.upperQuart = iCtrlNormUQ;
KKIntData.control.intStats = CIstats;
KKIntData.control.intBinSize = size(ctrlIntiMScaledBinned,1);
KKIntData.control.KKBinSize = pairsBinsCtrl;
CKKstats.median = MedianCtrlNormKK;
CKKstats.lowerCI = LCICtrlNormKK;
CKKstats.upperCI = UCICtrlNormKK;
CKKstats.lowerQuart = LQCtrlNormKK;
CKKstats.upperQuart = UQCtrlNormKK;
KKIntData.control.KKStats = CKKstats;

%[~, lMADNormKK, uMADNormKK] = getKKMedianMAD(KKnormCell); %in case I
%wanted to do MAD rather than quartiles - will stick with quartiles

if isTreat
    [MedianTreatNormKK, LQTreatNormKK, UQTreatNormKK, avInBinTreat, pairsBinsTreat] = getKKquartiles(treatKKnormCell);
    [~, LCITreatNormKK, UCITreatNormKK, ~, ~] = getKKConfInts(treatKKnormCell);
    [iTreatNormMed, iTreatNormLQ, iTreatNormUQ] = makeQuartiles(treatIntiMScaledBinned);
    [~, iTreatNormLCI, iTreatNormUCI] = makeMedianConfInts(treatIntiMScaledBinned);
    if isNorm
        iTreatNormMed = iTreatNormMed/maxInt;
        iTreatNormLQ = iTreatNormLQ/maxInt;
        iTreatNormUQ = iTreatNormUQ/maxInt;
        iTreatNormLCI = iTreatNormLCI/maxInt;
        iTreatNormUCI = iTreatNormUCI/maxInt;
        treatIntiMScaledBinned = treatIntiMScaledBinned/maxInt;
    end
    KKIntData.treatment.KK = treatKKnormCell;
    KKIntData.treatment.nKKCells = nTCells;
    KKIntData.treatment.nIntCells = getnCells(treat_iMintiMPairs(:,2));
    KKIntData.treatment.intensities = treatIntiMScaledBinned;
    TIstats.median = iTreatNormMed;
    TIstats.lowerCI = iTreatNormLCI;
    TIstats.upperCI = iTreatNormUCI;
    TIstats.lowerQuart = iTreatNormLQ;
    TIstats.upperQuart = iTreatNormUQ;
    KKIntData.treatment.intStats = TIstats;
    KKIntData.treatment.intBinSize = size(treatIntiMScaledBinned,1);
    KKIntData.treatment.KKBinSize = pairsBinsTreat;
    TKKstats.median = MedianTreatNormKK;
    TKKstats.lowerCI = LCITreatNormKK;
    TKKstats.upperCI = UCITreatNormKK;
    TKKstats.lowerQuart = LQTreatNormKK;
    TKKstats.upperQuart = UQTreatNormKK;
    KKIntData.treatment.KKStats = TKKstats;

    halfChanges = getHalfChanges(KKIntData, 'dataType', 'kk', 'isCtrl', 0);
    KKIntData.treatment.intStats.HCint = halfChanges.Mad2.HCintensity;
    KKIntData.treatment.intStats.HCbin = halfChanges.Mad2.HCbin;
    KKIntData.treatment.KKStats.HCKK = halfChanges.Other.HCintensityKKDelta;
    KKIntData.treatment.KKStats.HCbin = halfChanges.Other.HCbin;
    KKIntData.treatment.KKStats.HCMad2Int = halfChanges.Other.HCMad2Int;

end
halfChanges = getHalfChanges(KKIntData, 'dataType', 'kk', 'isCtrl', 1);
KKIntData.control.intStats.HCint = halfChanges.Mad2.HCintensity;
KKIntData.control.intStats.HCbin = halfChanges.Mad2.HCbin;
KKIntData.control.KKStats.HCKK = halfChanges.Other.HCintensityKKDelta;
KKIntData.control.KKStats.HCbin = halfChanges.Other.HCbin;
KKIntData.control.KKStats.HCMad2Int = halfChanges.Other.HCMad2Int;

xvals = 1:nBins;

minX = min(xvals);
maxX = max(xvals);

indCols = [0 0 0;...
    0.4 0.4 0.4];
%black and grey
%#691787
depCols = [0.8118 0.1333 0.5333;...
            0.4471 0.0667 0.4275;...
            0.9451 0.7176 0.8510;...
            0.3137 0.6667 0.1882;...
            0.1098 0.6745 0.3412;...
            0.7255 0.8784 0.5255];
%% control figure
ctrlFig = figure;
hold on
    yyaxis left

    if opts.swapBin
        wkxVals = xvals - (nBins+1); %shift to negative 
        finxVals = abs(wkxVals); %make positive
        xvals = finxVals; %this is a better option than fliplr because that would mess up sub-bin localisation - this is less like mirroring the x-axis and more like just replacing the axis tick numbers
    end

    patch([xvals fliplr(xvals)], ...
        [iCtrlNormLQ fliplr(iCtrlNormUQ)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles

    patch([xvals fliplr(xvals)], ...
        [iCtrlNormLCI fliplr(iCtrlNormUCI)], indCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence intervals

    plot(xvals, iCtrlNormMed, '-', 'Color', indCols(1,:), ...
        'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5); %median

    line('XData',[KKIntData.control.intStats.HCbin, KKIntData.control.intStats.HCbin], ... 
        'YData',[0 1], 'Color', indCols(1,:), 'LineStyle', '--'); %half changes

    yyaxis right
    patch([xvals fliplr(xvals)], ...
        [LQCtrlNormKK fliplr(UQCtrlNormKK)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles

    patch([xvals fliplr(xvals)], ...
        [LCICtrlNormKK fliplr(UCICtrlNormKK)], depCols(1,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence intervals

    plot(xvals, MedianCtrlNormKK, '-', 'Color', depCols(1,:), ...
        'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5); %median

    line('XData',[KKIntData.control.KKStats.HCbin, KKIntData.control.KKStats.HCbin], ...
        'YData',[0.9 1.2], 'Color', depCols(1,:), 'LineStyle', '--'); %half changes

    hold off

hold on
if opts.normIndInt
    indDivSignaller = '/CENP-C';
else
    indDivSignaller = '';
end
yyaxis left
if isNorm
    labelInd = sprintf('Venus-MAD2%s intensity, normalised (AU)', indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
else
    labelInd = sprintf('Venus-MAD2%s intensity (AU)', indDivSignaller);
    set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
end
ylabel(labelInd);

yyaxis right

labelDep = 'K-K distance (µm)';
set(gca, 'YColor', depCols(1,:), 'TickDir', 'out', 'YGrid', 'on')

set(gca, 'XTick', 1:nBins, 'XLim', [(minX-1) (maxX+1)], 'TickDir', 'out');


xlabel(sprintf('Venus-MAD2%s pseudo-timeline bins', indDivSignaller))

if opts.swapBin
    set(gca, 'XDir', 'reverse');
end
ylabel(labelDep);
ctrlax = gca;
yyaxis left
minIntCtrlTick = min(yticks);
maxIntCtrlTick = max(yticks);
yyaxis right 
minKKCtrlTick = min(yticks);
maxKKCtrlTick = max(yticks);
hold off
%% treatment figure
if isTreat
    treatFig = figure;
    hold on
        yyaxis left
        patch([xvals fliplr(xvals)], ...
            [iTreatNormLQ fliplr(iTreatNormUQ)], indCols(2,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xvals fliplr(xvals)], ...
            [iTreatNormLCI fliplr(iTreatNormUCI)], indCols(2,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence interval
        plot(xvals, iTreatNormMed, '-', 'Color', indCols(2,:), ...
            'MarkerFaceColor', indCols(1,:), 'LineWidth', 1.5);
        line('XData',[KKIntData.treatment.intStats.HCbin, KKIntData.treatment.intStats.HCbin], 'YData',[0 1], 'Color', indCols(2,:), 'LineStyle', '--');
        yyaxis right
        patch([xvals fliplr(xvals)], ...
            [LQTreatNormKK fliplr(UQTreatNormKK)], depCols(2,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %quartiles
        patch([xvals fliplr(xvals)], ...
            [LCITreatNormKK fliplr(UCITreatNormKK)], depCols(2,:), ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none') %95pct confidence interval
        plot(xvals, MedianTreatNormKK, '-', 'Color', depCols(2,:), ...
            'MarkerFaceColor', depCols(1,:), 'LineWidth', 1.5);
        line('XData',[KKIntData.treatment.KKStats.HCbin, KKIntData.treatment.KKStats.HCbin], 'YData',[0.7 1], 'Color', depCols(2,:), 'LineStyle', '--');
        hold off
    
    hold on
    
    yyaxis left
    if isNorm
        labelInd = sprintf('Venus-MAD2%s intensity, normalised (AU)', indDivSignaller);
        set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
    else
        labelInd = sprintf('Venus-MAD2%s intensity (AU)', indDivSignaller);
        set(gca, 'YColor', indCols(1,:), 'TickDir', 'out', 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 9);
    end
    ylabel(labelInd);
    
    yyaxis right
    
    labelDep = 'K-K distance (µm)';
    set(gca, 'YColor', depCols(1,:), 'TickDir', 'out', 'YGrid', 'on')
    
    set(gca, 'XTick', 1:nBins, 'XLim', [(minX-1) (maxX+1)], 'TickDir', 'out');

    xlabel(sprintf('Venus-MAD2%s pseudo-timeline bins', indDivSignaller))


    if opts.swapBin
        set(gca, 'XDir', 'reverse');
    end

    ylabel(labelDep);
    treatax = gca;
    yyaxis left
    minIntTreatTick = min(yticks);
    maxIntTreatTick = max(yticks);
    yyaxis right 
    minKKTreatTick = min(yticks);
    maxKKTreatTick = max(yticks);
    hold off
end
%%
if isNorm
    if isTreat
        %int first
        minIntValue = min([iCtrlNormLQ, iTreatNormLQ], [], 'all');
        minIntMinAx = min([minIntCtrlTick, minIntTreatTick], [], 'all');
        maxIntMinAx = max([minIntCtrlTick, minIntTreatTick], [], 'all');
        possIntVec = [maxIntMinAx, minIntMinAx, -0.2, -0.1, -0.05];
        %now KK
        minKKValue = min([LQCtrlNormKK, LQTreatNormKK], [], 'all');
        maxKKValue = max([UQCtrlNormKK, UQTreatNormKK], [], 'all');
        minKKMinAx = min([minKKCtrlTick, minKKTreatTick], [], 'all');
        maxKKMinAx = max([minKKCtrlTick, minKKTreatTick], [], 'all');
        maxKKMaxAx = max([maxKKCtrlTick, maxKKTreatTick], [], 'all');
        minKKMaxAx = min([maxKKCtrlTick, maxKKTreatTick], [], 'all');
    else
        minIntValue = min(iCtrlNormLQ);
        possIntVec = [maxIntCtrlTick, -0.5, -0.2, -0.1, -0.05];
    end
    
    lowerIntLimit = [];
    possIntVec = sort(possIntVec);
    if minIntValue < possIntVec(2)
        lowerIntLimit = possIntVec(1);
    elseif minIntValue < possIntVec(3)
        lowerIntLimit = possIntVec(2);
    elseif minIntValue < possIntVec(4)
        lowerIntLimit = possIntVec(3);
    elseif minIntValue < possIntVec(5)
        lowerIntLimit = possIntVec(4);
    else
        lowerIntLimit = possIntVec(5);
    end
    
    figure(ctrlFig);
    yyaxis left
    set(ctrlax, 'YLim', [lowerIntLimit maxIntCtrlTick]);
    
    yyaxis right
    %only change KK axis if there is treatment data
    if isTreat
        if minKKValue < maxKKMinAx
            lowerKKAx = minKKMinAx; %lower KK bound smallest value
        else
            lowerKKAx = maxKKMinAx; %lower KK bound second smallest value
        end
        if maxKKValue > minKKMaxAx
            upperKKAx = maxKKMaxAx; %upper KK bound highest value
        else
            upperKKAx = minKKMaxAx; %upper KK bound second highest value
        end
        %the two Mad2 axes will probably be scaled to different values so
        %I also have to scale the KK axes for the treatment/control graphs.
        %Easiest/best way to do this when Mad2 treatment < Mad2 ctrl is to
        %stretch Mad2 axis on treat. If Mad2 treatment > Mad2 ctrl,
        %increase upper bound of treatment KK.
        set(ctrlax, 'YLim', [lowerKKAx, upperKKAx]);
        ctrlIntRange = maxIntCtrlTick - lowerIntLimit;
        treatIntRange = maxIntTreatTick - lowerIntLimit;
        if ctrlIntRange > treatIntRange
            maxIntTreatTick = maxIntCtrlTick;
            upperTreatKK = upperKKAx;
        else
            treatCtrlAxRatio = treatIntRange/ctrlIntRange;
            KKRange = upperKKAx - lowerKKAx;
            treatKKRange = KKRange*treatCtrlAxRatio;
            upperTreatKK = lowerKKAx + treatKKRange;
        end
    end
    ctrlIntRange = maxIntCtrlTick - lowerIntLimit;
    set(ctrlax, 'Units', 'pixels')
    ctrlPos = get(ctrlax, 'Position');
    set(ctrlax, 'Position', [ctrlPos(1), ctrlPos(2)*1.25, ctrlPos(3), ctrlIntRange*200])
    set(ctrlFig, 'Position', [100 60 ctrlPos(3)+(ctrlPos(1)*2) (ctrlIntRange*200)+(ctrlPos(2)*2)])
    set(ctrlax, 'XTickLabelRotation', 0)
    
    if isTreat
        figure(treatFig);
        yyaxis left
        set(treatax, 'YLim', [lowerIntLimit maxIntTreatTick]);
        yyaxis right
        set(treatax, 'YLim', [lowerKKAx upperTreatKK]);
        treatIntRange = maxIntTreatTick - lowerIntLimit;
        set(treatax, 'Units', 'pixels')
        treatPos = get(treatax, 'Position');
        set(treatax, 'Position', [treatPos(1), treatPos(2)*1.25, treatPos(3), treatIntRange*200])
        set(treatFig, 'Position', [700 60 treatPos(3)+(treatPos(1)*2) (treatIntRange*200)+(treatPos(2)*2)])
        set(treatax, 'XTickLabelRotation', 0)
    end
end


end



%% getKKBins_withNormInt function
function [KKbins, nCellsInExpt] = getKKBins_withNormInt(iM, intiM, scaledBinnedintiMs, iMax, nBins, prevKKbins, indIntNorm)
if isempty(prevKKbins)
    KKbins = cell(1, nBins);
else
    KKbins = prevKKbins;
end
cellLabel = iM.label;
nPairsKTs = length(cellLabel);
intOrg = intiM.intensity.mean.outer;
if indIntNorm
    intOrg = intOrg ./ intiM.intensity.mean.inner;
end
intOrgScaled = intOrg/iMax;
CellIDs = [];
for iPair = 1:nPairsKTs
    if ~isnan(intOrgScaled(iPair, 1)) && ~isnan(intOrgScaled(iPair, 2))
        % find column (bin) and row (position in bin) in rebinned intiM of
        % both KTs in pair
        [KT1binR, KT1binC] = find(scaledBinnedintiMs==intOrgScaled(iPair,1));
        [KT2binR, KT2binC] = find(scaledBinnedintiMs==intOrgScaled(iPair,2));
        if ~isempty(KT1binC) && ~isempty(KT2binC)
            % continue if both KTs in pair have been accepted for further analysis
            if length(KT1binC) > 1
                % in case that intensity of KT = 1 (thanks to
                % normalisation), just select first intensity position
                % (will be maximum bin anyway)
                KT1binC = KT1binC(1);
                KT1binR = KT1binR(1);
            end
            if length(KT2binC) > 1
                % as above
                KT2binC = KT2binC(1);
                KT2binR = KT2binR(1);
            end
            if (KT1binC == KT2binC) || (KT1binC+1 == KT2binC) || (KT1binC == KT2binC+1)
                % if both KTs in pair are in same bin or consecutive bins,
                % then get KK distance
                workingKK = iM.microscope.sisSep.threeD(iPair);
                if workingKK > 2
                    % skip if KK > 2, probably paired incorrectly
                    continue
                end
                % get idx of cell with pair
                CellIDs = cat(2, CellIDs, str2num(cellLabel(iPair,3:4)));
                if KT1binC == KT2binC
                    % if both KTs in same bin, put KK value in their bin
                    binID = KT1binC;
                elseif (KT1binC+1 == KT2binC) || (KT1binC == KT2binC+1)
                    % otherwise, find where the "centre" of both KT
                    % positions are - i.e. if KT1 is bin 1, row 151 (of
                    % 200) and KT2 is bin 2, row 101 (of 200), the "point"
                    % position is 1.75 for KT1 and 2.5 for KT2, if you
                    % average those you get 2.125. therefore push KK value
                    % to bin 2
                    KT1val = KT1binC+((KT1binR-1)/size(scaledBinnedintiMs,1));
                    KT2val = KT2binC+((KT2binR-1)/size(scaledBinnedintiMs,1));
                    avKTval = mean([KT1val, KT2val]);
                    threshBin = max([KT1binC, KT2binC]);
                    % if "average bin" below threshold, push to lower bin,
                    % if above threshold, push to upper bin, if equal, push
                    % to lower bin.
                    if avKTval < threshBin
                        binID = min([KT1binC, KT2binC]);
                    elseif avKTval > threshBin
                        binID = max([KT1binC, KT2binC]);
                    else
                        binID = max([KT1binC, KT2binC]); %20250429
                    end
                end
                % add KK value to relevant bin in {cell} structure
                if ~isempty(KKbins{binID})
                    workingData = KKbins{binID};
                    pasteData = cat(1, workingData, workingKK);
                    KKbins{binID} = pasteData;
                else
                    KKbins{binID} = workingKK;
                end
                
            end
        end

    end
end
CellIDs = unique(CellIDs);
nCellsInExpt = length(CellIDs); %to print nCells that have at least one pair assessed in each iExpt

end


%% getKKquartiles function
function [MedianKK, lowerKK, upperKK, averageInBin, pairsInBins] = getKKquartiles(KKcell)

MedianKK = [];
lowerKK = [];
upperKK = [];
pairsInBins = [];

for iBin = 1:length(KKcell)
    CellMedian = median(KKcell{iBin});
    CellLQ = quantile(KKcell{iBin}, 0.25);
    CellUQ = quantile(KKcell{iBin}, 0.75);
    pairsThisBin = [iBin length(KKcell{iBin})]; %to print nPairs per bin
    pairsInBins = cat(1, pairsInBins, pairsThisBin(2));
    MedianKK = cat(2, MedianKK, CellMedian);
    lowerKK = cat(2, lowerKK, CellLQ);
    upperKK = cat(2, upperKK, CellUQ);
end

averageInBin = round(mean(pairsInBins));

end
%% getKKConfInts function, for confidence intervals 
function [MedianKK, lowerKK, upperKK, averageInBin, pairsInBins] = getKKConfInts(KKcell)

MedianKK = [];
lowerKK = [];
upperKK = [];
pairsInBins = [];

for iBin = 1:length(KKcell)
    CellMedian = median(KKcell{iBin});
    nVals = length(KKcell{iBin});
    lPos = (nVals*0.5) - (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101 - first two 0.5 values indicate we are assessing median quantile, third 0.5 is from 1-0.5 (probability of not 0.5).
    uPos = (nVals*0.5) + (1.96*sqrt(nVals*0.5*0.5)); %from doi 10.11613/BM.2019.010101
    lPct = lPos/nVals;
    uPct = uPos/nVals;
    lowerCI = quantile(KKcell{iBin}, lPct);
    upperCI = quantile(KKcell{iBin}, uPct);
    pairsThisBin = [iBin nVals]; %to print nPairs per bin
    pairsInBins = cat(1, pairsInBins, pairsThisBin(2));
    MedianKK = cat(2, MedianKK, CellMedian);
    lowerKK = cat(2, lowerKK, lowerCI);
    upperKK = cat(2, upperKK, upperCI);
end

averageInBin = round(mean(pairsInBins));
avBin = sum(pairsInBins);

end
%% makeQuartiles function
function [medianVal, lowerQuart, upperQuart] = makeQuartiles(binnedMatrix)
medianVal = median(binnedMatrix);
lowerQuart = quantile(binnedMatrix, 0.25);
upperQuart = quantile(binnedMatrix, 0.75);
end

%% makeMedianConfInts function, for confidence intervals
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
%% processOptions function

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

%% makeBins function
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
function nCells = getnCells(intiMs)
nExpts = size(intiMs);
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