function halfChanges = getHalfChanges(orgData, varargin)
%I want half change data to also tell me the maximum and minimum values and
%the ANOVAs of "similar" data (i.e. if max bin is bin 3 but ANOVA says bin
%1 through 6 are similar to 3 then I will record that data. I will also
%record the median value of the max bins and min bins. I will then find the
%halfway point, and calculate which bin position that is, and record. I
%could probably also create vertical lines to show where the half change
%is if I wanted to integrate into my graph plotting.
%
%
% Copyright (c) 2025 C. C. Conway



opts.isCtrl = 1; %other option is 0 for 'treatment'
opts.dataType = 'int'; %other options are 'kk' or 'delta'
opts = processOptions(opts,varargin{:});

dataType = opts.dataType;
switch dataType
    case 'int'
        sampleData = 1;
        if opts.isCtrl
            Mad2data = orgData.control.independent.intensities;
            OtherData = orgData.control.dependent.intensities;
        else
            Mad2data = orgData.treatment.independent.intensities;
            OtherData = orgData.treatment.dependent.intensities;
        end
        [~, ~, anovaOther] = anova1(OtherData, [], 'off');
    case 'kk'
        sampleData = 1;
        if opts.isCtrl
            Mad2data = orgData.control.intensities;
            OtherData = [];
            OtherDataGroups = [];
            for iBin = 1:length(orgData.control.KK)
                OtherData = cat(1,OtherData, orgData.control.KK{iBin});
                OtherDataGroups = cat(1, OtherDataGroups, ones([length(orgData.control.KK{iBin}),1])*iBin);
            end
        else
            Mad2data = orgData.treatment.intensities;
            OtherData = [];
            OtherDataGroups = [];
            for iBin = 1:length(orgData.treatment.KK)
                OtherData = cat(1,OtherData, orgData.treatment.KK{iBin});
                OtherDataGroups = cat(1, OtherDataGroups, ones([length(orgData.treatment.KK{iBin}),1])*iBin);
            end
        end
        [~, ~, anovaOther] = anova1(OtherData, OtherDataGroups, 'off');
    case 'delta'
        sampleData = 0;
        Mad2data = orgData.ints.data;
        % otherData will be ignored for now - going to have to not do
        % the multCompare method and will instead have to do a sample
        % draw method
end

[~, ~, anovaMad2] = anova1(Mad2data, [], 'off');

%% Mad2 data
[maxMad2, maxMad2Bin] = max(anovaMad2.means);
[minMad2, minMad2Bin] = min(anovaMad2.means);

cMad2 = multcompare(anovaMad2, 'Display', 'off');
maxPValsMad2 = [];
minPValsMad2 = [];

nRowsMad2 = size(cMad2, 1);

for iRow = 1:nRowsMad2
    pos1 = cMad2(iRow, 1);
    pos2 = cMad2(iRow, 2);
    if pos1 == maxMad2Bin
        wkMax = [pos1, pos2, cMad2(iRow, 6)]; %gives p val
        maxPValsMad2 = cat(1, maxPValsMad2, wkMax);
    elseif pos1 == minMad2Bin
        wkMin = [pos1, pos2, cMad2(iRow, 6)]; %gives p val
        minPValsMad2 = cat(1, minPValsMad2, wkMin);
    end
    if pos2 == maxMad2Bin
        wkMax = [pos2, pos1, cMad2(iRow, 6)]; %gives p val
        maxPValsMad2 = cat(1, maxPValsMad2, wkMax);
    elseif pos2 == minMad2Bin
        wkMin = [pos2, pos1, cMad2(iRow, 6)]; %gives p val
        minPValsMad2 = cat(1, minPValsMad2, wkMin);
    end
end

nRowsMad2 = size(maxPValsMad2, 1);

maxStatSimilarMad2 = maxMad2Bin;
minStatSimilarMad2 = minMad2Bin;


for iRow = 1:nRowsMad2
    if maxPValsMad2(iRow, 3) >= 0.05
        maxStatSimilarMad2 = cat(1, maxStatSimilarMad2, maxPValsMad2(iRow, 2));
    end
    if minPValsMad2(iRow, 3) >= 0.05
        minStatSimilarMad2 = cat(1, minStatSimilarMad2, minPValsMad2(iRow, 2));
    end
end

nMaxColsMad2 = length(maxStatSimilarMad2);
nMinColsMad2 = length(minStatSimilarMad2);

allMaxValsMad2 = [];
allMinValsMad2 = [];

for iCol = 1:nMaxColsMad2
    colID = maxStatSimilarMad2(iCol);
    wkData = Mad2data(:,colID);
    allMaxValsMad2 = cat(1, allMaxValsMad2, wkData);
end

for iCol = 1:nMinColsMad2
    colID = minStatSimilarMad2(iCol);
    wkData = Mad2data(:,colID);
    allMinValsMad2 = cat(1, allMinValsMad2, wkData);
end

maxMedianMad2 = median(allMaxValsMad2);
minMedianMad2 = median(allMinValsMad2);

halfDiffMad2 = (maxMedianMad2-minMedianMad2)/2;
halfChangeMad2 = minMedianMad2+halfDiffMad2;

Mad2Medians = median(Mad2data);

betweenPosMad2 = [];
nComps = length(Mad2Medians)-1;
for iComp = 1:nComps
    if halfChangeMad2 <= Mad2Medians(iComp) && Mad2Medians(iComp+1) <= halfChangeMad2
        % if col1 > halfChange > col2
        %calculate position from y=mx+c slope
        slope = Mad2Medians(iComp+1)-Mad2Medians(iComp); %we know x change is 1 so no need to divide
        mxC = Mad2Medians(iComp)-(iComp*slope);
        mXc = (halfChangeMad2-mxC)/slope;
        betweenPosMad2 = cat(1, betweenPosMad2, mXc);
    elseif halfChangeMad2 <= Mad2Medians(iComp+1) && Mad2Medians(iComp) <= halfChangeMad2
        % if col2 > halfChange > col1
        slope = Mad2Medians(iComp+1)-Mad2Medians(iComp); %we know x change is 1 so no need to divide
        mxC = Mad2Medians(iComp)-(iComp*slope);
        mXc = (halfChangeMad2-mxC)/slope;
        betweenPosMad2 = cat(1, betweenPosMad2, mXc);
    end
end

meanBetweenMad2 = mean(betweenPosMad2);

halfChanges.Mad2.HCintensity = halfChangeMad2;
halfChanges.Mad2.HCbin = meanBetweenMad2;


%% other data
if sampleData
    [maxOther, maxOtherBin] = max(anovaOther.means);
    [minOther, minOtherBin] = min(anovaOther.means);
    
    cOther = multcompare(anovaOther, 'Display', 'off');
    maxPValsOther = []; %to collate bins that are statistically "similar" to the maximum mean intensity/KK/Delta
    minPValsOther = []; %to collate bins that are statistically "similar" to the minimum mean intensity/KK/Delta
    
    nRowsOther = size(cOther, 1);
    
    for iRow = 1:nRowsOther
        pos1 = cOther(iRow, 1);
        pos2 = cOther(iRow, 2);
        if pos1 == maxOtherBin
            wkMax = [pos1, pos2, cOther(iRow, 6)]; %gives p val
            maxPValsOther = cat(1, maxPValsOther, wkMax);
        elseif pos1 == minOtherBin
            wkMin = [pos1, pos2, cOther(iRow, 6)]; %gives p val
            minPValsOther = cat(1, minPValsOther, wkMin);
        end
        if pos2 == maxOtherBin
            wkMax = [pos2, pos1, cOther(iRow, 6)]; %gives p val
            maxPValsOther = cat(1, maxPValsOther, wkMax);
        elseif pos2 == minOtherBin
            wkMin = [pos2, pos1, cOther(iRow, 6)]; %gives p val
            minPValsOther = cat(1, minPValsOther, wkMin);
        end
    end
    
    nRowsOther = size(maxPValsOther, 1); %this should be nBins-1
    
    maxStatSimilarOther = maxOtherBin;
    minStatSimilarOther = minOtherBin;
    
    
    for iRow = 1:nRowsOther
        if maxPValsOther(iRow, 3) >= 0.05
            maxStatSimilarOther = cat(1, maxStatSimilarOther, maxPValsOther(iRow, 2));
        end
        if minPValsOther(iRow, 3) >= 0.05
            minStatSimilarOther = cat(1, minStatSimilarOther, minPValsOther(iRow, 2));
        end
    end
    
    nMaxColsOther = length(maxStatSimilarOther);
    nMinColsOther = length(minStatSimilarOther);
    
    allMaxValsOther = [];
    allMinValsOther = [];
    
    for iCol = 1:nMaxColsOther
        colID = maxStatSimilarOther(iCol);
        switch dataType
            case 'kk'
                if opts.isCtrl
                    wkData = orgData.control.KK{colID};
                else
                    wkData = orgData.treatment.KK{colID};
                end
            case 'int'
                wkData = OtherData(:,colID);
        end
        allMaxValsOther = cat(1, allMaxValsOther, wkData);
    end
    
    for iCol = 1:nMinColsOther
        colID = minStatSimilarOther(iCol);
        switch dataType
            case 'kk'
                if opts.isCtrl
                    wkData = orgData.control.KK{colID};
                else
                    wkData = orgData.treatment.KK{colID};
                end
            case 'int'
                wkData = OtherData(:,colID);
        end
        allMinValsOther = cat(1, allMinValsOther, wkData);
    end
    
    maxMedianMeanOther = median(allMaxValsOther); %both intensity and KK data are plotted directly from real data and medians of these are plotted
    minMedianMeanOther = median(allMinValsOther);


else
    % first have to do a multi comparison between all bins with estimated values based on mean and std dev. Correct using sample size.

    contingencyTable = [];
    n = orgData.binSize;
    for i = 1:24
        for j = (i+1):25
            % following was generated by copilot

            % Define group statistics
            mean1 = orgData.delta.mean(i);
            std1 = orgData.delta.std(i);
            
            mean2 = orgData.delta.mean(j);
            std2 = orgData.delta.std(j);
            
            % Calculate pooled standard error
            se1 = std1^2 / n;
            se2 = std2^2 / n;
            se_pooled = sqrt(se1 + se2);
            
            % Calculate t-statistic
            t_stat = (mean1 - mean2) / se_pooled;
            
            % Degrees of freedom (Welch-Satterthwaite approximation)
            df_num = (se1 + se2)^2;
            df_den = (se1^2 / (n - 1)) + (se2^2 / (n - 1));
            df = df_num / df_den;
            
            % Calculate two-tailed p-value
            p_value = 2 * tcdf(-abs(t_stat), df);
            
            % copilot no longer making this code
            % this is bonferroni correction
            % nchoosek is number of tests I've made (300 for pairwise between all 25 bins)
            nRowsOther = nchoosek(25,2);
            if p_value < 0.05/nRowsOther
                significant = 1;
            else
                significant = 0;
            end
            currentData = [i j p_value significant];
            contingencyTable = cat(1, contingencyTable, currentData);
        end
    
    end

    [maxOther, maxOtherBin] = max(orgData.delta.mean);
    [minOther, minOtherBin] = min(orgData.delta.mean);

    maxBinsOther = [maxOtherBin, maxOther]; %to collate bins that are statistically "similar" to the maximum mean Delta
    minBinsOther = [minOtherBin, minOther]; %to collate bins that are statistically "similar" to the minimum mean Delta
    
    
    for iRow = 1:nRowsOther
        pos1 = contingencyTable(iRow, 1);
        pos2 = contingencyTable(iRow, 2);
        isSig = contingencyTable(iRow, 4);
        if pos1 == maxOtherBin
            if ~isSig
                meanOther = orgData.delta.mean(pos2);
                maxBinsOther = cat(1, maxBinsOther, [pos2, meanOther]);
            end
        elseif pos1 == minOtherBin
            if ~isSig
                meanOther = orgData.delta.mean(pos2);
                minBinsOther = cat(1, minBinsOther, [pos2, meanOther]);
            end
        end
        if pos2 == maxOtherBin
            if ~isSig
                meanOther = orgData.delta.mean(pos1);
                maxBinsOther = cat(1, maxBinsOther, [pos1, meanOther]);
            end
        elseif pos2 == minOtherBin
            if ~isSig
                meanOther = orgData.delta.mean(pos1);
                minBinsOther = cat(1, minBinsOther, [pos1, meanOther]);
            end
        end
    end

    maxMedianMeanOther = mean(maxBinsOther(:,2)); %this is because Delta values are all plotted based on mean values and has no "real data" to draw from, more like fitting a model
    minMedianMeanOther = mean(minBinsOther(:,2));
end


halfDiffOther = (maxMedianMeanOther-minMedianMeanOther)/2;
halfChangeOther = minMedianMeanOther+halfDiffOther;
switch dataType
    case 'kk'
        if opts.isCtrl
            OtherMeanMedians = orgData.control.KKStats.median;
        else
            OtherMeanMedians = orgData.treatment.KKStats.median;
        end
    case 'int'
        OtherMeanMedians = median(OtherData);
    case 'delta'
        OtherMeanMedians = orgData.delta.mean;
end


betweenPosOther = [];
nComps = length(OtherMeanMedians)-1;
for iComp = 1:nComps
    if halfChangeOther <= OtherMeanMedians(iComp) && OtherMeanMedians(iComp+1) <= halfChangeOther
        % if col1 > halfChange > col2
        %calculate position from y=mx+c slope
        slope = OtherMeanMedians(iComp+1)-OtherMeanMedians(iComp); %we know x change is 1 so no need to divide
        mxC = OtherMeanMedians(iComp)-(iComp*slope);
        mXc = (halfChangeOther-mxC)/slope;
        betweenPosOther = cat(1, betweenPosOther, mXc);
    elseif halfChangeOther <= OtherMeanMedians(iComp+1) && OtherMeanMedians(iComp) <= halfChangeOther
        % if col2 > halfChange > col1
        slope = OtherMeanMedians(iComp+1)-OtherMeanMedians(iComp); %we know x change is 1 so no need to divide
        mxC = OtherMeanMedians(iComp)-(iComp*slope);
        mXc = (halfChangeOther-mxC)/slope;
        betweenPosOther = cat(1, betweenPosOther, mXc);
    end
end

meanBetweenOther = mean(betweenPosOther); %just in case there are multiple "up-down-up" sections so the half-change is traversed multiple times

halfChanges.Other.HCintensityKKDelta = halfChangeOther;
halfChanges.Other.HCbin = meanBetweenOther;

%now want to get Mad2 intensity on same half change
lowMad2Bin = ceil(meanBetweenOther); %because higher bin number = lower Mad2 int
highMad2Bin = floor(meanBetweenOther); %because lower bin number = higher Mad2 int
lowMad2Val = Mad2Medians(lowMad2Bin);
highMad2Val = Mad2Medians(highMad2Bin);

slope = lowMad2Val-highMad2Val; %m = dy/dx (dx = lowMad2Bin-highMad2Bin = 1 therefore can simplify)
constant = lowMad2Val-(slope*lowMad2Bin); %c=y-mx
Mad2Y = (slope*meanBetweenOther)+constant; %y=mx+c

halfChanges.Other.HCMad2Int = Mad2Y;

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