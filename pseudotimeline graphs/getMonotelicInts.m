function MonotelicInts = getMonotelicInts(iMintiMSets, ctrlIntiMs, varargin)
% GETMONOTELICINTS Produces a structure with statistics of second channel
% intensities from monotelic pairs (one low Mad2 KT, one high Mad2 KT).
% 
%   GETMONOTELICINTS(iMintiMSets, dividingIntiMs)
%   iMintiMSets must be organised as {iM, intiMInd, intiMDep; iM2, 
%   intiM2Ind, intiM2Dep} and so on. ctrlIntiMs must be organised as
%   {unpairedintiMInd, unpairedintiMDep; unpairedintiMInd2,
%   unpairedintiMDep2} and so on. Unpaired data should be from DMSO repeat
%   of data. If only paired data is available, use same intiMs as used in
%   iMintiMSets. Must have same number of expts in both iMintiMSets and
%   ctrlIntiMs.
%   
%   Options are available:
%
%   Options, defaults in {}:-
%
%   nBins: {25} or other integer. The number of bins to split each intiM
%       structure into for normalisation.
%
%   normIndInt: 0 or {1}. Whether to normalise the independent intensity
%       marker (i.e. Mad2) to CenpC or not.
%
%   normDepInt: 0 or {1}. Whether to normalise the dependent intensity
%       marker to CenpC or not.
%
%   normalise: 0 or {1}. Whether to normalise the median of the maximum bin
%       for each intensity class to one. 
%
%   minHighMad2Int: {0.4} or another value. The minimum normalised
%       intensity a KT must have to be considered high Mad2.
%
%   maxLowMad2Int: {0.01} or another value. The maximum normalised
%       intensity a KT must have to be considered low Mad2.
%
% Copyright (c) 2025 C. C. Conway

opts.nBins = 25;
opts.normIndInt = 1;
opts.normDepInt = 1;
opts.normalise = 1;
opts.minHighMad2Int = 0.4; %minimum KT value to count as Mad2 +ve
opts.maxLowMad2Int = 0.01; %maximum KT value to count as Mad2 -ve
opts = processOptions(opts,varargin{:});

nBins = opts.nBins;
isNorm = opts.normalise;
nExpts = size(iMintiMSets, 1);


nCells = 0;
ctrlIndIntiMScaled = [];
ctrlDepIntiMScaled = [];
allciMax = [];
allcdMax = [];

for iExpt = 1:nExpts
    ctrlIntiMs{iExpt,1}.intensity.mean.inner = rmLowCenpC(ctrlIntiMs{iExpt,1});
    ctrlIntiMs{iExpt,2}.intensity.mean.inner = rmLowCenpC(ctrlIntiMs{iExpt,2});
   
    cIndIntiM = ctrlIntiMs{iExpt,1}.intensity.mean.outer;
    if opts.normIndInt
        cIndIntiM = cIndIntiM ./ ctrlIntiMs{iExpt,1}.intensity.mean.inner;
    end

    cDepIntiM = ctrlIntiMs{iExpt,2}.intensity.mean.outer;
    if opts.normDepInt
        cDepIntiM = cDepIntiM ./ ctrlIntiMs{iExpt,2}.intensity.mean.inner;
    end

    cIndIntiM = rmmissing(cIndIntiM(:));
    cDepIntiM = rmmissing(cDepIntiM(:));
    [ctrlIndIntiM, ctrlDepIntiM] = rebinIntensities(cIndIntiM, cDepIntiM, 'nBins', nBins);

    ciMaxMed = max(median(ctrlIndIntiM));
    cdMaxMed = max(median(ctrlDepIntiM));
    
    if isNorm
        allciMax = cat(2, allciMax, ciMaxMed);
        wkCtrlIndIntiM = ctrlIndIntiM/ciMaxMed;
        ctrlIndIntiMScaled = cat(1, ctrlIndIntiMScaled, wkCtrlIndIntiM);
        
        allcdMax = cat(2, allcdMax, cdMaxMed);
        wkCtrlDepIntiM = ctrlDepIntiM/cdMaxMed;
        ctrlDepIntiMScaled = cat(1, ctrlDepIntiMScaled, wkCtrlDepIntiM);

    end
end

[ctrlIndIntiMScaledBinned, ctrlDepIntiMScaledBinned] = rebinIntensities(ctrlIndIntiMScaled, ctrlDepIntiMScaled, 'nBins', nBins);
intData = [];

ciNormMax = max(median(ctrlIndIntiMScaledBinned));
cdNormMax = max(median(ctrlDepIntiMScaledBinned));

if isNorm
    allciMax = allciMax/ciNormMax;
    allcdMax = allcdMax/cdNormMax;
end


for iExpt = 1:nExpts
    iMintiMSets{iExpt,2}.intensity.mean.inner = rmLowCenpC(iMintiMSets{iExpt,2});
    iMintiMSets{iExpt,3}.intensity.mean.inner = rmLowCenpC(iMintiMSets{iExpt,3});

    [intData, wknCells] = getIntData(iMintiMSets{iExpt,1}, iMintiMSets{iExpt,2}, iMintiMSets{iExpt,3}, ...
        allciMax(iExpt), allcdMax(iExpt), intData, opts.normIndInt, opts.normDepInt, opts.minHighMad2Int, opts.maxLowMad2Int);

        nCells = nCells + wknCells;
end
MonotelicInts.intData = intData;
MonotelicInts.nCells = nCells;



end



%% getIntData function
function [intData, nCellsInExpt] = getIntData(iM, indIntiM, depIntiM, iMax, dMax, prevData, indIntNorm, depIntNorm, minHighMad2, maxLowMad2)
if isempty(prevData)
    intData = [];
else
    intData = prevData;
end
cellLabel = iM.label;
nPairsKTs = length(cellLabel);
iIntOrg = indIntiM.intensity.mean.outer;
if indIntNorm
    iIntOrg = iIntOrg ./ indIntiM.intensity.mean.inner;
end
iIntOrgScaled = iIntOrg/iMax;

dIntOrg = depIntiM.intensity.mean.outer;
if depIntNorm
    dIntOrg = dIntOrg ./ depIntiM.intensity.mean.inner;
end
dIntOrgScaled = dIntOrg/dMax;
CellIDs = [];
Mad2PosDepInts = [];
Mad2NegDepInts = [];
for iPair = 1:nPairsKTs
    if ~isnan(iIntOrgScaled(iPair, 1)) && ~isnan(iIntOrgScaled(iPair, 2))
        % find Mad2 intensity in divided intiM of
        % both KTs in pair
        KT1Int = iIntOrgScaled(iPair, 1);
        KT2Int = iIntOrgScaled(iPair, 2);
        if KT1Int >= minHighMad2
            KT1Marker = -1; %if KT is high intensity give it value -1
        elseif KT1Int <= maxLowMad2
            KT1Marker = 1; %if KT is low intensity give it value 1
        else
            KT1Marker = 2; %if KT is in between high and low give it value 2
        end

        if KT2Int >= minHighMad2
            KT2Marker = -1; %if KT is high intensity give it value -1
        elseif KT2Int <= maxLowMad2
            KT2Marker = 1; %if KT is low intensity give it value 1
        else
            KT2Marker = 2; %if KT is in between high and low give it value 2
        end
        
        KTPairVal = KT1Marker + KT2Marker; %if we have high and low KT we will have value of 0
        if KTPairVal == 0
            workingKK = iM.microscope.sisSep.threeD(iPair);
            if workingKK > 2
                % skip if KK > 2, probably paired incorrectly
                continue
            end
            if KT1Marker == -1
                hKT = 1; %high Mad2 KT
                lKT = 2; %low Mad2 KT
            else
                hKT = 2;
                lKT = 1;
            end
            wkPos = dIntOrgScaled(iPair, hKT);
            wkNeg = dIntOrgScaled(iPair, lKT);
            %Mad2PosDepInts = cat(1, Mad2PosDepInts, wkPos);
            %Mad2NegDepInts = cat(1, Mad2NegDepInts, wkNeg);
            %get cell ID
            CellIDs = cat(2, CellIDs, str2num(cellLabel(iPair,3:4)));
            % add KK value to structure
            if ~isempty(intData)
                workingKKData = intData.KKdist;
                pasteData = cat(1, workingKKData, workingKK);
                intData.KKdist = pasteData;
                wkHMad2Data = intData.Mad2PosInts;
                wkHMad2Data = cat(1, wkHMad2Data, wkPos);
                intData.Mad2PosInts = wkHMad2Data;
                wkLMad2Data = intData.Mad2NegInts;
                wkLMad2Data = cat(1, wkLMad2Data, wkNeg);
                intData.Mad2NegInts = wkLMad2Data;
            else
                intData.KKdist = workingKK;
                intData.Mad2PosInts = wkPos;
                intData.Mad2NegInts = wkNeg;
            end
        end
          

    end
end
CellIDs = unique(CellIDs);
nCellsInExpt = length(CellIDs); %to print nCells that have at least one pair assessed in each iExpt

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


%% rmLowCenpC function
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

%% rebinIntensities function
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