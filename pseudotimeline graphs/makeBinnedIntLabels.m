function [IndIntsBinned, LabelsBinned] = makeBinnedIntLabels(BinByIntiM, intiMlabels, numBins, iExpt, swapBin)
% MAKEBINNEDINTLABELS is a function to bin kinetochore labels present in an
% intiM structure based on the intensity of a given marker.
%
%
%   MAKEBINNEDINTLABELS(BinByIntiM, intiMlabels, numBins, iExpt, swapBin)
%     - BinByIntiM is actually a matrix of intensities that should be the
%       same length as intiMlabels
%     - intiMLabels should be a vector of labels originating from an intiM
%       structure
%     - numBins is how many bins you want to split your labels into
%     - iExpt is the number that the first two numbers of the label will be
%       changed to permit experiment concatenation in other functions
%     - swapBin is whether bin 1 or bin 25 should correspond to the maximum
%       intensity dataset
%
%
% Copyright (c) 2026 C. C. Conway

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
end