function [medianVal, lowerCI, upperCI, lowerQuart, upperQuart] = makeMedianConfIntsQuartiles(binnedMatrix)
%
% gets the median, quartiles, and 95% confidence intervals of a binned
% matrix (usually intensity)
%
% Copyright (c) 2024 C. C. Conway
%

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