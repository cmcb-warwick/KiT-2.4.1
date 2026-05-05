function [IndIntsBinned, DepIntsBinned] = binIntiMs(indIntiM, depIntiM, varargin)
% Function to enable binning of two intiMs. The order of intiM binning
% will be shared across vectors.
%
% Copyright (c) 2023 C. C. Conway
%
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
    
    [minIndSubset, minIdx] = mink(indIntsCol, n); %gets lowest n values and their positions (ordered lowest to highest)
    [maxIndSubset, maxIdx] = maxk(minIndSubset, k); %gets highest k values from lowest n values and their positions (ordered lowest to highest)
    
    minDepSubset = depIntsCol(minIdx); %gets corresponding depInt values from lowest n indInt values
    maxDepSubset = minDepSubset(maxIdx); %gets corresponding depInt values from highest k of lowest n indInt values

    IndIntsBinned(:, bin) = maxIndSubset;
    DepIntsBinned(:, bin) = maxDepSubset;
    
    n = n + k;
end

if opts.minIntFirstCol
    IndIntsBinned = flipud(IndIntsBinned); %if first column is min, then top left should be min value
    DepIntsBinned = flipud(DepIntsBinned); %flips correspondingly
else
    IndIntsBinned = fliplr(IndIntsBinned); %if last column is min, flip matrix accordingly
    DepIntsBinned = fliplr(DepIntsBinned); %flips correspondingly
end


end


%% processOptions
function options=processOptions(defaults,varargin)
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