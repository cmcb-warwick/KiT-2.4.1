function [rebinnedIndInts, rebinnedDepInts] = rebinIntensities(IndIntsOrg, DepIntsOrg, varargin)
%
% Function to take binned normalised data (independent and dependent, i.e.
% Mad2 and other marker) from multiple experiments and re-bin on the
% independent datasets.
%
% Copyright (c) 2025 C. C. Conway
%
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
    
    [minIndSubset, minIdx] = mink(indIntsCol, n); %gets lowest n values and their positions (ordered lowest to highest)
    [maxIndSubset, maxIdx] = maxk(minIndSubset, k); %gets highest k values from lowest n values and their positions (ordered lowest to highest)
    
    minDepSubset = depIntsCol(minIdx); %gets corresponding depInt values from lowest n indInt values
    maxDepSubset = minDepSubset(maxIdx); %gets corresponding depInt values from highest k of lowest n indInt values

    IndIntsBinned(:, bin) = maxIndSubset;
    DepIntsBinned(:, bin) = maxDepSubset;
    
    n = n + k;
end

if opts.minIntFirstCol
    rebinnedIndInts = flipud(IndIntsBinned); %if first column is min, then top left should be min value
    rebinnedDepInts = flipud(DepIntsBinned); %flips correspondingly
else
    rebinnedIndInts = fliplr(IndIntsBinned); %if last column is min, flip matrix accordingly
    rebinnedDepInts = fliplr(DepIntsBinned); %flips correspondingly
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