function dublShowMultiKTs(job,varargin)
% DUBLSHOWMULTIKTS Plots image of dual or triple-channel movie with coordinates 
% for a given kinetochore.
%
%    DUBLSHOWMULTIKTS(JOB,...) Plots coordinates in one, two or three channels
%    over a still image for a given kinetochore at a given timepoint.
%
%    Options, defaults in {}:-
%
%    chanOrder: {[2 1 3]} or some perturbation. Order in which the titled
%    channels are presented in the figure, where typically:
%    1=red, 2=green, 3=blue.
%
%    contrast: {[0.1 0.9995] [0.1 0.9995] [0.1 1]}, 'help', or similar. 
%       Upper and lower contrast limits for each channel. Values must be 
%       [0 1]. 'help' outputs the minimum and maximum intensity values as 
%       in range guidance, then requests values from the user.
%       Tips: - Increase the lower limit to remove background noise.
%             - Decrease the upper limit to increase brightness.
%
%    newFig: 0 or {1}. Whether or not to show the kinetochore in a new
%           figure.
%
%    plotChannels: {[1 2 3]} or some subset of [1 2 3]. Vector of channels
%           for plotting.
%
%    plotTitles: {'Mad2','CenpC','SKAP'} or similar. Titles of
%           each plot denoting which channel is which.
%
%    kinetochore: {1} or vector of integers. Kinetochore(s) within JOB
%           being plotted.
%
%    subpixelate: {9} or odd number. Number of pixels of accuracy within
%           which to correct for chromatic shift.
%
%    timePoint: {1} or number. Timepoint at which to plot the kinetochore.
%
%    transpose: {0} or 1. Whether to tranpose the image. The program will
%           attempt to auto-correct for transposition if found, but cannot
%           auto-correct if the image has the same crop size in x and y.
%           If the crosses in the image do not align with the spots, try
%           setting transpose to 1. 
%
%    zoom: 0 or {1}. Whether or not to zoom into the specific kinetochore.
%           A value of 0 plots the whole cell, with just the coordinates
%           of the kinetochore plotted.
%
%    zProjRange: {2} or distance in pixels. Total distance about the
%           original spots over which to project in the z-coordinate.
%
%    KTborder: {1} or distance in um. Space either side of the KT. Here, 1
%           would result in an image with 1.1um space either side of the
%           extreme points of the KTs.
%
%    sepChannels: {0} or 1. Whether to do inset pictures of individual
%           channels rather than just one large picture.
% 
% 
% Copyright (c) 2024 C. C. Conway
% Code edited from dublShowSisterPair, copyright 2014 C. A. Smith and
% pairSpots, copyright 2018 C. A. Smith

if nargin<1
  error('Must supply a job.');
end
if all(job.options.jobProcess == 'zandt')
    error('Movie capability not yet available.')
    return
end

% set default options
opts.chanOrder =    [2 1 3]; % green, red, blue
opts.contrast =     {[0.1 0.9995] [0.1 0.9995] [0.1 0.9995]};
opts.newFig =       1;
opts.plotChannels = 1:3;
opts.plotTitles =   {'Mad2','CenpC','SKAP'};
opts.kinetochore =  1;
opts.subpixelate =  9;
opts.timePoint =    1;
opts.transpose =    0;
opts.zoom =         0;
opts.zProjRange =   2;
opts.KTborder =     1;
opts.sepChannels =  0;


% process options
opts = processOptions(opts, varargin{:});


%% Image and coordinate acquisition

metadataValidated = job.metadata.validated;

if metadataValidated
    [md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);
end

% get coordinate system and plot channels
coordSysChan = job.options.coordSystemChannel;
plotChans = sort(opts.plotChannels);
nChans = length(plotChans);
timePoint = opts.timePoint;

for iCh = plotChans
    opts.coordChans(iCh) = isfield(job.dataStruct{iCh},'initCoord');
end

opts.coordChans = find(opts.coordChans);
% get crop information, if any
crop = job.ROI.crop;
cropSize = job.ROI.cropSize;
if isempty(crop)
    cropSize = md.frameSize;
end
% specify RGB channel order [R G B]
chanOrder = opts.chanOrder;
    try
        % produce image file
        fullImg = zeros([cropSize([1 2]) 3]);
        for iChan = opts.plotChannels
            img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
            fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
            irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
            fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
        end
        isTransposed = 0; %use this to check later which case we are in
    catch ME %image may be transposed
        % produce image file
        fullImg = zeros([cropSize([2 1]) 3]);
        for iChan = opts.plotChannels
            img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
            fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
            irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
            fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
        end
        isTransposed = 1; %use this to check later
    end

    % produce figure environment
    figure(1);
    clf


    if length(opts.plotChannels) == 1
        imshow(fullImg(:,:,chanOrder(opts.plotChannels)),'InitialMagnification','fit');

    else
        imshow(fullImg,'InitialMagnification','fit'); %CCC test 2023.07.20
    hold on
    end
% calculate size of cropped region in pixels based on maxSisSep
pixelSize = job.metadata.pixelSize(1:3);

chrShift = job.options.chrShift.result;


%% GET IMAGE AND COORDINATE INFORMATION

% 

    if opts.transpose == 0
        isTransposed = any(size(kitReadImageStack(reader,md,1,2,crop,0)) ~= cropSize);
    else
        isTransposed = 1;
    end
    if isTransposed
        xyz_ind = [1 2 3];
    else
        xyz_ind = [2 1 3];
    end
    coordsPix = job.dataStruct{coordSysChan}.initCoord.allCoordPix(:,xyz_ind);


% get KT information
iCoord = opts.kinetochore;

% get origin pixel-coordinates
iCoordsPix = [];
for idx = iCoord
    wkCoordsPix = coordsPix(idx,:);
    iCoordsPix = cat(1, iCoordsPix, wkCoordsPix);
end

for iPos = 1:length(iCoord)
    scatter(iCoordsPix(iPos, 1), iCoordsPix(iPos, 2), 'xk', 'sizeData', 320, 'LineWidth', 3.5)
    scatter(iCoordsPix(iPos, 1), iCoordsPix(iPos, 2), 'xw', 'sizeData', 300, 'LineWidth', 2.5)
    text(iCoordsPix(iPos, 1), iCoordsPix(iPos, 2), string(iCoord(iPos)),'HorizontalAlignment','center','Color','w','FontSize',11);
end





    

    
    

    
%     text(cenKTx, cenKTy, string(iCoord),'HorizontalAlignment','center','Color','k','FontWeight','bold','FontSize',11);
%     text(cenKTx, cenKTy, string(iCoord),'HorizontalAlignment','center','Color','w','FontSize',11);


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

function [scatterx, scattery] = getScatterCoords(transposed, plotPix, coordRange, subpix)

    scatterx = [];
    scattery = [];
    if ~isempty(plotPix)
        if transposed
            scatterx = subpix*(plotPix(:,2) - coordRange(2,1)+0.5);
            scattery = subpix*(plotPix(:,1) - coordRange(1,1)+0.5);
        else
            scatterx = subpix*(plotPix(:,1) - coordRange(1,1)+0.5);
            scattery = subpix*(plotPix(:,2) - coordRange(2,1)+0.5);
        end
    end

end