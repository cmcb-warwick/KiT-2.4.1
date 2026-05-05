function dublShowSingleKT(job,varargin)
% DUBLSHOWSINGLEKT Plots image of dual or triple-channel movie with coordinates 
% for a given kinetochore.
%
%    DUBLSHOWSINGLEKT(JOB,...) Plots coordinates in one, two or three channels
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
%    kinetochore: {1} or number. Kinetochore within JOB being plotted.
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
%           original spot over which to project in the z-coordinate.
%
%    KTborder: {1} or distance in um. Space either side of the KT. Here, 1
%           would result in an image of 2.2um x 2.2um (with 10% overage).
%
%    sepChannels: 0 or {1}. Whether to do inset pictures of individual
%           channels rather than just one large picture.
% 
%    maxZProject: {0} or 1. Whether to do a maximum intensity projection of
%           whole cell or not.
% 
% 
% Copyright (c) 2024 C. C. Conway
% Code edited from dublShowSisterPair, copyright 2014 C. A. Smith and
% pairSpots, copyright 2018 C. A. Smith

if nargin<1
  error('Must supply a job.');
end
if strcmp(job.options.jobProcess,'zandt')
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
opts.zoom =         1;
opts.zProjRange =   2;
opts.KTborder =     1;
opts.sepChannels =  1;
opts.maxZProject =  0;

%% TO DO
%make multiKTs plottable
%add KT numbers to plot

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

% calculate size of cropped region in pixels based on maxSisSep
pixelSize = job.metadata.pixelSize(1:3);
if opts.maxZProject
    zCropRange = cropSize(3);
else
    zCropRange = 1 + (opts.zProjRange*2);
end
if opts.zoom
    cropRange = (1.1*[opts.KTborder*2 opts.KTborder*2 zCropRange])./pixelSize;
    cropRange = round(cropRange);
end
chrShift = job.options.chrShift.result;





%% GET IMAGE AND COORDINATE INFORMATION

% 

if ~isfield(job.dataStruct{coordSysChan},'initCoord')
  kitLog('No initCoord present.');
  coordsPix = [];
% check whether or not this movie has an initCoord
elseif isfield(job.dataStruct{coordSysChan},'failed') && job.dataStruct{coordSysChan}.failed
  kitLog('Tracking failed to produce an initCoord.');
  coordsPix = [];
else
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
end

if length(opts.kinetochore)>1
    multiKTs = 1;
else
    multiKTs = 0;
end
% get KT information
iCoord = opts.kinetochore;

% get origin pixel-coordinates
iCoordsPix = coordsPix(iCoord,:);
centreCoords = round(iCoordsPix);
if multiKTs
    minX = min(centreCoords(:,1));
    maxX = max(centreCoords(:,1));
    minY = min(centreCoords(:,2));
    maxY = max(centreCoords(:,2));
    minZ = min(centreCoords(:,3));
    maxZ = max(centreCoords(:,3));
end
% IMAGE PROCESSING
% predesignate cropImg structure
if opts.zoom
    if multiKTs
        coordRange = [max(1,minX-cropRange(1)) min(maxX+cropRange(1),cropSize(1));...
                      max(1,minY-cropRange(2)) min(maxY+cropRange(2),cropSize(2));...
                      max(1,minZ-opts.zProjRange) min(maxZ+opts.zProjRange,cropSize(3))];
    else
        coordRange = [max(1,centreCoords(1)-cropRange(1)) min(centreCoords(1)+cropRange(1),cropSize(1));...
                      max(1,centreCoords(2)-cropRange(2)) min(centreCoords(2)+cropRange(2),cropSize(2));...
                      max(1,centreCoords(3)-opts.zProjRange) min(centreCoords(3)+opts.zProjRange,cropSize(3))];
    end
else
    coordRange = [1 cropSize(1); 1 cropSize(2); 1 cropSize(3)];
end

if isTransposed
    dims = [1 2];
else
    dims = [2 1];
end

subpix = opts.subpixelate;

rgbImg = zeros([job.ROI.cropSize(dims), 3]);
rgbImgCS = zeros([job.ROI.cropSize(dims)*subpix, 3]);

if ~isTransposed
      rgbCrpd = zeros([coordRange(2,2)-coordRange(2,1)+1,...
        coordRange(1,2)-coordRange(1,1)+1, 3]);

      rgbCrpdCS= zeros([subpix*(coordRange(2,2)-coordRange(2,1)+1) ...
            subpix*(coordRange(1,2)-coordRange(1,1)+1), 3]);
else
      rgbCrpd = zeros(coordRange(1,2)-coordRange(1,1)+1,...
        coordRange(2,2)-coordRange(2,1)+1, 3);

      rgbCrpdCS= zeros([subpix*(coordRange(1,2)-coordRange(1,1)+1) ...
            subpix*(coordRange(2,2)-coordRange(2,1)+1), 3]);
end
for c = plotChans

    % read stack
    img = kitReadImageStack(reader,md,timePoint,c,job.ROI.crop,0);
    
    % project over chosen number of z-slices around point
    if opts.maxZProject
        img = max(img, [], 3);
    elseif opts.zProjRange == 0
        img = img(:,:,centreCoords(3));
    else
        if multiKTs
            img = max(img(:,:,minZ-opts.zProjRange:maxZ+opts.zProjRange), [], 3);
        else
            img = max(img(:,:,centreCoords(3)-opts.zProjRange:centreCoords(3)+opts.zProjRange), [], 3);
        end
    end
    
    if isTransposed
        img = img';
    end
    
    rgbImg(:,:,chanOrder(c)) = img;

end

% produce cropped image around spot centre
    
% non-chromatic shifted regions
xReg = coordRange(1,:);
yReg = coordRange(2,:);

for c = plotChans
    if isTransposed
        rgbCrpd(:,:,chanOrder(c)) = rgbImg(xReg(1):xReg(2),yReg(1):yReg(2),chanOrder(c));
    else
        rgbCrpd(:,:,chanOrder(c)) = rgbImg(yReg(1):yReg(2),xReg(1):xReg(2),chanOrder(c));
    end
end

% produce chromatic shifted image
if nChans > 1
    first = 1;
    for c = setdiff(plotChans,coordSysChan) 
        if first
        [img1,img2] = ... 
            chrsComputeCorrectedImage(rgbCrpd(:,:,chanOrder(coordSysChan)),rgbCrpd(:,:,chanOrder(c)),chrShift{coordSysChan,c}, ...
            'subpixelate',subpix);
    
        % give the subpixelated images to the image structure
        rgbCrpdCS(:,:,chanOrder(coordSysChan)) = img1;
        rgbCrpdCS(:,:,chanOrder(c)) = img2;
    
        first = 0;
    
        else
        [~,img2] = ...
            chrsComputeCorrectedImage(rgbCrpd(:,:,chanOrder(coordSysChan)),rgbCrpd(:,:,chanOrder(c)),chrShift{coordSysChan,c}, ...
            'subpixelate',subpix);
    
        % give the new subpixelated images to the image structure
        rgbCrpdCS(:,:,chanOrder(c)) = img2;
    
        end
    end
else
    extraChan = setdiff(1:3, plotChans);
    [img1,~] = ... 
            chrsComputeCorrectedImage(rgbCrpd(:,:,chanOrder(plotChans(1))),rgbCrpd(:,:,chanOrder(extraChan(1))),chrShift{plotChans(1),extraChan(1)}, ...
            'subpixelate',subpix);
    rgbCrpdCS(:,:,chanOrder(plotChans(1))) = img1;
end
    
% define contrast stretch for shifted cropped, and apply
for c = plotChans
    if iscell(opts.contrast)
        irange = stretchlim(rgbCrpdCS(:,:,chanOrder(c)),opts.contrast{c});
    elseif strcmp(opts.contrast,'help')
        % get maximum and minimum intensities of image
        intensityRange(1) = min(min(rgbCrpdCS(:,:,chanOrder(c))));
        intensityRange(2) = max(max(rgbCrpdCS(:,:,chanOrder(c))));
        % output possible range 
        fprintf('Channel %i. Range in third coordinate: [%i %i].\n',c,intensityRange(1),intensityRange(2));
        % request input from user
        userRange = input('Please provide range of intensities to image: ');
        while userRange(1)>userRange(2) 
            userRange = input('The maximum cannot be smaller than the minimum. Please try again: ');
        end
        irange = [userRange(1) userRange(2)];
        fprintf('\n');
    else
        irange = opts.contrast(c,:);
    end
    rgbCrpdCS(:,:,chanOrder(c)) = imadjust(rgbCrpdCS(:,:,chanOrder(c)), irange, []);
end
    

    
    
% plot image in figure 1
if opts.newFig
    f = figure;
else
    f = figure(1);
end
clf
set(f,'Position',[200 20 1000 800]);
C = [0  1  0; 1  0  0; 0  0  1];

[cenKTx, cenKTy] = getScatterCoords(isTransposed, iCoordsPix, coordRange, subpix);


if nChans == 1
    imshow(im2gray(rgbCrpdCS(:,:,chanOrder(plotChans(1)))));
    hold on

    
    % origin coordinates in white
    scatter(cenKTx, cenKTy,'xk','sizeData',320,'LineWidth',3.5);
    scatter(cenKTx, cenKTy,'xw','sizeData',300,'LineWidth',2.5);
    if isTransposed
        xlim([0 size(rgbCrpdCS,2)])
        ylim([0 size(rgbCrpdCS,1)])
    else
        xlim([0 size(rgbCrpdCS,1)])
        ylim([0 size(rgbCrpdCS,2)])
    end
    XAxLim = get(gca,'XLim'); %CCC
    YAxLim = get(gca,'YLim'); %CCC
    XAxPos = [(XAxLim(2)*0.95)-(1/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
    YAxPos = [YAxLim(2)*0.95 YAxLim(2)*0.95]; %CCC
    line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3) %CCC making 1um scalebar
    title(sprintf('KT %d', iCoord),'FontSize',12);

elseif nChans > 1
    if opts.sepChannels
        tiledlayout(nChans,nChans+1,'TileSpacing','tight','Padding','tight')
        for c = 1:nChans
            nexttile(c*(nChans+1))
            hold on
            chanID = plotChans(c);
            title(opts.plotTitles{chanID},'FontSize',10,'Color',C(plotChans(c),:),'HorizontalAlignment','center','Margin',1);
            imshow(im2gray(rgbCrpdCS(:,:,chanOrder(chanID))));
    
            % origin coordinates in white
            scatter(cenKTx, cenKTy,'xk','sizeData',110,'LineWidth',2);
            scatter(cenKTx, cenKTy,'xw','sizeData',100,'LineWidth',1.5);
    
            if isTransposed
                xlim([0 size(rgbCrpdCS,2)])
                ylim([0 size(rgbCrpdCS,1)])
            else
                xlim([0 size(rgbCrpdCS,1)])
                ylim([0 size(rgbCrpdCS,2)])
            end
            XAxLim = get(gca,'XLim'); %CCC
            YAxLim = get(gca,'YLim'); %CCC
            XAxPos = [(XAxLim(2)*0.95)-(1*opts.subpixelate/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
            YAxPos = [YAxLim(2)*0.95 XAxLim(2)*0.95]; %CCC
            line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3) %CCC making 1um scalebar
        end
    
        nexttile(1,[nChans nChans]); %large multicolour image
        hold on
        title('All channels','FontSize',12,'HorizontalAlignment','center');
        imshow(rgbCrpdCS);
    
        
        % origin coordinates in white
        scatter(cenKTx, cenKTy,'xk','sizeData',320,'LineWidth',3.5);
        scatter(cenKTx, cenKTy,'xw','sizeData',300,'LineWidth',2.5);
        
        
        if isTransposed
            xlim([0 size(rgbCrpdCS,2)])
            ylim([0 size(rgbCrpdCS,1)])
        else
            xlim([0 size(rgbCrpdCS,1)])
            ylim([0 size(rgbCrpdCS,2)])
        end
        XAxLim = get(gca,'XLim'); %CCC
        YAxLim = get(gca,'YLim'); %CCC
        XAxPos = [(XAxLim(2)*0.95)-(1*opts.subpixelate/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
        YAxPos = [YAxLim(2)*0.95 YAxLim(2)*0.95]; %CCC
        line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3) %CCC making 1um scalebar
        sgtitle(sprintf('KT %d', iCoord),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    else
        title(sprintf('KT %d', iCoord),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
        imshow(rgbCrpdCS, 'InitialMagnification', 'fit');
    
        hold on

        % origin coordinates in white
        scatter(cenKTx, cenKTy,'xk','sizeData',320,'LineWidth',3.5);
        scatter(cenKTx, cenKTy,'xw','sizeData',300,'LineWidth',2.5);
        
        if isTransposed
            xlim([0 size(rgbCrpdCS,2)])
            ylim([0 size(rgbCrpdCS,1)])
        else
            xlim([0 size(rgbCrpdCS,1)])
            ylim([0 size(rgbCrpdCS,2)])
        end
        XAxLim = get(gca,'XLim'); %CCC
        YAxLim = get(gca,'YLim'); %CCC
        XAxPos = [(XAxLim(2)*0.95)-(1*opts.subpixelate/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
        YAxPos = [YAxLim(2)*0.95 YAxLim(2)*0.95]; %CCC
        line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3) %CCC making 1um scalebar
        sgtitle(sprintf('KT %d', iCoord),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    end
  
    hold off         
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