function dublShowSisterPair(job,varargin)
% DUBLSHOWSISTERPAIR Plots image of dual-channel movie with coordinates
% for a given sister pair for a movie in a jobset created with KiT version
% 2.4.0 or later.
%
%    DUBLSHOWSISTERPAIR(JOB,...) Plots coordinates in two channels over
%    the movie image for a given sister pair at a given timepoint.
%
%    Options, defaults in {}:-
%
%    channelMap: {[2 1 3]} or some perturbation. Order in which the titled
%    channels are presented in the figure, where typically:
%    1=red, 2=green, 3=blue.
%
%    contrast: {{[0.1 1],[0.1 1],[0.1 1]}}, 'help', or similar. Upper and
%       lower contrast limits for each channel. Values must be in range
%       [0 1]. 'help' outputs the minimum and maximum intensity values as 
%       guidance, then requests values from the user.
%       Tips: - Increase the lower limit to remove background noise.
%             - Decrease the upper limit to increase brightness.
%
%    newFig: {0} or 1. Whether or not to show the sister pair in a new
%           figure.
%
%    plotChannels: {[1 2]} or some subset of [1 2 3]. Vector of channels
%           for plotting, where typically:
%               1=red, 2=green, 3=blue.
%
%    plotTitles: {'mNeonGreen','tagRFP','JF-646'} or similar. Titles of
%           each plot denoting which channel is which.
%
%    sisterPair: {1} or number. Sister pair within JOB being plotted.
%
%    subpixelate: {9} or odd number. Number of pixels of accuracy within
%           which to correct for chromatic shift.
%
%    timePoint: {1} or number. Timepoint at which to plot the sister pair.
%
%    transpose: {0} or 1. Whether to tranpose the image.
%
%    zoomScale: {1} or number. Magnification of a zoomed image as a
%           proportion of the default, so that 0.5 will zoom out, and 1.5
%           will zoom in.
%
%    zoom: 0 or {1}. Whether or not to zoom into the specific sister
%           pair. A value of 0 plots the whole cell, with just the
%           coordinates of the sister pair plotted.
%
%    zProject: 0, 1 or {-1}. Whether or not to project in the z-direction.
%           -1 will project in the 5 z-slices surrounding the sister pair.
%
% Copyright (c) 2014 C. A. Smith

if nargin<1
  error('Must supply a job.');
end

% set default options
opts.channelMap = [2 1 3]; % green, red, blue
opts.contrast = repmat({[0.1 1]},1,3);
opts.newFig = 0;
opts.plotChannels = 1:2;
opts.plotTitles = {'mNeonGreen','tagRFP','JF-646'};
opts.sisterPair = 1;
opts.subpixelate = 9;
opts.timePoint = 1;
opts.transpose = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.zProject = -1;

% process options
opts = processOptions(opts, varargin{:});
while length(opts.plotTitles) < length(opts.plotChannels)
    opts.plotTitles{length(opts.plotTitles)+1} = ''; %CCC edit 2023.06.09 - was stopping plots of 2 channels
end

%% Image and coordinate acquisition

% open movie
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);

% get coordinate system and plot channels
coordSysChan = job.options.coordSystemChannel;
plotChans = sort(opts.plotChannels);
nChans = length(plotChans);

% get sister information
sisPair = opts.sisterPair;
sisterList = job.dataStruct{coordSysChan}.sisterList;

% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% if only one channel, run the single channel version
if nChans == 1
    kitLog('Only one channel being shown. Running kitShowSisterPair instead.');
    kitShowSisterPair(job,'channel',plotChans,'contrast',opts.contrast,'newFig',opts.newFig,...
        'sisterPair',sisPair,'timePoint',timePoint,'title',opts.plotTitles{1},'transpose',opts.transpose,...
        'withinFig',0,'zoomScale',opts.zoomScale,'zoom',opts.zoom,'zProject',opts.zProject);
    return
end

% get pixel resolution
pixelSize = job.metadata.pixelSize;

% accumulate track information by channel and sister
trackCoord = nan(nChans,3,2);
for c = 1:nChans
    iChan = plotChans(c);
    for iSis = 1:2
        tk = trackIDs(iSis);
        track = job.dataStruct{iChan}.tracks(tk);
        
        startTime = track.seqOfEvents(1,1);
        if timePoint < startTime
            trackCoord(c,:,iSis) = nan(1,3);
        else
            trackCoord(c,:,iSis) = ...
                track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        end
    end
end

% calculate pair centre - check whether coordinate system is reference
if ismember(coordSysChan,plotChans)
    pos = find(plotChans==coordSysChan);
    centrePoint = nanmean(trackCoord(pos,:,:),3);
else
    % if not, then use first available channel
    coordSysChan = 1;
    centrePoint = nanmean(trackCoord(coordSysChan,:,:),3);
end
% convert to pixels
centrePoint = centrePoint./pixelSize;
centrePxl = round(centrePoint);

mapChans = opts.channelMap;
%if opts.transpose
    %dims = [2 1]; %CCC edit 2022/08/31 - ends up double transposing!
%else
    dims = [2 1];
%end

%% RGB image production

% predesignation of images
rgbImg = zeros([job.ROI.cropSize(dims), 3]);
rgbImgCS = zeros([job.ROI.cropSize(dims)*opts.subpixelate, 3]);
            
if opts.zoom
    % calculate spread about the centre pixels
    cropSpread = opts.zoomScale*(2./pixelSize);
	cropSpread = ceil(cropSpread);
    rgbCrpd= zeros([2*cropSpread(1)+1  2*cropSpread(2)+1  3]);
    rgbCrpdCS= zeros([opts.subpixelate*(2*cropSpread(1)+1) ...
                opts.subpixelate*(2*cropSpread(2)+1) ...
                3]);
end

% produce raw image
for c = 1:nChans
    iChan = plotChans(c);
    
    % read stack
    img = kitReadImageStack(reader,md,timePoint,iChan,job.ROI.crop,0);
    
    % max project over 5 z-slices around point
    if opts.zProject == 1
        img = max(img, [], 3);
    elseif opts.zProject == -1
        img = max(img(:,:,centrePxl(3)-2:centrePxl(3)+2), [], 3);
    else
        img = img(:,:,centrePxl(3));
    end
    if opts.transpose
        img = img';
    end
    rgbImg(:,:,mapChans(iChan)) = img;
    
end

% produce cropped image around pair centre
if opts.zoom
    
    % non-chromatic shifted regions
    xReg = [centrePxl(1)-cropSpread(1) centrePxl(1)+cropSpread(1)];
    yReg = [centrePxl(2)-cropSpread(2) centrePxl(2)+cropSpread(2)];
    
    for iChan = plotChans
        if opts.transpose
            rgbCrpd(:,:,mapChans(iChan)) = rgbImg(xReg(1):xReg(2),yReg(1):yReg(2),mapChans(iChan));
        else
            rgbCrpd(:,:,mapChans(iChan)) = rgbImg(yReg(1):yReg(2),xReg(1):xReg(2),mapChans(iChan));
        end
    end
    
    % produce chromatic shifted image
    first = 1;
    for iChan = setdiff(plotChans,coordSysChan) 
      if first

        [img1,img2] = ... 
            chrsComputeCorrectedImage(rgbCrpd(:,:,mapChans(coordSysChan)),rgbCrpd(:,:,mapChans(iChan)),job.options.chrShift.result{coordSysChan,iChan}, ...
            'subpixelate',opts.subpixelate);

        % give the subpixelated images to the image structure
        rgbCrpdCS(:,:,mapChans(coordSysChan)) = img1;
        rgbCrpdCS(:,:,mapChans(iChan)) = img2;

        first = 0;

      else
        [~,img2] = ...
            chrsComputeCorrectedImage(rgbCrpd(:,:,mapChans(coordSysChan)),rgbCrpd(:,:,mapChans(iChan)),job.options.chrShift.result{coordSysChan,iChan}, ...
            'subpixelate',opts.subpixelate);

        % give the new subpixelated images to the image structure
        rgbCrpdCS(:,:,mapChans(iChan)) = img2;

      end
    end
    
    % define contrast stretch for shifted cropped, and apply
    for iChan = plotChans
        if iscell(opts.contrast)
            irange = stretchlim(rgbCrpdCS(:,:,mapChans(iChan)),opts.contrast{iChan});
        elseif strcmp(opts.contrast,'help')
            % get maximum and minimum intensities of image
            intensityRange(1) = min(min(rgbCrpdCS(:,:,mapChans(iChan))));
            intensityRange(2) = max(max(rgbCrpdCS(:,:,mapChans(iChan))));
            % output possible range 
            fprintf('Channel %i. Range in third coordinate: [%i %i].\n',iChan,intensityRange(1),intensityRange(2));
            % request input from user
            userRange = input('Please provide range of intensities to image: ');
            while userRange(1)>userRange(2) 
                userRange = input('The maximum cannot be smaller than the minimum. Please try again: ');
            end
            irange = [userRange(1) userRange(2)];
            fprintf('\n');
        else
            irange = opts.contrast(iChan,:);
        end
        rgbCrpdCS(:,:,mapChans(iChan)) = imadjust(rgbCrpdCS(:,:,mapChans(iChan)), irange, []);
    end
    
else
    
    % produce chromatic shifted image
    first = 1;
    for iChan = setdiff(plotChans,coordSysChan) 
      if first
        
        [img1,img2] = ... 
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(iChan)),job.options.chrShift.result{coordSysChan,iChan}, ...
            'subpixelate',opts.subpixelate);
        
        % give the subpixelated images to the image structure
        rgbImgCS(:,:,mapChans(coordSysChan)) = img1;
        rgbImgCS(:,:,mapChans(iChan)) = img2;
        
        first = 0;
    
      else
        [~,img2] = ...
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(iChan)),job.options.chrShift.result{coordSysChan,iChan}, ...
            'subpixelate',opts.subpixelate);
        
        % give the new subpixelated images to the image structure
        rgbImgCS(:,:,mapChans(iChan)) = img2;
        
      end
    end
    
    % define contrast stretch for cropped images, and apply
    for iChan = plotChans
        if iscell(opts.contrast)
            irange = stretchlim(rgbImgCS(:,:,mapChans(iChan)),opts.contrast{iChan});
        elseif strcmp(opts.contrast,'help')
            % get maximum and minimum intensities of image
            intensityRange(1) = min(rgbImgCS(:,:,mapChans(iChan)));
            intensityRange(2) = max(rgbImgCS(:,:,mapChans(iChan)));
            % output possible range 
            fprintf('Channel %i. Range in third coordinate: [%i %i].\n',iChan,intensityRange(1),intensityRange(2));
            % request input from user
            userRange = input('Please provide range of intensities to image: ');
            while userRange(1)>userRange(2) 
                userRange = input('The maximum cannot be smaller than the minimum. Please try again: ');
            end
            irange = [userRange(1) userRange(2)];
            fprintf('\n');
        else
            irange = opts.contrast(iChan,:);
        end
        rgbImgCS(:,:,mapChans(iChan)) = imadjust(rgbImgCS(:,:,mapChans(iChan)), irange, []);
    end
    
end


%% Find coordinates to plot for each RGB image

% THIS LINE NOT REALLY NECESSARY
coord = trackCoord;

% adjust to pixels
for i = 1:2
    coord(:,:,i) = coord(:,:,i)./repmat(pixelSize,nChans,1);
end
coordCS = coord*opts.subpixelate;    

% correct coordinates for region position
if opts.zoom
    coord(:,1,:) = coord(:,1,:) - (xReg(1)+1);
    coord(:,2,:) = coord(:,2,:) - (yReg(1)-1);
    coordCS(:,1,:) = coordCS(:,1,:) - xReg(1)*opts.subpixelate + (opts.subpixelate+1)/2; % for subpix=9, +5; subpix=3, +2; subpix=5, +3 
    coordCS(:,2,:) = coordCS(:,2,:) - (yReg(1)*opts.subpixelate - (opts.subpixelate+1)/2); % CCC test edit from + to correlate with other y
end

%% Producing the figure

if nChans == 1
    plotImg = rgb2gray(rgbImg);
    if opts.zoom; plotImgCrpd = rgb2gray(rgbCrpd); end
    plotCoord = coord;
else
    plotImg = rgbImgCS;
    if opts.zoom; plotImgCrpd = rgbCrpdCS; end
    plotCoord = coordCS;
end 

% prepare figure environment
if opts.newFig
    figure
else
    figure(1)
end
clf

C = [ 0.9  0  0;
      0  0.9  0;
      0  0  0.9];
plotStyle = ['x' 'x' 'x'];

if opts.zoom
    tl = tiledlayout(nChans,nChans+1,'TileSpacing','tight','Padding','tight'); %CCC
    % plot individual channels
    bigImgInd = [];
    if nChans > 1
        for c = 1:nChans
            iChan = plotChans(c);
            nexttile(c*(nChans+1)) %CCC
            hold on
            imshow(plotImgCrpd(:,:,mapChans(iChan)))
            title(opts.plotTitles{c},'FontSize',16,'Color',C(mapChans(iChan),:))
            
            bigImgInd = [bigImgInd, (c-1)*(nChans+1)+1 : c*(nChans+1)-1 ];
        end
    end
    % plot dual image
    if isempty(bigImgInd)
        bigImgInd = 1;
    end
    nexttile(1, [nChans, nChans]) %CCC
    imshow(plotImgCrpd)
    title('All channels','FontSize',20)

    hold on

    XAxLim = get(gca,'XLim'); %CCC
    YAxLim = get(gca,'YLim'); %CCC
    if nChans>1 %CCC
        XAxPos = [(XAxLim(2)*0.95)-(1*opts.subpixelate/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
    else %CCC
        XAxPos = [(XAxLim(2)*0.95)-(1/pixelSize(1)) XAxLim(2)*0.95]; %CCC making 1um scalebar
    end %CCC
    YAxPos = [YAxLim(2)*0.95 YAxLim(2)*0.95]; %CCC
    line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3) %CCC making 1um scalebar

    % plot coordinates
    for c = 1:nChans
        iChan = plotChans(c);
        for i = 1:2
            nexttile(1, [nChans, nChans]) %CCC
            if opts.transpose
                plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                    'Color',C(mapChans(iChan),:),'Marker',plotStyle(mapChans(iChan)),'MarkerSize',30,'LineWidth',1) %CCC
            else
                plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                    'Color',C(mapChans(iChan),:),'Marker',plotStyle(mapChans(iChan)),'MarkerSize',30,'LineWidth',1) %CCC
            end
            if nChans > 1
                nexttile(c*(nChans+1)) %CCC
                if opts.transpose
                    plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                        'Color','k','Marker',plotStyle(mapChans(iChan)),'MarkerSize',15)
                else
                    plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                        'Color','k','Marker',plotStyle(mapChans(iChan)),'MarkerSize',15)
                end
            end
        end
        
            line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 2) %CCC adding 1um scale bar to mini images
    end
    
    %adding a scale bar trial
%     nexttile(1, [nChans, nChans])
%     XAxLim = get(gca,'XLim');
%     YAxLim = get(gca,'YLim');
%     if nChans>1
%         XAxPos = [(XAxLim(2)*0.95)-(1*opts.subpixelate/pixelSize(1)) XAxLim(2)*0.95];
%     else
%         XAxPos = [(XAxLim(2)*0.95)-(1/pixelSize(1)) XAxLim(2)*0.95];
%     end
%     YAxPos = [YAxLim(2)*0.95 YAxLim(2)*0.95];
%     line(XAxPos,YAxPos, 'Color', 'w', 'LineWidth', 3)

else
    % plot image
    imshow(plotImg)
    
    hold on
    % plot tracked channel's coordinates (points too close together on full image)
    for c = 1:nChans
        iChan = plotChans(c);
        for i = 1:2
            if opts.transpose
                plot(plotCoord(iChan,2,i),plotCoord(iChan,1,i),...
                    'Color',C(mapChans(iChan),:),'Marker',plotStyle(mapChans(iChan)),'MarkerSize',15)
            else
                plot(plotCoord(iChan,1,i),plotCoord(iChan,2,i),...
                    'Color',C(mapChans(iChan),:),'Marker',plotStyle(mapChans(iChan)),'MarkerSize',15)
            end
        end
    end
end

hold off

end
