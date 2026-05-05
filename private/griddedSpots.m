function imgDims = griddedSpots(job,varargin)
% IMGDIMS = GRIDDEDSPOTS Plots images of each spot in a grid for use by
% kitFilterSpots.
%
%    IMGDIMS = GRIDDEDSPOTS(JOB,...) Plots coordinates over images of each
%    spot localised in a given channel, outputting dimensions and positions
%    of each element of the grid.
%
%    Options, defaults in {}:-
%
%    channel: {1}, 2, 3 or 4. Which channel to show.
%
%    rawData: 0 or {1}. Whether or not to show raw data.
%
%    zoomed: 0 or {1}. Whether to show spot in centre of 8x8 µm or 1x1 µm
%            image.
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1; %default options
opts.rawData = 1;
opts.zoomed = 1;
opts = processOptions(opts,varargin{:});


% some variables
gridsep = 2;
if opts.zoomed
    imgHalfWidth = 0.5; %um
    opts.zProject = 0; %CCC addition for edited indivSpot 2023.10.18
else
    imgHalfWidth = 4; %um
    opts.zProject = 1; %CCC addition for edited indivSpot 2023.10.18
end
opts.bgcol = dot([0.94 0.94 0.94],[0.2989 0.5870 0.1140]); %background colour

% calculate imgWidth in pixels
pixelSize = job.metadata.pixelSize;
opts.imgHalfWidth = ceil(imgHalfWidth/pixelSize(1));
imgWidth = 2*opts.imgHalfWidth+1;

% suppress warnings
w = warning;
warning('off','all');

% get important information
dS = job.dataStruct{opts.channel};
nSpots = size(dS.initCoord.allCoord,1);
opts.imageSize = job.ROI.cropSize;
if isfield(job,'nROIs')
    nROIs = [' of ' num2str(job.nROIs)];
else
    nROIs = '';
end

% set up figure
figure; clf
fig_n=ceil(sqrt(nSpots));
fig_m=ceil(nSpots/fig_n);

pause(0.05);

% open movie and read stack
if length(job.ROI)>1
    [~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI(job.index).movie),'ROI');
    img = kitReadImageStack(reader,job.metadata,1,opts.channel,job.ROI(job.index).crop,0);
else
    [~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'ROI');
    img = kitReadImageStack(reader,job.metadata,1,opts.channel,job.ROI.crop,0);
end

% make empty gridded image
gridw = fig_m*imgWidth + (fig_m+1)*gridsep;
gridh = fig_n*imgWidth + (fig_n+1)*gridsep;
gridImg = opts.bgcol*ones(gridh,gridw);

% get all coordinate
if isfield(dS,'rawData') && opts.rawData
  iC = dS.rawData.initCoord;
else
  iC = dS.initCoord;
end
allCoords = iC.allCoordPix(:,1:3);

% form all spot positions
spotpos = 1:nSpots;
spotpos = spotpos(:);
spotpos = [mod(spotpos,fig_m) ceil(spotpos./fig_m)];
spotpos(spotpos(:,1)==0,1) = fig_m;

% preset allRnge array
rnge = [gridsep+(gridsep+imgWidth).*(spotpos(:,1)-1)+1 ...
        gridsep+(gridsep+imgWidth).*(spotpos(:,2)-1)+1];
% show each sister
for iSpot=1:nSpots
    % check whether any coordinates have been found, do nothing if so
    coords = allCoords(iSpot,:);
    if any(isnan(coords))
        allCoords(iSpot,:) = NaN;
        continue
    end
    
    % calculate the pixels in which to push new image
    [gridImg(rnge(iSpot,2):rnge(iSpot,2)+imgWidth-1, ...
        rnge(iSpot,1):rnge(iSpot,1)+imgWidth-1),coords] ...
        = indivSpot(img,coords,opts);
    
    % calculate position of coordinates to be plotted
    allCoords(iSpot,1:2) = coords(:,1:2)+rnge(iSpot,:);
    
end

% plot the full image and coordinates
pause(0.05);
imshow(gridImg,'Border','tight');
hold on
pause(0.05);
if opts.zoomed
    scatter(allCoords(:,1),allCoords(:,2),15*fig_m,'b','x')
else
    scatter(allCoords(:,1),allCoords(:,2),15*fig_m,'c','o')
end

markerSize = ceil(job.options.intensity.maskRadius / pixelSize(1));
for iCoord = 1:size(allCoords,1)
    drawCircle(allCoords(iCoord,1),allCoords(iCoord,2),markerSize,'w');
end

figtit = sprintf('Spot filtering: Image %i%s, channel %i',job.index,nROIs,opts.channel);
screenSz = get(0, 'MonitorPositions');
if size(screenSz,1) == 1
    set(gcf,'Resize','off','Name',figtit,'Units','characters',...
        'Position',[50 4.2 85 49],'NumberTitle','off'); %CCC change for Mac screen
elseif screenSz(2,1) < 0
    if opts.zoomed
        set(gcf,'Resize','off','Name',figtit,'Units','characters',...
            'Position',[-250 30 102 59],'NumberTitle','off'); %CCC change for work screen
    else
        set(gcf,'Resize','off','Name',figtit,'Units','characters',...
             'Position',[-290 2 155 87],'NumberTitle','off'); %CCC change for work screen
    end
else
    set(gcf,'Resize','off','Name',figtit,'Units','characters',...
        'Position',[250 0.5 115 64],'NumberTitle','off'); %CCC change for home screen
end
% elseif screenSz(2,1) < 0
%     set(gcf,'Resize','off','Name',figtit,'Units','characters',...
%         'Position',[-250 0 160 89],'NumberTitle','off'); %CCC change for work screen
% else
%     if max(spotpos(:,1)) < max(spotpos(:,2))
%         set(gcf,'Resize','off','Name',figtit,'Units','characters',...
%             'Position',[250 0 115 65],'NumberTitle','off'); %CCC change for home screen
%     else
%         set(gcf,'Resize','off','Name',figtit,'Units','characters',...
%             'Position',[250 0 125 65],'NumberTitle','off'); %CCC change for home screen
%     end
% end
%movegui(gcf,'center');

% save image dimensions and positions
imgDims = [rnge repmat(imgWidth,nSpots,1)];

% close the reader
close(reader);

% reset warnings
warning(w);

end

function drawCircle(x,y,r,color)
% Draws circle.

% Estimate pixels in circumference.
c = 2*pi*r;
theta = linspace(0,2*pi,ceil(c));
cx = x + r*cos(theta);
cy = y + r*sin(theta);
plot(cx, cy, [color '-'], 'Visible', 0);

end
