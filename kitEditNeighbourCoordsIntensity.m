function editedJob = kitEditNeighbourCoordsIntensity(job, varargin)
%this works only on individual jobs, not whole jobsets, so either have to
%load in and save individual jobs or get a mS and edit each cell in this
%
% This enables you to rotate a channel 180 degrees and re-track to find new
% initialCoords for that channel and new intensity measurements for all
% channels in the vicinity of the new spots, and also measure the intensity
% of the rotated channel in the vicinity of the other channel's spots.
%
% Largely created from kitTrackMovie.
% C. C. Conway, original 2022, edited to enable rotation 2025


% my head can't understand the case of rotating the coordinate system
% channel so please don't do that, rotate the others instead.
opts.coordChans = []; % the channel(s) that you want to re-track coordinates for
opts.intensityChans = []; %the channel(s) that you want to re-measure intensity for
opts.maskRadius = job.options.intensity.maskRadius;
opts.chanFlip180 = []; %vector of channels to rotate 180deg 20250828

% process user-defined options
opts = processOptions(opts,varargin{:});
if opts.maskRadius > 1
    opts.maskRadius = opts.maskRadius/1000;
end

[~, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
editedJob = job;
cChans = opts.coordChans;
iChans = opts.intensityChans;

if ~isempty(cChans)
    channels = setdiff(editedJob.analyzedChannels, editedJob.options.coordSystemChannel);
    chan = intersect(channels, cChans);
    for c = chan
        if isempty(intersect(opts.chanFlip180,c))
            flip180 = 0;
        else
            flip180 = 1;
        end
        % Find particle coordinates 
        kitLog('Finding particle coordinates in channel %d',c);
        editedJob = kitFindCoords_edit(editedJob, reader, c, flip180);
    
        % Transform coordinates into plane.
        planeChan = editedJob.options.coordSystemChannel;
        kitLog('Transforming coordinates in channel %d to plane from channel %d',c,planeChan);
        editedJob.dataStruct{c}.planeFit = editedJob.dataStruct{planeChan}.planeFit;
        editedJob.dataStruct{c} = kitFitPlane(editedJob,reader,editedJob.dataStruct{c},c,1); %this doesn't need changing I think

        if isfield(editedJob.dataStruct{editedJob.options.coordSystemChannel}, 'tracks')
        % Update coords obtained from pairing
            kitLog('Updating pairing coordinates in channel %d',c);
            editedJob = updateTracks(editedJob, c);
        end
        
        if ~ismember(c, iChans)
            % Measure intensity
            kitLog('Measure particle intensity in channel %d',c);
            editedJob.options.intensity.maskRadius = opts.maskRadius;
            editedJob = kitLocalIntensity_edit(editedJob, reader, editedJob.metadata, c, editedJob.options.intensity,opts.chanFlip180);
        end
    end
end

if ~isempty(iChans)
    for c = iChans
        % Measure intensity
        kitLog('Measure particle intensity in channel %d',c);
        editedJob.options.intensity.maskRadius = opts.maskRadius;
        editedJob = kitLocalIntensity_edit(editedJob, reader, editedJob.metadata, c, editedJob.options.intensity,opts.chanFlip180);
    end
end






end
%%
function job = updateTracks(job, iChan)

% Get sisterIdxArray as in pairSpots
nSis = size(job.dataStruct{1,2}.sisterList,2);

sisterIdxArray = nan(nSis,2);

for i = 1:nSis
    loc1 = i;
    loc2 = i + nSis;
    idx1 = job.dataStruct{1,2}.trackList(loc1).featIndx;
    idx2 = job.dataStruct{1,2}.trackList(loc2).featIndx;
    sisterIdxArray(i,1) = idx1;
    sisterIdxArray(i,2) = idx2;
end

% set up empty sisterLists
emptySisterList = struct('trackPairs',[],...
    'coords1',[],...
    'coords2',[],...
    'sisterVectors',[],...
    'distances',[]);
emptyTracks = struct('tracksFeatIndxCG',0,...
    'tracksCoordAmpCG',[],...
    'seqOfEvents', [1 1 1; 1 2 1],...
    'coordAmp4Tracking',[]);



% construct new sisterList
sisterList = emptySisterList;
sisterList(1).trackPairs(:,1) = 1:size(sisterIdxArray,1);
sisterList(1).trackPairs(:,2) = sisterList(1).trackPairs(:,1)+size(sisterIdxArray,1);

% get microscope coordinates and standard deviations in µm
if isfield(job.dataStruct{iChan},'initCoord')
    coords = job.dataStruct{iChan}.initCoord(1).allCoord;
    amps = job.dataStruct{iChan}.initCoord(1).amp;
else
    coords = [];
    amps = [];
end
% get plane coordinates if applicable
if isfield(job.dataStruct{iChan},'planeFit') && ~isempty(job.dataStruct{iChan}.planeFit(1).plane)
    rotCoords = job.dataStruct{iChan}.planeFit(1).rotatedCoord;
else
    rotCoords = coords;
end

if ~isempty(rotCoords)
    for iSis = 1:size(sisterIdxArray,1)
        sisterList(iSis).coords1 = rotCoords(sisterIdxArray(iSis,1),:);
        sisterList(iSis).coords2 = rotCoords(sisterIdxArray(iSis,2),:);
        sisterList(iSis).distances(:,1) = sqrt(sum((sisterList(iSis).coords1(:,1:3) - sisterList(iSis).coords2(:,1:3)).^2,2));
        sisterList(iSis).distances(:,2) = sqrt(sum((sisterList(iSis).coords1(:,4:6) - sisterList(iSis).coords2(:,4:6)).^2,2));
    end
else
    for iSis = 1:size(sisterIdxArray,1)
        sisterList(iSis).coords1 = NaN(1,6);
        sisterList(iSis).coords2 = NaN(1,6);
        sisterList(iSis).distances = NaN(1,2);
    end
end

if isempty(sisterIdxArray)
    % assign in the case no sisters found
    tracks = emptyTracks;
end

% construct new tracks
if ~isempty(rotCoords)
    for iTrack = 1:length(sisterIdxArray(:))

        % ensure tracks are re-allocated for new cell
        tracks(iTrack) = emptyTracks;
        tracks(iTrack).tracksFeatIndxCG = sisterIdxArray(iTrack);
        tracks(iTrack).tracksCoordAmpCG(:,1:3) = coords(sisterIdxArray(iTrack),1:3);
        tracks(iTrack).tracksCoordAmpCG(:,[4 8]) = amps(sisterIdxArray(iTrack),1:2);
        tracks(iTrack).tracksCoordAmpCG(:,5:7) = coords(sisterIdxArray(iTrack),4:6);
        tracks(iTrack).coordAmp4Tracking(:,1:3) = rotCoords(sisterIdxArray(iTrack),1:3);
        tracks(iTrack).coordAmp4Tracking(:,[4 8]) = amps(sisterIdxArray(iTrack),1:2);
        tracks(iTrack).coordAmp4Tracking(:,5:7) = rotCoords(sisterIdxArray(iTrack),4:6);

    end
else
    for iTrack = 1:length(sisterIdxArray(:))

        % ensure tracks are re-allocated for new cell
        tracks(iTrack) = emptyTracks;
        tracks(iTrack).tracksFeatIndxCG = sisterIdxArray(iTrack);
        tracks(iTrack).tracksCoordAmpCG = NaN(1,8);
        tracks(iTrack).coordAmp4Tracking = NaN(1,8);

    end
end

job.dataStruct{iChan}.tracks = tracks;
job.dataStruct{iChan}.sisterList = sisterList;
job = kitExtractTracks(job,iChan);

end

%%
function job=kitFindCoords_edit(job, reader, channel, flip180)
%KITFINDCOORD Find kinetochore coordinates
%
% SYNOPSIS: job=kitFindCoords(job, raw, channel)
%
% INPUT job: Struct containing tracking job setup options.
%            Requires at least the following fields:
%
%       reader: BioFormats reader.
%
%       channel: Channel to find coords in.
%
%       flip180: whether to rotate channel 180deg or not CCC 20250828 for
%       MPS1 analysis
%
% OUTPUT job: as input but with updated values.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2016 C. A. Smith
% Modified by C. C. Conway (2022 2025)

% Method of fixing spot locations: centroid or Gaussian MMF, or none.
method = job.options.coordMode{channel};
% Handle old jobset versions.
if ~isfield(job.options,'spotMode')
  spotMode = 'histcut';
else
  % Method of identify first spot candidates: Histogram cut 'histcut',
  % adaptive threshold 'adaptive', or  multiscale wavelet product 'wavelet'.
  spotMode = job.options.spotMode{channel};
end

% Set up data struct.
options = job.options;

nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
filters = createFilters(ndims,job.dataStruct{channel}.dataProperties);

% Read image
movie = kitReadWholeMovie_edit(reader,job.metadata,channel,job.ROI.crop,0,1,flip180); %CCC 20250828
[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% Initialize output structure
localMaxima = repmat(struct('cands',[]),nFrames,1);

% Find candidate spots.
switch spotMode
  case 'histcut'
    kitLog('Detecting particle candidates using unimodal histogram threshold');
    spots = cell(nFrames,1);
    for i=1:nFrames
      img = movie(:,:,:,i);
      spots{i} = histcutSpots(img,options,job.dataStruct{channel}.dataProperties);
    end

  case 'adaptive'
    kitLog('Detecting particle candidates using adaptive thresholding');
    if ~isfield(options,'adaptiveLambda')
      options.adaptiveLambda = 1; % if lambda not set, then set to 1 as default
    end
    if ~isfield(options,'realisticNumSpots')
      options.realisticNumSpots = 92;
    end
    if ~isfield(options, 'globalBackground')
      options.globalBackground = 1;
    end
    if ~isfield(options, 'robustObjectiveFn')
      options.robustObjectiveFn = 1;
    end
    spots = adaptiveSpots(movie,options.adaptiveLambda,...
                          options.realisticNumSpots,options.globalBackground,...
			  options.robustObjectiveFn,...
                          options.debug.showAdaptive);
    
  case 'CFAR'
    kitLog('Detecting particle candidates using CFAR detection');
    if ~isfield(options,'realisticNumSpots')
      options.realisticNumSpots = 92;
    end
    spots = spotsDetectCFAR(movie, options.realisticNumSpots,...
                            options.debug.showAdaptive);

  case 'wavelet'
    kitLog('Detecting particle candidates using multiscale wavelet product');
     options.waveletLevelThresh = waveletAdapt(movie,options);
     job.options = options;
     kitSaveJob(job); % Record used value.

    for i=1:nFrames
      img = movie(:,:,:,i);
      spots{i} = waveletSpots(img,options);
    end
    
  case 'neighbour'
    kitLog('Detecting particle candidates using neighbouring channel');
    refDataStruct = job.dataStruct{options.coordSystemChannel};
    [spots,spotIDs] = neighbourSpots_edit(movie,refDataStruct,channel,job.metadata,options, flip180); %CCC edit from neighbourSpots 2022/10/19
    
  case 'manual'
    kitLog('Detecting particle candidates using manual detection');
    spots = manualDetection(movie,job.metadata,options);
    
  otherwise
    error('Unknown particle detector: %s',spotMode);
end

nSpots = zeros(nFrames,1);
for i=1:nFrames
  nSpots(i) = size(spots{i},1);

  % Round spots to nearest pixel and limit to image bounds.
  if nSpots(i) > 1
    spots{i} = bsxfun(@min,bsxfun(@max,round(spots{i}),1),[imageSizeX,imageSizeY,imageSizeZ]);
  end

  % Store the cands of the current image
  % TODO this is computed in both spot detectors, just return it.
  img = movie(:,:,:,i);
  if verLessThan('images','9.2')
    background = fastGauss3D(img,filters.backgroundP(1:3),filters.backgroundP(4:6));
  else
    background = imgaussfilt3(img,filters.backgroundP(1:3),'FilterSize',filters.backgroundP(4:6));
  end
  localMaxima(i).cands = spots{i};
  if strcmp(spotMode,'neighbour')
      localMaxima(i).spotID = spotIDs{i};
  end
  if nSpots(i) > 1
    spots1D = sub2ind(size(img),spots{i}(:,1),spots{i}(:,2),spots{i}(:,3));
  else
    spots1D = [];
  end
  localMaxima(i).candsAmp = img(spots1D);
  localMaxima(i).candsBg = background(spots1D);

  % Visualize candidates.
  if options.debug.showMmfCands ~= 0
    showSpots(img,spots{i});
    title(['Local maxima cands n=' num2str(size(spots{i},1))]);
    drawnow;
    switch options.debug.showMmfCands
      case -1
        pause;
      case -2
        keyboard;
    end
  end
end
kitLog('Average particles per frame: %.1f +/- %.1f',mean(nSpots),std(nSpots));

% Refine spot candidates.
switch method
  case 'centroid'
    job = kitCentroid(job,movie,localMaxima,channel);
  case 'gaussian'
    job = kitMixtureModel(job,movie,localMaxima,channel);
    case 'none' % Edit from 'norefine' to 'none' A.D. 07/02/2023
    % No refinement. Copy localMaxima to initCoords.
    initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0, ...
                                  'amp',[],'bg',[]);
    initCoord(1).localMaxima = localMaxima;
    for i=1:nFrames
      initCoord(i).nSpots = size(localMaxima(i).cands,1);
      initCoord(i).allCoordPix = [localMaxima(i).cands(:,[2 1 3]) ...
                          0.25*ones(initCoord(i).nSpots,3)];
      initCoord(i).allCoord = bsxfun(@times, initCoord(i).allCoordPix,...
        repmat(job.metadata.pixelSize,[1 2]));
      initCoord(i).amp = [localMaxima(i).candsAmp zeros(initCoord(i).nSpots,1)];
    end
    % Store data.
    job.dataStruct{channel}.initCoord = initCoord;
    job.dataStruct{channel}.failed = 0;
  otherwise
    error(['Unknown coordinate finding mode: ' job.coordMode]);
end
end
%%
function [spots,spotIDs] = neighbourSpots_edit(movie,refDataStruct,channel,metadata,options, flip180)
% NEIGHBOURSPOTS Finds spots in 3D in neighbouring channels using Gaussian
% mixture model fitting
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2017 C. A. Smith
% Edited 2022 C. C. Conway

% Input + initialization

% get all required coordinate information
refInitCoord = refDataStruct.initCoord;

% get pixel and chromatic shift information
pixelSize = metadata.pixelSize(1:3);
chrShift = options.chrShift.result{options.coordSystemChannel,channel}(:,1:3);
chrShift = chrShift./pixelSize; chrShift(1:2) = chrShift([2 1]);
% 20250828
if flip180
    chrShift(1) = chrShift(1)*(-1);
    chrShift(2) = chrShift(2)*(-1);
end
% image size
imageSize = size(movie);
nFrames = metadata.nFrames;
nPlanes = imageSize(3);
% get channel orientation (1 if channel outside reference, 0 otherwise)
refChanPos = find(options.neighbourSpots.channelOrientation==options.coordSystemChannel);
chanPos = find(options.neighbourSpots.channelOrientation==channel);
chanOrient = (chanPos>refChanPos);

% get list of frames over which to look for neighbours
if isempty(options.neighbourSpots.timePoints{channel})
    timePoints = 1:nFrames;
else
    timePoints = options.neighbourSpots.timePoints{channel};
end
if isfield(refDataStruct.initCoord(1),'exceptions')
    emptyFrames = refDataStruct.initCoord(1).exceptions.emptyFrames;
else
    emptyFrames = [];
end
emptyFrames = union(emptyFrames, setxor(timePoints,1:nFrames));
timePoints = setxor(emptyFrames,1:nFrames); % final list of frames
timePoints = timePoints(:)'; % ensure is a row vector - setxor has undergone a change over MATLAB versions

%get number of z-slices
if isempty(options.neighbourSpots.zSlices{channel})
    zSlices = 1:nPlanes;
else
    zSlices = options.neighbourSpots.zSlices{channel};
end

%turn warnings off
warningState = warning();

% Get frame with plane fits

if ~strcmp(options.jobProcess,'chrshift')
  % get plane fit information from reference channel
  planeFit = refDataStruct.planeFit;
 
  % find frames with and without plane fit
  framesNoPlane = [];
  for iFrame = 1:nFrames
    if isempty(planeFit(iFrame).plane)
      framesNoPlane = [framesNoPlane iFrame];
    end
  end
else
  framesNoPlane = 1;
end

% Create Mask

r = round(options.neighbourSpots.maskRadius/pixelSize(1));
if ~strcmp(options.neighbourSpots.maskShape, 'cone')
  se = strel('disk', r, 0);
  asymMask = 0;
  mask = double(se.getnhood());
  mask(mask == 0) = nan;
  if strcmp(options.neighbourSpots.maskShape, 'semicirc')
    mask(r+2:end,:) = nan;
    asymMask = 1;
  end
else
  mask = ones(r);
  % Make conical mask.
  ycut = tand(options.neighbourSpots.maskConeAngle)*(0:r-1);
  for i=1:r
    ycutFlr = floor(ycut(i));
    if ycutFlr>=1
      mask(1:ycutFlr,i) = nan;
    end
  end
  % Add asymmetry and orient.
  mask = [flipud(mask); nan(r-1,r)];
  % Mirror
  mask = [fliplr(mask) mask(:,2:r)];
  asymMask = 1;
end

if asymMask == 1
  lMask = mask; % For -ve x pole attached.
  rMask = flipud(mask); % For +ve x pole attached.
else 
  lMask = mask;
  rMask = mask;
end
if ~isempty(framesNoPlane)
    se = strel('disk', r, 0);
    noPlaneMask = double(se.getnhood());
    noPlaneMask(noPlaneMask == 0) = nan;
end

maskWarning=0;

% Local maxima detection relative to reference channel, based on job process

switch options.jobProcess
    
  case 'zandt'
      
    % get trackList and sisterList from reference channel
    refTrackList = refDataStruct.trackList;
    refSisterList = refDataStruct.sisterList;

    % get number of sisters
    if ~isempty(refSisterList(1).trackPairs)
        nSisters = length(refSisterList);
    else
        warning('No sisterList found for channel %d. Tracking failed.',options.coordSystemChannel)
        spots = repmat({[]},nFrames,1);
        spotIDs = repmat({[]},nFrames,1);
        return
    end
    
    % create structure to store sister number and track numbers numbers
    referenceIDs = nan(nSisters,3);
    % create spots and spotIDs structures
    spots = repmat({[]},nFrames,1);
    spotIDs = repmat({[]},nFrames,1);

    for iSisPair = 1:nSisters

        % get sisterList, and its track and spot IDs
        sisList = refSisterList(iSisPair);
        iTracks = refSisterList(1).trackPairs(iSisPair,1:2);
        spotIDbySister{iSisPair}  = [refTrackList(iTracks(1)).featIndx refTrackList(iTracks(2)).featIndx];
        referenceIDs(iSisPair,:) = [iSisPair iTracks];

        for iFrame = timePoints
          % get initCoord for this timepoint
          initCoord = refInitCoord(iFrame);

          for iSis = 1:2

              trackID = iTracks(iSis);
              spotID = spotIDbySister{iSisPair}(iFrame,iSis);
              % check that this coordinate has a spotID
              if isnan(spotID)  
                continue
              end

              % get coordinate information for this each sister, then transpose
              % (the results of MMF are transposed, so need to be transposed back)
              coords = initCoord.allCoordPix(spotID,[2 1 3]);
              
              % chromatic shift coordinates to new channel's frame, then
              % find nearest whole pixel
              coords = coords + chrShift;
              coords = round(coords);
              
              % define and check the image range
              range = [coords(1)-r coords(1)+r;...
                       coords(2)-r coords(2)+r;...
                       coords(3)   coords(3) ];
              % CCC edit 2022/10/19 from coords(2)-r coords(2)-r to coords(2)-r coords(2)+r

              if any(isnan(range(:))) || any(range(:,1)<=0) || any(range(:,2)>imageSize(1:3)')
                  continue
              elseif range(3,1)<min(zSlices) || range(3,2)>max(zSlices)
                  continue
              end
                  
              % realign mask into correct orientation for this track
              if sum(framesNoPlane == iFrame) > 0
                  mask = noPlaneMask;
              else
                  if refTrackList(trackID).attach > 0
                    if chanOrient
                      mask = rMask;
                    else
                      mask = lMask;
                    end
                  elseif ~chanOrient
                      mask = rMask;
                  else
                      mask = lMask;
                  end
                  
                  if ~isempty(planeFit(iFrame).planeVectors)
                    % rotate mask into coordinate system
                    angle = acos(planeFit(iFrame).planeVectors(1,1))*180/pi;
                    mask = imrotate(mask, angle, 'nearest', 'crop');
                  else
                    if ~maskWarning
                      warning('No rotation vector for mask.');
                      maskWarning=1;
                    end
                  end
              end
                  
              % get frame
              image = movie(:,:,:,iFrame);
              
              % get intensity values of mask-derived spot vicinity
              imageMask = mask .* image(range(1,1):range(1,2), ...
                  range(2,1):range(2,2), range(3,1):range(3,2));

              % find local maximum intensity
              locMax1DIndx = find(imageMask==nanmax(imageMask(:)));
              if isempty(locMax1DIndx)
                  continue
              elseif length(locMax1DIndx)>1
                  locMax1DIndx = locMax1DIndx(1); % IDEALLY SHOULD BE THE SPOT CLOSEST TO THE ORIGINAL SPOT
              end

              [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx);
              % correct coordinates to full image
              locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
              locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
              locMaxCrd(3) = coords(3);

              if any(locMaxCrd>imageSize(1:3))
                  continue
              end

              % compile coordinates
              spots{iFrame} = [spots{iFrame}; locMaxCrd];
              spotIDs{iFrame} = [spotIDs{iFrame}; spotID];

          end %iSis
        end %iFrame

    end %iSisPair
    
  case {'zonly','chrshift'}
    
    % get number of spots from reference initCoord
    %20250828
    nSpots = length(refInitCoord.allCoord);
    if nSpots == 1
        nSpots = refInitCoord.nSpots;
    end
    % create spots and spotIDs structures
    spots = {[]};
    spotIDs = {[]};

    for iSpot = 1:nSpots

      % get coordinate information
      coords = refInitCoord(1).allCoordPix(iSpot,[2 1 3]);
      
      % chromatic shift coordinates to new channel's frame, then
      % find nearest whole pixel
      coords = coords + chrShift;
      coords = round(coords);

      % define and check the image range
      range = [coords(1)-r coords(1)+r;...
               coords(2)-r coords(2)+r;...
               coords(3)   coords(3) ];
     % CCC edit 2022/10/19 from coords(2)-r coords(2)-r to coords(2)-r coords(2)+r
           
      if any(isnan(range(:))) || any(range(:,1)<=0) || any(range(:,2)>imageSize(1:3)')
          continue
      end
  
      % get frame
      image = movie(:,:,:,1);
        
      % get intensity values of mask-derived spot vicinity
      imageMask = mask .* image(range(1,1):range(1,2), ...
          range(2,1):range(2,2), range(3,1):range(3,2));

      % find local maximum intensity
      locMax1DIndx = find(imageMask==nanmax(imageMask(:)));
      if isempty(locMax1DIndx)
          continue
      elseif length(locMax1DIndx)>1
          locMax1DIndx = locMax1DIndx(1); % IDEALLY SHOULD BE THE SPOT CLOSEST TO THE ORIGINAL SPOT
      end
        
      [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx);
      % correct coordinates to full image
      locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
      locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
      locMaxCrd(3) = coords(3);
      
      if any(locMaxCrd>imageSize(1:3))
          continue
      end

      % compile coordinates for MMF
      spots{1} = [spots{1}; locMaxCrd];
      spotIDs{1} = [spotIDs{1}; iSpot];

    end %iSpot
    
end

%go back to original warnings state
warning(warningState);
end
%%
function job=kitLocalIntensity_edit(job, reader, metadata, channel, opts, chanFlip180)
% KITLOCALINTENSITY_EDIT Measure intensity local to kinetochores with
% option to rotate coords
%
% Copyright (c) 2018 C A Smith
% Modified by C. C. Conway (2022 & 2025)

% Preparation: get metadata, produce structures

% get some metadata
nFrames = metadata.nFrames;

% check if this channel contains coordinate or track information, or neither
if length(job.dataStruct)<channel || ~isfield(job.dataStruct{channel},'initCoord')
  refChan = job.options.coordSystemChannel;
else
  refChan = channel;
end
%chrShift = job.options.chrShift.result{job.options.coordSystemChannel,channel};
%chrShift = chrShift(1:3)./metadata.pixelSize(1:3);
% 2022.10.06 CCC removed
initCoord = job.dataStruct{refChan}.initCoord;
planeFit = job.dataStruct{refChan}.planeFit;
if nFrames==1 || ~isfield(job.dataStruct{refChan},'trackList')
  useTracks = 0;
  nKTs = size(initCoord(1).allCoord,1);
else
  useTracks = 1;
  trackList = job.dataStruct{refChan}.trackList;
  nKTs = length(trackList);
end
nChans = 4;

% predesignate intensity structure
intStruct(1:nFrames) = struct(...
    'intensity',nan(nKTs,nChans),...
    'intensity_median',nan(nKTs,nChans),...
    'intensity_min',nan(nKTs,nChans),...
    'intensity_max',nan(nKTs,nChans),...
    'intensity_ratio',nan(nKTs,nChans),...
    'maskCoord',nan(nKTs,3),...
    'maxCoord',nan(nKTs,3),...
    'angleToMax',nan(nKTs,nChans),...
    'distToMax',nan(nKTs,nChans));

% convert pole shift to pixels
poleShiftPixels = round(opts.poleShift / metadata.pixelSize(1));

% Gaussian filter
hgauss = fspecial('gauss');


% Create mask
r = round(opts.maskRadius / metadata.pixelSize(1));
if ~strcmp(opts.maskShape, 'cone')
  se = strel('disk', r, 0);
  asymMask = 0;
  mask = double(se.getnhood());
  mask(mask == 0) = nan;
  if strcmp(opts.maskShape, 'semicirc')
    mask(r+2:end,:) = nan;
    asymMask = 1;
  end
else
  mask = ones(r);
  % Make conical mask.
  ycut = tand(opts.maskConeAngle)*(0:r-1);
  for i=1:r
    ycutFlr = floor(ycut(i));
    if ycutFlr>=1
      mask(1:ycutFlr,i) = nan;
    end
  end
  % Add asymmetry and orient.
  mask = [flipud(mask); nan(r-1,r)];
  % Mirror
  mask = [fliplr(mask) mask(:,2:r)];
  asymMask = 1;
end

if asymMask == 1
  lMask = mask; % For -ve x pole attached.
  rMask = flipud(mask); % For +ve x pole attached.
end

maskWarning=0;


% Get intensities

prog = kitProgress(0);
chans = find(opts.execute);
nChans = length(chans);
for iChan = chans
    % get chromatic shift - CCC addition 2022.10.06
    % chromatic shift from refChan to coordinate channel
    chrShift_coordSys = job.options.chrShift.result{refChan, job.options.coordSystemChannel};
    %don't rotate the coordSysChan!! That's not the point of this!!! I'm not checking for it
    %here so if you have and the data doesn't work then that's your problem
    if ~isempty(intersect(chanFlip180,refChan))
        chrShift_coordSys(1) = chrShift_coordSys(1)*(-1);
        chrShift_coordSys(2) = chrShift_coordSys(2)*(-1);
    end
    % chromatic shift from coordinate channel to channel under consideration
    chrShift_iChan = job.options.chrShift.result{job.options.coordSystemChannel, iChan};
    if ~isempty(intersect(chanFlip180,iChan))
        flip180 = 1;
        chrShift_iChan(1) = chrShift_iChan(1)*(-1);
        chrShift_iChan(2) = chrShift_iChan(2)*(-1);
    else
        flip180 = 0;
    end 
    chrShift = chrShift_coordSys + chrShift_iChan;
    chrShift = chrShift(1:3)./metadata.pixelSize(1:3);
    
    % read whole movie
    if length(job.ROI)>1
      movie = kitReadWholeMovie_edit(reader, metadata, iChan, job.ROI(job.index).crop, 0, 0, flip180); %20250828
    else
      movie = kitReadWholeMovie_edit(reader, metadata, iChan, job.ROI.crop, 0, 0, flip180); %20250828
    end

    for t=1:nFrames

        stack = movie(:,:,:,t);

        if opts.gaussFilterSpots
          % Gaussian filter each z-plane.
          for z=1:size(stack,3)
            stack(:,:,z) = imfilter(stack(:,:,z), hgauss, 'replicate');
          end
        end

        if t==1 && iChan==channel
          % Take overall fluorescence estimate.
          back = mean(stack(:));

          % Estimate background and signal.
          [backMode,sigVal,imgHist] = estBackgroundAndSignal(stack);
          intensityDistF(:) = imgHist(:,2);
          intensityDistX(:) = imgHist(:,1);

          % Estimate background as difference between signal mode and background mode.
          backDiff = sigVal - backMode;
        end

        for j=1:nKTs
            
          if useTracks
            % Map track back to coords. FIXME use trackList instead
            pixIdx = trackList(j).featIndx(t);
          else
            pixIdx = j;
          end

          if ismember(pixIdx,1:nKTs)
            % Read off intensity in radius r around spot.
            pixCoords = initCoord(t).allCoordPix(pixIdx,[2 1 3]);

            if sum(isnan(pixCoords)) == 0

              if useTracks && trackList(j).attach ~= 0 && ...
                    size(pixCoords,2)==size(planeFit(t).planeVectors,2) && ...
                    refChan ~= channel
                % Shift mask toward pole if measuring intensity of
                % non-localised channel, if known sister.
                maskCoords = pixCoords + [sign(trackList(j).attach),0,0] * ...
                    poleShiftPixels*planeFit(t).planeVectors;
              else
                maskCoords = pixCoords;
              end
              % chromatic shift from coordSysChan to the channel being measured
              maskCoords = maskCoords + chrShift([2 1 3]);

              % Extract intensity.
              x = max(1,round(maskCoords(1)));
              y = max(1,round(maskCoords(2)));
              z = max(1,round(maskCoords(3)));
              if z > size(stack,3)
                warning('Spot outside frame boundaries');
                continue
              end
              if size(stack,3)>1
                imgPlane = stack(:,:,z);
              else
                imgPlane = stack;
              end
              [mx,my] = size(imgPlane);
              % If too close to edge, skip.
              if x <= r || y <= r || x >= mx-r || y >= my-r
                continue
              end

              % If asymmetric mask, choose left/right mask based on pole
              % attachment.
              if useTracks && asymMask
                if trackList(j).attach > 0
                  mask = rMask;
                else
                  mask = lMask;
                end
                if ~isempty(planeFit(t).planeVectors)
                  % Rotate mask into coordinate system.
                  angle = acos(planeFit(t).planeVectors(1,1))*180/pi;
                  mask = imrotate(mask, angle, 'nearest', 'crop');
                else
                  if ~maskWarning
                    warning('No rotation vector for mask');
                    maskWarning=1;
                  end
                end
              end
              % Coordinates are in image system, so plot(x,y) should draw the
              % spots in the correct place over the image. However, the image
              % matrix is indexed by (row,col) => (y,x).
              maskImg = mask .* imgPlane(x-r:x+r, y-r:y+r);
              nonNanPix = maskImg(~isnan(maskImg));
              intStruct(t).intensity(j,iChan) = mean(nonNanPix);
              intStruct(t).intensity_median(j,iChan) = median(nonNanPix);
              intStruct(t).intensity_max(j,iChan) = max(nonNanPix);
              intStruct(t).intensity_min(j,iChan) = min(nonNanPix);
              intStruct(t).intensity_ratio(j,iChan) = ...
                  intStruct(t).intensity_max(j,iChan)/intStruct(t).intensity_min(j,iChan);
              intStruct(t).maskCoord(j,:) = maskCoords;
              %[maxX,maxY] = ind2sub(size(maskImg),maxIdx);
              [~,maxIdx] = max(maskImg(:));
              [maxY,maxX] = ind2sub(size(maskImg),maxIdx);
              intStruct(t).maxCoord(j,:) = [maxX+x-r-1,maxY+y-r-1,z];

              % Calculate angle between spot and max point.
              vector = pixCoords-intStruct(t).maxCoord(j,:);
              intStruct(t).distToMax(j) = norm(vector(1:2),2)*metadata.pixelSize(1);
              % Rotate angle into coordinate system.
              if size(pixCoords,2)==size(planeFit(t).planeVectors,2)
                vector = vector*planeFit(t).planeVectors;
              end
              angle = -atan2(vector(2),-vector(1)); % See doc atan2 for diagram.
              intStruct(t).angleToMax(j) = angle;
              intStruct(t).referenceChannel = refChan;
              if useTracks
                intStruct(t).referenceStruct = 'trackList';
              else
                intStruct(t).referenceStruct = 'initCoord';
              end
            end
          end
        end

      % Report progress.
      prog = kitProgress((t/nFrames)*(iChan/nChans), prog);
    end

    if opts.photobleachCorrect && nFrames>1
      % Compute photobleach from entire image.
      pbProfile=kitIntensityDistn(job,reader,metadata,channel,[],[],1,job.ROI.crop);
      t1=((1:size(pbProfile,1))-1)';

      % Normalize PB if necessary.
      if pbProfile(1) ~= 1 || ~all(pbProfile(:) < 1)
        pbProfile(:) = pbProfile(:)/pbProfile(1);
      end

      if numel(t1) > 4 && license('test','Curve_Fitting_Toolbox')
        % Fit double exp.
        pbF = fit(t1,pbProfile,'exp2');
        % Correct.
        for t=1:nFrames
          intStruct(t).intensity(:,iChan) = intStruct(t).intensity(:,iChan)/pbF(t1(t));
          intStruct(t).intensity_median(:,iChan) = intStruct(t).intensity_median(:,iChan)/pbF(t1(t));
          intStruct(t).intensity_max(:,iChan) = intStruct(t).intensity_max(:,iChan)/pbF(t1(t));
          intStruct(t).intensity_min(:,iChan) = intStruct(t).intensity_min(:,iChan)/pbF(t1(t));
        end
      else
        kitLog('Warning: Not correcting for photobleach');
      end

    end
end
    
% record background intensity
cellInt.back = back;
cellInt.backMode = backMode;
cellInt.backDiff = backDiff;
cellInt.maskPixelRadius = r;
cellInt.intensityDistF = intensityDistF;
cellInt.intensityDistX = intensityDistX;

% return data
job.dataStruct{channel}.cellInt = cellInt;
job.dataStruct{channel}.spotInt = intStruct;

end

%%
function movie=kitReadWholeMovie_edit(imageReader,metadata,c,crop,normalizePlanes,normalize,flip180)
% KITREADIMAGESTACK Read a whole movie from a movie file
%
%    MOVIE = KITREADIMAGESTACK(IMAGEREADER,METADATA,C,CROP,NORMALIZE) Read a
%    whole movie in channel C from IMAGEREADER described by METADATA.
%
%    CROP Optional, vector of [XMIN,YMIN,WIDTH,HEIGHT] for cropping stack.
%
%    NORMALIZEPLANES Optional, 0, 1 or -1. Normalize by maximum pixel value. Defaults
%    to 0. If -1, no normalization is performed and image is returned in
%    original datatype, otherwise it is converted to double.
%
%    NORMALIZE Optional, 0 or 1. Normalize by maximum pixel entire movie. Defaults to 0.
%
%    flip180 whether to rotate channel 180deg (CCC for MPS1 data 20250828)
%
%    Alternatively, IMAGEREADER can be a string filename as a shortcut.
%
% Copyright (c) 2013 Jonathan W. Armond
if nargin<3
  c = 1;
end

if nargin<4
  crop = [];
end

if nargin<5
  normalizePlanes = 0;
end
if nargin<6
  normalize = 0;
end

if ischar(imageReader)
  % imageReader is a filename.
  [metadata,rdr] = kitOpenMovie(imageReader);
  movie = kitReadWholeMovie_edit(rdr,metadata,c,crop,normalizePlanes,normalize,flip180);
  rdr.close();
  return
end

if normalizePlanes == -1
  dataType = metadata.dataType;
else
  dataType = 'double';
end

stackSize = kitComputeStackSize(crop,metadata.frameSize);

movie = zeros([stackSize, metadata.nFrames], dataType);
for t = 1:metadata.nFrames
  movie(:,:,:,t) = kitReadImageStack_edit(imageReader, metadata, t, c, crop, normalizePlanes, flip180);
end

if normalize>0
  movie = movie/max(movie(:));
end
end

%%
function stack=kitReadImageStack_edit(imageReader,metadata,t,c,crop,normalize, flip180)
% KITREADIMAGESTACK Read a single image frame stack from a movie file
%
%    STACK = KITREADIMAGESTACK(IMAGEREADER,METADATA,T,C,CROP,NORMALIZE) Read
%    a single image frame stack (i.e. all z-planes) at time T in channel C from
%    IMAGEREADER described by METADATA.
%
%    CROP Optional, vector of [XMIN,YMIN,WIDTH,HEIGHT] for cropping stack.
%
%    NORMALIZE Optional, 0, 1 or -1. Normalize by maximum pixel value. Defaults
%    to 1. If -1, no normalization is performed and image is returned in
%    original datatype, otherwise it is converted to double.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<5
  crop = [];
end

if nargin<6
  normalize = 0;
end

if normalize == -1
  dataType = metadata.dataType;
else
  dataType = 'double';
end

stackSize = kitComputeStackSize(crop,metadata.frameSize);

stack = zeros(stackSize, dataType);
for z = 1:metadata.frameSize(3)
  stack(:,:,z) = kitReadImagePlane_edit(imageReader, metadata, t, c, z, crop, normalize, flip180);
end

end

%%
function plane=kitReadImagePlane_edit(imageReader,metadata,t,c,z,crop,normalize, flip180)
% KITREADIMAGEPLANE Read a single image plane from a movie file.
%
%    PLANE = KITREADIMAGEPLANE(IMAGEREADER,METADATA,T,C,Z,CROP,NORMALIZE) Read a
%    single image plane from a movie file associated with open IMAGEREADER
%    and described by METADATA. Returns plane Z, at timepoint T, in channel C
%    (all one-based indexes).
%
%    CROP Optional, vector of [XMIN,YMIN,WIDTH,HEIGHT] for cropping.
%
%    NORMALIZE Optional, 0, 1 or -1. Normalize by maximum pixel value. Defaults
%    to 0, which converts datatype to double but doesn't normalize. If -1, no
%    normalization is performed and image is returned in original datatype.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<6
  crop = [];
end

if nargin<7
  normalize = 0;
end

iPlane = imageReader.getIndex(z-1,c-1,t-1) + 1;
plane = bfGetPlane(imageReader,iPlane)';

% Crop if requested.
if ~isempty(crop)
  plane = imcrop(plane,crop);
end

% Convert to double.
if normalize >= 0
  if strcmp(metadata.dataType,'int8')
    % im2double doesn't support int8 for some reason.
    plane = int16(plane)*256;
  end
  plane = im2double(plane);
  if metadata.isFloatingPoint
    % Scale to [0,1] by dividing by 16-bit integer range.
    plane = plane / (2^16-1);
  end

  % Normalize by max pixel.
  if normalize == 1
    plane = plane / max(plane(:));
  end
end

if flip180
    plane = rot90(plane,2);
    %rotates plane by 90deg twice ie by 180deg
end

end