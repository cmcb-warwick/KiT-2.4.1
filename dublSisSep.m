function compiledIntra = dublSisSep(movies,varargin)
% DUBLSISSEP Produces a structure of population-level inter-sister
% distance measurements over multiple experiments. Use if 1) only a single
% channel has detectable spots in all cells/kinetochores or 2) if you want
% to measure inter-sister distance when your reference kinetochore marker
% is not the 'inner' channel.
%
%    DUBLSISSEP(MOVIES,...) Provides coordinates and inter-sister
%    measurements for all sisters within all cells across all experiments.
%    The resulting structure provides the tools for plotting inter-sister
%    pseudo-timeline graphs. Options are available.
%
%    Options, defaults in {}:-
%
%    category: {[]} or string. The category from which to calculate intra-
%       measurements.
%
%    channel: {[2]} or number from 1 to 3. The channel that inter-sister 
%       distance will be measured in.
%
%    paired: 0 or {1}. Whether or not to take paired measurements, or raw
%       spot-by-spot measurements.
%
%    prevMeas: {[]} or a structure previously generated. A structure of
%       results from previous experiments to allow new experiment data to
%       be appended.
%
%    spotSelection: {[]} or output from kitSelectData. A structure
%       containing a selection of either sister pair or track IDs per
%       movie, for each experiment. Allows for only specific sisters or
%       spots to be included in the data collection.
%
%
% Copyright (c) 2018 C. A. Smith, edited 2026 C. C. Conway

% default options
opts.category = [];
opts.centralise = 1;
opts.channel = 2;
opts.paired = 1;
opts.prevMeas = [];
opts.spotSelection = [];
% user options
opts = processOptions(opts,varargin{:});

%% Pre-processing input structure

%check structure of movies
if ~iscell(movies{1})
    movies = {movies};
    kitLog('Movie structure provided implies only one experiment. Assuming only one experiment.');
end
%find number of movies
numExpts1 = length(movies);

%process input so that all structs are in cell format
if isempty(opts.spotSelection)
  subset = repmat({[]},numExpts1,1);
  selType = 0;
elseif isstruct(opts.spotSelection) && isfield(opts.spotSelection,'dataType')
  subset = opts.spotSelection.selection;
  switch opts.spotSelection.dataType
    case 'spots' %tracks
      selType = 1;
    case 'sisters'
      selType = 2;
    case 'initCoord'
      selType = 3;
  end
else
  kitLog('Provided spotSelection structure was not derived from kitSelectData. No selection will be imposed.');
  subset = repmat({[]},numExpts1,1);
  selType = 0;
end

%find number of movies and sisters, and ensure they match
numExpts2 = length(subset);
if numExpts1 ~= numExpts2
  error('Have %i spot selections for %i experiments. Please provide spot selection for each experiment.',numExpts2,numExpts1)
end
numExpts = numExpts1;

%% Pre-processing output structure

% assume no pairing until we find one movie that is
paired = 0;
if opts.paired
  for iExpt = 1:numExpts
    for iMov = 1:length(movies{iExpt})
      paired = isfield(movies{iExpt}{iMov}.dataStruct{opts.channel},'sisterList');
      if paired; break; end
    end
    if paired; break; end
  end
end
if selType==3
    if paired
        subset = opts.spotSelection.selection;
        selType = 1;
    else
        subset = opts.spotSelection.rawSelection;
    end
end

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = dublSingleChanIntraStructure(paired);
    allData = struct2strForm(allData);
    
else
    % get all old data
    allData = struct2strForm(opts.prevMeas);
    
end

% predesignate error arrays
noFail = [];
noSis = [];
noSkip = [];

%% Compiling measurements

for iExpt = 1:numExpts
    
    % get movie and sister list
    theseMovies = movies{iExpt};
    iSubset = subset{iExpt};
    

    
    for iMov = 1:length(theseMovies)
        
      % get the movie index
      movNum = theseMovies{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseMovies{iMov},'dataStruct')
        noFail = [noFail; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dSinner = theseMovies{iMov}.dataStruct{opts.channel};
      
      % check whether the movie failed
      if (isfield(dSinner,'failed') && dSinner.failed) || ~isfield(dSinner,'initCoord')
        noFail = [noFail; iExpt movNum];
        continue
      end
      
      % check whether the user skipped this movie
      if (isfield(theseMovies{iMov},'keep') && ~theseMovies{iMov}.keep)
        noSkip = [noSkip; iExpt movNum];
        continue
      end
      
      % get basic metadata
      nFrames = theseMovies{iMov}.metadata.nFrames;
      pixelSize = theseMovies{iMov}.metadata.pixelSize;
      % check if there is a plate fit
      plane = (isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit) && ~isempty(dSinner.planeFit.planeVectors));
      
      % get initCoord structure
      
      iCinner = dSinner.initCoord;  

      if paired
        
        % check whether a sisterList is present, and if it contains any sisters
        if ~isfield(dSinner,'sisterList') || isempty(dSinner.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % if no sisters given, go through all sisters in movie
        switch selType
            case 0 % no spot selection
                iSubset = 1:length(theseMovies{iMov}.dataStruct{opts.channel}.sisterList);
            case 1 % using spots/tracks
                iSubset = 1:length(theseMovies{iMov}.dataStruct{opts.channel}.sisterList);
                theseTracks = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
            case 2 % using sisters
                iSubset = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
        end

        % check that there are sisters
        if isempty(dSinner.sisterList(1).trackPairs)
            continue
        end
        
        for iSis = iSubset
            
            % construct sister pair label
            label = sprintf('%02d%02d%03d',iExpt,iMov,iSis);
            
            % start counter for storing data
            c=1;

            % get sisterLists
            sLinner = dSinner.sisterList(iSis);
            
            % get trackID and spotIDs, make spotIDs nan if deselected
            trackIDs = dSinner.sisterList(1).trackPairs(iSis,1:2);
            spotIDs = nan(nFrames,2);
            for iTrack = 1:2
                if selType~=1 %none or sisters
                    spotIDs(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).featIndx;
                else %tracks
                    if ismember(trackIDs(iTrack),theseTracks)
                        spotIDs(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).featIndx;
                    else
                        trackIDs(iTrack) = NaN;
                        switch iTrack
                            case 1
                                sLinner.coords1(:,1:3) = NaN;
                            case 2
                                sLinner.coords2(:,1:3) = NaN;
                        end
                    end
                end 
            end
            % filter based on chosen category
            if ~isempty(opts.category) && nFrames==1
                if isfield(theseMovies{iMov},'categories') && ...
                        isfield(theseMovies{iMov}.categories,opts.category)
                    spotIDs = intersect(spotIDs,theseMovies{iMov}.categories.(opts.category));
                else
                    trackIDs = [NaN NaN];
                end
            end
            % if both spots skipped
            if all(isnan(trackIDs))
                continue;
            end

            %% Kinetochore positions
            
            % get microscope coordinates, and plane coordinates if present
            mCoordsInner = nan(nFrames,6);
            pCoordsInner = nan(nFrames,6);
            for iFrame = 1:nFrames
                for iTrk = 1:2
                    % track one stored in (:,1:3), two in (:,4:6)
                    rng = (3*(iTrk-1)+1):3*iTrk;
                    if ~isnan(spotIDs(iFrame,iTrk))
                        mCoordsInner(iFrame,rng) = iCinner(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3);
                        % check whether or not this movie has a planeFit
                        if plane
                            % rotate coordinates into plane
                            pF = dSinner.planeFit;
                            if ~isempty(pF(iFrame).planeVectors)
                                coordSystem = pF(iFrame).planeVectors;
                                pCoordsInner(iFrame,:) = mCoordsInner(iFrame,:) - repmat(pF.planeOrigin,1,2);
                                
                                pCoordsInner(iFrame,1:3) = (coordSystem\(pCoordsInner(iFrame,1:3)'))';
                                pCoordsInner(iFrame,4:6) = (coordSystem\(pCoordsInner(iFrame,4:6)'))';

                            end
                        end
                    end
                end
            end
            
            % put data into string format
            newData(c,:) = {'label',label}; c=c+1;
            
            % get microscope coordinates of each spot
            mCoords_x = [mCoordsInner(:,[1 4])];
            mCoords_y = [mCoordsInner(:,[2 5])];
            mCoords_z = [mCoordsInner(:,[3 6])];
            % put data into string format
            newData(c,:) = {'microscope.coords.x',mCoords_x}; c=c+1;
            newData(c,:) = {'microscope.coords.y',mCoords_y}; c=c+1;
            newData(c,:) = {'microscope.coords.z',mCoords_z}; c=c+1;
            
            % get plate coordinates of each spot
            pCoords_x = [pCoordsInner(:,[1 4])];
            pCoords_y = [pCoordsInner(:,[2 5])];
            pCoords_z = [pCoordsInner(:,[3 6])];
            % put data into string format
            newData(c,:) = {'plate.coords.x',pCoords_x}; c=c+1;
            newData(c,:) = {'plate.coords.y',pCoords_y}; c=c+1;
            newData(c,:) = {'plate.coords.z',pCoords_z}; c=c+1;
            
            %% Inter- and intra-kinetochore measurements
            
            micrData = pairedMeasurements(mCoordsInner,0);
            % put data into string format
            newData(c,:) = {'microscope.sisSep.x',micrData.sisSep_x}; c=c+1;
            newData(c,:) = {'microscope.sisSep.y',micrData.sisSep_y}; c=c+1;
            newData(c,:) = {'microscope.sisSep.z',micrData.sisSep_z}; c=c+1;
            newData(c,:) = {'microscope.sisSep.twoD',micrData.sisSep_2D}; c=c+1;
            newData(c,:) = {'microscope.sisSep.threeD',micrData.sisSep_3D}; c=c+1;

            
            plateData = pairedMeasurements(pCoordsInner,1);
            % put data into string format
            newData(c,:) = {'plate.sisSep.x',plateData.sisSep_x}; c=c+1;
            newData(c,:) = {'plate.sisSep.y',plateData.sisSep_y}; c=c+1;
            newData(c,:) = {'plate.sisSep.z',plateData.sisSep_z}; c=c+1;
            newData(c,:) = {'plate.sisSep.twoD',plateData.sisSep_2D}; c=c+1;
            newData(c,:) = {'plate.sisSep.threeD',plateData.sisSep_3D}; c=c+1;
            newData(c,:) = {'plate.twist.y',plateData.twist_y}; c=c+1;
            newData(c,:) = {'plate.twist.z',plateData.twist_z}; c=c+1;
            newData(c,:) = {'plate.twist.threeD',plateData.twist_3D}; c=c+1;
            newData(c,:) = {'plate.sisterCentreSpeed',plateData.sisCentreSpeed}; c=c+1;

            % get plate thickness measurements
            sisCentre_x = [];
            if length(dSinner.sisterList) == 1
              plateThickness = nan(1,nFrames);
            else
              for jSis = 1:length(dSinner.sisterList)
                sisCentre_x = [sisCentre_x nanmean([dSinner.sisterList(jSis).coords1(:,1) dSinner.sisterList(jSis).coords2(:,1)],2)];
              end
              plateThickness = nanstd(sisCentre_x,[],2);
            end
            % put data into string format
            newData(c,:) = {'plate.plateThickness',plateThickness}; c=c+1;


            
         
            
            %% Directional information
            
            % get direction of movement
            direc = [];
            for iTrack = 1:2
              if ~isnan(trackIDs(iTrack)) && ~isempty(dSinner.trackList(trackIDs(iTrack)).direction)
                direc(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).direction;
              else
                direc(:,iTrack) = nan(nFrames,1);
              end
            end
            direc(end+1,:) = NaN;
            
            % get potential switch events (i.e. individual timepoints between P and AP)
            switchBuffer = 4;
            switchEvent = zeros(size(direc));
            switchDirec = [diff(direc); NaN NaN];
            for iPoint = 1:2:size(switchEvent,1)-(switchBuffer+1)
                for jSis = 1:2
                    tempDirec = switchDirec(iPoint:iPoint+(switchBuffer-1),jSis);
                    idx = find(tempDirec(2:(switchBuffer-1))==0);
                    if max(abs(tempDirec))==1 && abs(sum(tempDirec))>1 && ~isempty(idx)
                        switchEvent(iPoint+idx(1):iPoint+idx(end),jSis) = 1;
                    end
                end
            end
            % calculate directional information
            direc_P = +(direc==1);               direc_P(direc_P==0) = NaN;
            direc_AP = +(direc==-1);             direc_AP(direc_AP==0) = NaN;
            direc_S = +(switchEvent==1);         direc_S(direc_S==0) = NaN;
            direc_N = +((direc+switchEvent)==0); direc_N(direc_N==0) = NaN;
            % put data into string format
            dirLabel = {'P','AP','S','N'};
            for iDir = 1:4
                eval(['newData(c,:) = {''direction.' dirLabel{iDir} ''',direc_' dirLabel{iDir} '}; c=c+1;']);
            end
        
            % compile new data with original
            allData = combineStrForms(allData,newData);
            
            % clear some data to ensure no overlap on next loop
            clear spotIDs newData direc
        
        end % sisters
        
      else
          
          % start counter for storing data
          c=1;
          
          % if no sisters given, go through all sisters in movie
          switch selType
            case 0 % no spot selection
                spotIDs = 1:size(iCinner(1).allCoord,1);
            case 1 % using spots/tracks
                trackIDs = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
                spotIDs = cat(2,dSinner.trackList(trackIDs).featIndx);
            case 3 % using initCoord
                spotIDs = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
          end
          
          % filter based on chosen category
          if ~isempty(opts.category)
              if isfield(theseMovies{iMov},'categories') && ...
                      isfield(theseMovies{iMov}.categories,opts.category)
                  spotIDs = intersect(spotIDs,theseMovies{iMov}.categories.(opts.category));
              else
                  spotIDs = [];
              end
          end
          
          % check number of spots
          nSpots = length(spotIDs);
          if nSpots == 0
              continue
          end
          
          % construct spot label
          labels = '';
          for iSpot = 1:nSpots
            labels(iSpot,:) = sprintf('%02d%02d%03d',iExpt,iMov,iSpot);
          end
          
          % put data into string format
          newData(c,:) = {'label',labels}; c=c+1;
          
          %% Kinetochore positions
          
          % get microscope coordinates, and plane coordinates if present
          mCoordsInner = iCinner(1).allCoord(spotIDs,1:3);
          % check whether or not this movie has a planeFit
          if plane
              % rotate coordinates into plane
              pF = dSinner.planeFit;
              if ~isempty(pF(1).planeVectors)
                  coordSystem = pF(1).planeVectors;
                  pCoordsInner = (coordSystem\(mCoordsInner'))';
              end
          else
              % give empty datasets the size of microscopy coordinates
              pCoordsInner = nan(size(mCoordsInner));
          end
            
          % get microscope coordinates of each spot
          mCoords_x = [mCoordsInner(:,1)];
          mCoords_y = [mCoordsInner(:,2)];
          mCoords_z = [mCoordsInner(:,3)];
          % put data into string format
          newData(c,:) = {'microscope.coords.x',mCoords_x}; c=c+1;
          newData(c,:) = {'microscope.coords.y',mCoords_y}; c=c+1;
          newData(c,:) = {'microscope.coords.z',mCoords_z}; c=c+1;
          
          % get plate coordinates of each spot
          pCoords_x = [pCoordsInner(:,1)];
          pCoords_y = [pCoordsInner(:,2)];
          pCoords_z = [pCoordsInner(:,3)];
          % put data into string format
          newData(c,:) = {'plate.coords.x',pCoords_x}; c=c+1;
          newData(c,:) = {'plate.coords.y',pCoords_y}; c=c+1;
          newData(c,:) = {'plate.coords.z',pCoords_z}; c=c+1;
          
          
          % put data into string format
 
          

          
          %% Quality control

          
          % compile new data with original
          allData = combineStrForms(allData,newData);  
          clear newData spotIDs
            
      end % paired
    end % movies     
end % expts

%% Save results to structure

compiledIntra = strForm2struct(allData);

%% Output any error information

if ~isempty(noSkip)
  fprintf('\nThe following cells were skipped by the user:\n');
  for iCell = 1:size(noSkip,1)
    fprintf('    Exp %i, Mov %i\n',noSkip(iCell,1),noSkip(iCell,2));
  end
end
if ~isempty(noFail)
  fprintf('\nThe following cells failed during spot detection:\n');
  for iCell = 1:size(noFail,1)
    fprintf('    Exp %i, Mov %i\n',noFail(iCell,1),noFail(iCell,2));
  end
end
if ~isempty(noFail)
  fprintf('\nThe following cells found no spots:\n');
  for iCell = 1:size(noFail,1)
    fprintf('    Exp %i, Mov %i\n',noFail(iCell,1),noFail(iCell,2));
  end
end
if ~isempty(noSis)
  fprintf('\nThe following cells contain no sisterList:\n');
  for iCell = 1:size(noSis,1)
    fprintf('    Exp %i, Mov %i\n',noSis(iCell,1),noSis(iCell,2));
  end
end
fprintf('\n');

end
%%
function measurements = pairedMeasurements(coordsInner,plane)
    
    if nargin<2 || isempty(plane)
        plane = 0;
    end
        
    % Inter-kinetochore: separation, twist, and velocities
    
    % coordinate-specific sister separation (using inner-kinetochore)
    sisSep_xyz = diff(reshape(coordsInner,size(coordsInner,1),3,2),[],3);
    measurements.sisSep_x = sisSep_xyz(:,1);
    measurements.sisSep_y = sisSep_xyz(:,2);
    measurements.sisSep_z = sisSep_xyz(:,3);

    % 3D sister separation
    sisSep_3D = sqrt(sum(sisSep_xyz.^2,2));
    measurements.sisSep_3D = sisSep_3D;

    % 2D sister separation
    sisSep_2D = sqrt(sum(sisSep_xyz(:,1:2).^2,2));
    measurements.sisSep_2D = sisSep_2D;

    if plane

        % coordinate-specific twist
        twist_y = sisSep_xyz(:,2)./sisSep_xyz(:,1);
        twist_y = atand(twist_y);
        twist_z = sisSep_xyz(:,3)./sisSep_xyz(:,1);
        twist_z = atand(twist_z);
        measurements.twist_y = twist_y;
        measurements.twist_z = twist_z;

        % 3D twist (dot product with the x-axis with length 1)
        xAxis = repmat([1 0 0],size(sisSep_xyz,1),1);
        twist_3D = dot(sisSep_xyz,xAxis,2);
        twist_3D = twist_3D./(sisSep_3D(:,1));
        twist_3D = acosd(twist_3D);
        twist_3D(twist_3D>90) = 180-twist_3D(twist_3D>90);
        measurements.twist_3D = twist_3D;

        % sister centre velocities in x-coordinate
        sisCentre_x = sum(coordsInner(:,[1 4]),2)/2;
        sisCentreSpeed_x = [diff(sisCentre_x); NaN];
        measurements.sisCentreSpeed = sisCentreSpeed_x;

    end


    
    

end %pairedMeasurements subfunction




%%
function intraStruct = dublSingleChanIntraStructure(paired)
% dublSingleChanIntraStructure Produces an empty structure into which
% intra-kinetochore measurements are compiled using dublIntraMeasurements.
%
%    DUBLMAKEINTRASTRUCTURE() Produces a structure allowing for both inter-
%    and intra-kinetochore measurements to be compiled in order to draw
%    population-scale analyses. No input is required.
%
%
% Copyright (c) 2016 C. A. Smith

if nargin<1 || isempty(paired)
  paired = 1;
end

intraStruct.kitVersion = kitVersion;
intraStruct.label = [];

% make substructures for all pair-derived measurements
direction = struct('P',[],'AP',[],'N',[],'S',[]);
sisSep = struct('x',[],'y',[],'z',[],...
               'twoD',[],'threeD',[]);

% make substructure for each microscope and plate coordinate systems

coords = struct('x',[],'y',[],'z',[]);

twist  = struct('y',[],'z',[],'threeD',[]);
sisterCentreSpeed = [];
plateThickness = [];

% produce generic substructure for microscope coordinate system
microscope.coords             = coords;
if paired
microscope.sisSep             = sisSep;
end

% produce same again for plate coordinate system
plate = microscope;
if paired
  plate.twist         = twist;
  plate.sisterCentreSpeed = sisterCentreSpeed;
  plate.plateThickness    = plateThickness;
  
  intraStruct.direction = direction;
end
  intraStruct.microscope = microscope;
  intraStruct.plate = plate;

end
