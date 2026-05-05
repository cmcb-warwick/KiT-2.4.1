function kitFilterSpots(jobset,varargin)
% KITFILTERSPOTS(JOBSET,...) Displays a GUI to allow removal of erroneous
% spots from a JOBSET.
%
% Copyright (c) 2018 C. A. Smith
% Modified by C. C. Conway (2022)

% Define colours for rectangles.
scriptOpts.startMovie = 1; %CCC 2024.01.17
scriptOpts.suggestExcl = 0;
scriptOpts = processOptions(scriptOpts,varargin{:}); %CCC 2024.01.17

handles.col = [1 0   0;...
               0 0.75 0];
           
% Get the data.
job = kitLoadAllJobs(jobset);
handles.nMovs = length(job);
handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),jobset.options.spotMode));

% Get basic information from the jobset.
opts = jobset.options;
handles.chanID = opts.coordSystemChannel;

nFrames = job{1}.metadata.nFrames;
if nFrames > 1
  error('This function cannot yet be used for movies.')
end

% Predefine some handles required during looping
handles.movID = scriptOpts.startMovie; %CCC 2024.01.17, change from 1
if handles.movID == 1 %CCC 2024.01.17
    handles.prevEnable = 'off';
else %CCC 2024.01.17
    handles.prevEnable = 'on'; %CCC 2024.01.17
end %CCC 2024.01.17
handles.nextEnable = 'on';
handles.nextChan = 0;
handles.changeZoom = 0;
handles.zoom = 1;

% Start progress.
prog = kitProgress(0);

% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
        
    % get this image ID
    iMov = handles.movID;
    % get this channel ID
    iChan = handles.chanID;

    % check whether there is any data contained within this movie
    if ~isfield(job{iMov},'dataStruct') || ~isfield(job{iMov}.dataStruct{iChan},'failed') || job{iMov}.dataStruct{iChan}.failed
        fprintf('Analysis failed, or no data contained in this movie, so skipping to next movie \n');
        if handles.movID == handles.nMovs %CCC 2024.01.17
            handles.stop = 1; %CCC 2024.01.17
        else %CCC 2024.01.17
            handles.movID = handles.movID+1;
        end %CCC 2024.01.17
        continue
    end
    % get dataStruct
    dS = job{iMov}.dataStruct{iChan};
    % get the full initCoord and spotInt
    iC = dS.initCoord;
    % get IDs for all spots not previously filtered
    nonNaNs = find(~isnan(iC.allCoord(:,1)));
    
    % back up the full initCoord and spotInt
    if isfield(dS,'rawData')
      iC = dS.rawData.initCoord;
      if isfield(dS.rawData,'spotInt')
        sI = dS.rawData.spotInt;
      end
    elseif isfield(dS,'spotInt')
        sI = dS.spotInt;
    end
    raw.initCoord = iC;
    if exist('sI','var')
        raw.spotInt = sI;
    end
    dS.rawData = raw;
    job{iMov}.dataStruct{iChan} = dS;

    % if wanting to suggest spots to omit based on less than 25% of max 20%
    % of KT amplitude intensity
    if scriptOpts.suggestExcl
        if handles.chanID == opts.coordSystemChannel
            sortedAmps = sort(rmmissing(dS.rawData.initCoord.amp(:,1)), 'descend');
            nTop20pct = ceil(length(sortedAmps)/5);
            intTop20pct = mean(sortedAmps(1:nTop20pct), 'omitnan');
            minIntensity_25pct = intTop20pct*0.25;
            belowThreshold = find(dS.rawData.initCoord.amp(:,1) < minIntensity_25pct);
            intFilteredNonNaNs = setdiff(nonNaNs, belowThreshold);
        end
    end
    
    % get number of spots
    nSpots = size(dS.initCoord.allCoord,1);
    handles.nSpots = nSpots;
    
    % show all spots - defined using tracks
    rectDims = griddedSpots(job{iMov},'channel',handles.chanID,'zoomed',handles.zoom);
    
    % get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    if handles.nextChan || handles.changeZoom
      keptSpots = handles.keep .* [1:nSpots];
      handles.keep = ismember(keptSpots,nonNaNs);
    elseif scriptOpts.suggestExcl
      handles.keep = ismember(1:nSpots,intFilteredNonNaNs);
    else
      handles.keep = ismember(1:nSpots,nonNaNs);
    end
    %handles.keep
    % reset channel information if necessary
    handles.nextChan = 0;
    handles.changeZoom = 0;

    % draw rectangles
    hold on
    for iSpot = 1:nSpots
        % get the colour for this spot
        keep = handles.keep(iSpot);
        icol = handles.col(keep+1,:);
        % draw the rectangle
        rectangle('Position',[rectPos(iSpot,:)-0.5 rectWid rectWid],...
            'EdgeColor',icol,'LineWidth',3);
    end
    
    % Buttons and labels.
    btnw = [12 7]; btnh = 2; h = 1.5;
    figpos = get(gcf,'Position');
    dx = 2.5; ddx = 1;
    % make label at top left for instructions
    x = dx; y = figpos(4)-(btnh+ddx/2);
    labw = 60;
    handles.instructions = label(gcf,'Click on spots to keep (green) or remove (red).',[x y labw h],12);
    % add all buttons: finish, next and previous
    x = figpos(3)-(btnw(1)+dx);
    handles.finishBtn = button(gcf,'Finish',[x y btnw(1) btnh],@finishCB);
    x = x-(btnw(2)+ddx);
    handles.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextMovCB);
    handles.nextBtn.Enable = handles.nextEnable;
    x = x-(btnw(2)+ddx/2);
    handles.prevBtn = button(gcf,'Prev',[x y btnw(2) btnh],@prevMovCB);
    handles.prevBtn.Enable = handles.prevEnable;
    % deselect all
    x = figpos(3)-(btnw(1)+dx); y = ddx;
%     handles.deselectBtn = button(gcf,'Deselect all',[x y btnw(1) btnh],@deselectAllCB);
    handles.invertBtn = button(gcf,'Invert all',[x y btnw(1) btnh],@invertCB);
    x = x-(btnw(1)+dx); y = ddx;
    handles.nextChanBtn = button(gcf,'Next chan',[x y btnw(1) btnh],@nextChanCB);
    x = x-((btnw(1)+dx)*2); y = ddx;
    handles.intCircleBtn = button(gcf,'Show/hide intensity radius',[x y btnw(1)*2+dx btnh],@showIntRadCB);
    x = x-(btnw(1)+dx); y = ddx;
    handles.changeZoomBtn = button(gcf,'Change zoom',[x y btnw(1) btnh],@changeZoomCB);
    x = x-(btnw(1)+dx); y = ddx;
    handles.nKTs = label(gcf,sprintf('KTs = %d', length(find(handles.keep))),[x y btnw(1) h],12);


    
    % set up remove environment
    set(get(gca,'Children'),'ButtonDownFcn',@rmvCB);
    
    % GUI now set up, wait for user
    uiwait(gcf);
    
    % check whether user asked to switch channel
    if handles.nextChan
        % close the figure for the next channel
        close(gcf);
        continue
    end

    if handles.changeZoom
        close(gcf);
        continue
    end
    % get the spots requiring removal
    rmvList = find(~handles.keep);
    
    % check for any sisters not in the list, provide warning and
    % remove if so
    incorrect = setdiff(rmvList,1:nSpots);
    if ~isempty(incorrect)
        warning('The following selected spots do not exist: %s. Will ignore.',num2str(incorrect));
    end

    % process the list
    if ~isempty(rmvList)
      for jChan = handles.chans
        % push all removed spots to NaNs
        dS = job{iMov}.dataStruct{jChan};
        if isfield(dS,'rawData')
            iC = dS.rawData.initCoord;
        else
            iC = dS.initCoord;
        end
        dS.rawData.initCoord = iC;
        % push to initCoord
        iC.allCoord(rmvList,:) = NaN;
        iC.allCoordPix(rmvList,:) = NaN;
        iC.nSpots = sum(~isnan(iC.allCoord(:,1)));
        iC.amp(rmvList,:) = NaN;
        iC.bg(rmvList,:) = NaN;
        
        % back up results
        dS.initCoord = iC;
        job{iMov}.dataStruct{jChan} = dS;
      end
    
      % process the list
      for jChan = find(job{iMov}.options.intensity.execute)
        % push all removed spots to NaNs
        dS = job{iMov}.dataStruct{jChan};
        if isfield(dS,'rawData')
            if isfield(dS.rawData,'spotInt')
                sI = dS.rawData.spotInt;
            else
                sI = dS.spotInt;
            end
        else
            sI = dS.spotInt;
        end
        dS.rawData.spotInt = sI;
        % push to spotInt
        sI.intensity(rmvList,:) = NaN;
        %sI.intensity_mean(rmvList,:) = NaN;
        sI.intensity_median(rmvList,:) = NaN;
        sI.intensity_min(rmvList,:) = NaN;
        sI.intensity_max(rmvList,:) = NaN;
        sI.intensity_ratio(rmvList,:) = NaN;
        
        % back up results
        dS.spotInt = sI;
        job{iMov}.dataStruct{jChan} = dS;
      end
      
    end
    
    % save results
    job{iMov} = kitSaveJob(job{iMov});

    % update progress
    prog = kitProgress(iMov/handles.nMovs,prog);
    
    % close the figure for the next movie
    close(gcf);

end

kitLog('Manual filtering complete.');

% Re-run plane fitting for jobset with new filtered data.
kitLog('Re-fitting planes to filtered data.');
kitRunJob(jobset,'existing',1,'tasks',[2 6]);

%% Callback functions

function rmvCB(hObj,event)
  
  % get the position of the click
  pos=get(gca,'CurrentPoint');
  xpos = pos(1,1); ypos = pos(1,2);

  % get all positions
  allPos = handles.rectPos;
  
  % get candidates using click's row position
  diffs = xpos-allPos(:,1);
  diffs(diffs<0) = NaN;
  xidx = find(diffs == nanmin(diffs));
  % get candidates using click's column position
  diffs = ypos-allPos(:,2);
  diffs(diffs<0) = NaN;
  yidx = find(diffs == nanmin(diffs));
  % get the common candidate
  idx = intersect(xidx,yidx);

  % if a click is made elsewhere, remind user how to select images
  if isempty(idx)
    handles.instructions.String = 'Click on the images to select/deselect.';
    return
  end
  
  % get the colour for this spot
  keepStat = ~handles.keep(idx);
  icol = handles.col(keepStat+1,:);

  % draw the rectangle
  rectangle('Position',[handles.rectPos(idx,:)-0.5 handles.rectWid handles.rectWid],...
      'EdgeColor',icol,'LineWidth',3);

  handles.keep(idx) = keepStat;
  handles.nKTs.String = sprintf('KTs = %d', length(find(handles.keep)));
end

function deselectAllCB(hObj,event)
  hs = handles;
  % force all stops to be ignored
  handles.keep(1,:) = 0;
  for i = 1:hs.nSpots
    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',hs.col(1,:),'LineWidth',3);
  end
end

function invertCB(hObj,event)
  hs = handles;
  % force all stops to be ignored
  handles.keep = ~handles.keep;
  for i = 1:hs.nSpots
    
    % get the colour for this spot
    keepStat = handles.keep(i);
    jcol = handles.col(keepStat+1,:);
    
    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',jcol,'LineWidth',3);
  end
  handles.nKTs.String = sprintf('KTs = %d', length(find(handles.keep)));
end

function prevMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID-1;
  if handles.movID == 1
    handles.prevEnable = 'off';
  end
  handles.nextEnable = 'on';
  handles.nextChan = 0;
  handles.changeZoom = 0;
  handles.zoom = 1;
  % continue the function
  uiresume(gcf);
end

function nextMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID+1;
  handles.prevEnable = 'on';
  if handles.movID == handles.nMovs
    handles.nextEnable = 'off';
  end
  handles.nextChan = 0;
  handles.chanID = opts.coordSystemChannel;
  handles.changeZoom = 0;
  handles.zoom = 1;
  % continue the function
  uiresume(gcf);
end

function nextChanCB(hObj,event)
  % update the handles
  idx = find(handles.chanID==handles.chans);
  if idx == length(handles.chans)
    handles.chanID = handles.chans(1);
  else
    handles.chanID = handles.chans(idx+1);
  end
  handles.nextChan = 1;
  handles.changeZoom = 0;
  handles.zoom = 1;
  % continue the function
  uiresume(gcf);
end

function changeZoomCB(hObj,event)
  % update the handles
  idx = handles.zoom;
  handles.zoom = abs(idx-1);
  handles.changeZoom = 1;
  % continue the function
  uiresume(gcf);
end

function showIntRadCB(hObj,event)
  % update the handles
  lineobjs = findobj('Type', 'Line');
  for i = 1:length(lineobjs)
    lineobjs(i).Visible = ~lineobjs(i).Visible;
  end
  % continue the function
  %uiresume(gcf);
end

function finishCB(hObj,event)
  % force stop
  handles.stop = 1;
  % continue the function
  uiresume(gcf);
end

end


    
