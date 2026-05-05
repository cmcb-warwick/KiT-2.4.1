function selectedData = kitSelectData_SelChan_sS(job,varargin)
% KITSELECTDATA_SELCHAN_SS lets you select spots from a movie structure
% with an existing spot selection in one channel.
%
%   KITSELECTDATA_SELCHAN_SS(JOB, ...) creates a spot selection structure
%   using an initial spot selection to pre-exclude spots from
%   consideration. This creates a sS_sel structure. JOB should be full
%   movieStructure. Sister pairing must be completed for JOB.
%   
%   Options, defaults in {}:
%   
%   addMore: {[]} or an existing spot selection (sS_sel) for this channel
%       that you want to add more spots to. Only use if you're part way
%       through creating a sS_sel for this channel. Distinct from prevSel.
%       Should have been created using full movieStructure.
%
%   prevSel: {[]}, a spot selection (sS_sel) that you created in full for
%       another channel. Should have been created using full
%       movieStructure (haven't tested using mS_sel).
% 
%   startMovie: {1} or number. Which image to start selections from. Best
%       used in combination with an addMore structure created up to the
%       previous image.
% 
%   channel: 1, 2, or {3}. Which channel to perform this round of spot
%       selection in.
% 
%   dataType: {'spots'} or 'sisters'. Whether to assess individual spots or
%       sister pairs for inclusion. 'sisters' hasn't been tested yet.
% 
%   contrast: {[0.1 1]} or similar two-element vector. Range over which to
%       contrast images. Tips:
%           - Increase brightness by changing to [0.1 0.9]
%           - Decrease background noise by changing to [0.5 1]
%
%   Adapted from kitSelectData from KiT v2.4.0 and kitSelectData in KiT
%   v2.1.10 (both C. A. Smith).
%
%   Copyright (c) 2023 C. C. Conway


% Get the data.
opts.addMore = [];
opts.prevSel = [];
opts.contrast = [0.1 1];
opts.dataType = 'spots'; % can also be 'sisters'
opts.startMovie = 1;
opts.channel = 3;
opts = processOptions(opts,varargin{:});

% check for previous selections, set up new structure if none given
if isempty(opts.addMore)
    selectedData.dataType = opts.dataType;
    selectedData.selection = {[]};
else
    selectedData = opts.addMore;
    % ensure that dataTypes match
    if ~strcmp(selectedData.dataType,opts.dataType)
      error('Cannot add more data to the provided selection due to clash in dataType.');
    end
end

allSels = selectedData.selection;

if isempty(opts.prevSel)
    orgSels = [];
else
    orgSels = opts.prevSel.selection{1};
end

nMovs = length(job);
for iMov = opts.startMovie:nMovs
    movtxt = sprintf('Cell %d. Keep selecting?', iMov);
    keepSelecting = questdlg(movtxt, 'Keep selecting?', ...
	'Yes','Skip','No','Yes');
    % Handle response
    switch keepSelecting
        case 'Yes'
            kitLog('Spot selection in movie %i', iMov)
            orgCellIdx = find(orgSels(:,1) == iMov);
            if ~isempty(orgCellIdx)
                orgAccepted = orgSels(orgCellIdx,2);
            else
                orgAccepted = [];
            end
            tempList = selectSpots(job{iMov}, opts.channel, orgAccepted);
            if isempty(tempList)
                continue
            else
                nData = length(tempList);
                tempStruct = ones(nData,2)*iMov;
                tempStruct(:,2) = tempList;
                % collate selections
                allSels{1} = [allSels{1}; tempStruct];
            end
        case 'Skip'
            kitLog('Skipping spot selection in movie %i', iMov)
            continue
        case 'No'
            kitLog('Quitting. Spot selection saved to movie %i', iMov-1)
            break
    end

    assignin('base', 'workingSpotSelection', allSels)
    selectedData.selection = allSels;

    currentFolder = pwd;
    folderSeparator = filesep;
    todayDate = datestr(datetime('now'), 'yyyymmdd');

    filename = sprintf('%s%s%s_workingSpotSelection.mat', currentFolder, folderSeparator, todayDate);
    save(filename,'selectedData');
    %this is just to keep your spot selections safe in case MATLAB crashes
    %before you're finished (this has happened to me too often...)
    
    kitLog('Updated save file up to experiment %i of %i.',iMov,nMovs);
end

selectedData.selection = allSels;

end





%%
function tempList = selectSpots(job,channel,orgSelection)
% SELECTSPOTS(JOBSET,...) Displays a GUI to allow exclusion of erroneous
% spots from a spot selection.
%
% Copyright (c) 2018 C. A. Smith
% Modified by C. C. Conway (2022)

% Define colours for rectangles.
handles.col = [1 0   0;...
               0 0.75 0];
           
handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),job.options.spotMode));



% Get basic information from the jobset.
opts_more = job.options;
handles.chanID = channel;

nFrames = job.metadata.nFrames;
if nFrames > 1
  error('This function cannot yet be used for movies.')
end

% Predefine some handles required during looping
handles.nextChan = 0;



% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
        
    % get this channel ID
    iChan = handles.chanID;

    % check whether there is any data contained within this movie
    if ~isfield(job,'dataStruct') || ~isfield(job.dataStruct{iChan},'failed') || job.dataStruct{iChan}.failed
        fprintf('Analysis failed, or no data contained in this movie, so skipping to next movie \n');
        break
    end
    % get dataStruct
    dS = job.dataStruct{iChan};

    
    % get number of spots
    nSpots = length(dS.tracks);
    handles.nSpots = nSpots;
    if ~isfield(job.dataStruct{iChan}, 'tracks') || job.dataStruct{iChan}.tracks(1).tracksFeatIndxCG == 0
        kitLog('No pairs found')
        break
    end
    % show all spots - defined using tracks
    [rectDims, gsCoords] = griddedSpots_sS(job,'channel',iChan);
    
    % get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    nonNaNs = find(~isnan(gsCoords(:,1)));
    if handles.nextChan
      keptSpots = handles.keep .* [1:nSpots];
      handles.keep = ismember(keptSpots,nonNaNs);
    else
        nonNaNMembers = ismember([1:nSpots],nonNaNs);
        %nonNaNSpots = nonNaNMembers .* [1:nSpots];
        includedMembers = ismember([1:nSpots],orgSelection);
        %includedSpots = includedMembers .* [1:nSpots];
      handles.keep = nonNaNMembers .* includedMembers;
    end
    %handles.keep
    % reset channel information if necessary
    handles.nextChan = 0;


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

    % deselect all
    x = figpos(3)-(btnw(1)+dx); y = ddx;
%     handles.deselectBtn = button(gcf,'Deselect all',[x y btnw(1) btnh],@deselectAllCB);
    handles.invertBtn = button(gcf,'Invert all',[x y btnw(1) btnh],@invertCB);
    x = x-(btnw(1)+dx); y = ddx;
    handles.nextChanBtn = button(gcf,'Next chan',[x y btnw(1) btnh],@nextChanCB);
    x = x-((btnw(1)+dx)*2); y = ddx;
    handles.intCircleBtn = button(gcf,'Show/hide intensity radius',[x y btnw(1)*2+dx btnh],@showIntRadCB);

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

end
close(gcf);

if ~isfield(job,'dataStruct') || ~isfield(job.dataStruct{iChan},'failed') ...
        || job.dataStruct{iChan}.failed || ~isfield(job.dataStruct{iChan}, 'tracks') ...
        || job.dataStruct{iChan}.tracks(1).tracksFeatIndxCG == 0
        tempList = [];
else
    keptList = handles.keep .* [1:nSpots];
    tempList = intersect(keptList,[1:nSpots]);
end
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



