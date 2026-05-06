function [selectedData, spotSel] = kitSelectData(expts,varargin)
% KITSELECTDATA(EXPTS,...) Allows the user to choose specific data from a
% given set of EXPTS for use in downstream analysis tools.
%
% Copyright (c) 2017 C. A. Smith
% Modified by C. C. Conway (2022)

opts.addMore = [];
opts.channel = 1;
opts.contrast = [0.1 1];
opts.dataType = 'spots'; % can also be 'sisters'
opts.lineProfile = 0; %whether to plot line profiles in x and y across spots
opts.method = 'deselect'; % can also be 'select'
opts.startMovie = 1;
opts.zProject = -1;
opts.autoSize = 0;
opts = processOptions(opts,varargin{:});

% check for previous selections, set up new structure if none given
if isempty(opts.addMore)
    selectedData.dataType = opts.dataType;
    selectedData.selection = {[]};
    preprocExpts = 0;
else
    selectedData = opts.addMore;
    preprocExpts = length(selectedData.selection);
    % ensure that dataTypes match
    if ~strcmp(selectedData.dataType,opts.dataType)
      error('Cannot add more data to the provided selection due to clash in dataType.');
    end
end

% get a filename and save directory if not already provided
[filename,filepath] = uiputfile('*.mat','Save selection file','selectedSpots.mat');

% check whether input is from one or more experiments, change accordingly
if ~iscell(expts{1})
  expts = {expts};
end
nExpts = length(expts);

% get basic information
c = opts.channel;

% get previous results, if any
allSels = selectedData.selection;

autoResize = opts.autoSize;
if autoResize
    % set up figure
    f1 = figure(1); %CCC addition to presize figure
    szScreens = get(0, 'MonitorPositions');
    if size(szScreens,1) == 1
        pause(0.001);
        f1.WindowState = 'maximized';
        pause(0.001);
    else
        szScreens(2,4) = szScreens(2,4) - 30;
        pause(0.001);
        f1.WindowState = 'normal';
        pause(0.001);
        f1.OuterPosition = szScreens(2,:);
        pause(0.001);
        f1.Resize = 'off';
    end
end

% loop over all movies
for iExpt = 1:nExpts
    
    % get the jobs for this experiment
    jobs = expts{iExpt};
    nJobs = length(jobs);
    
    kitLog('Checking experiment %i of %i...',iExpt,nExpts);
    % start progress bar
    prog = kitProgress(0);
    
    % make empty allSels cell
    allSels{iExpt+preprocExpts} = [];

    % for easy access of spot selections to refer to later - make empty
    % spotSel structure (CCC addition 2022/10/20)
    spotSel = [];

    for iJob = opts.startMovie:nJobs
        % Cell number ID
        kitLog('Spots in movie %i (%i of %i)', jobs{iJob}.index, iJob, nJobs);
        % check whether there is any data contained within this movie
        if ~isfield(jobs{iJob},'dataStruct') || ~isfield(jobs{iJob}.dataStruct{c},'failed') || jobs{iJob}.dataStruct{c}.failed
            continue
        end
        % check whether there are any tracks contained within this movie
        if ~isfield(jobs{iJob}.dataStruct{c},'tracks') || jobs{iJob}.dataStruct{c}.tracks(1).tracksFeatIndxCG == 0
            continue
        end
        % get dataStruct
        dS = jobs{iJob}.dataStruct{c};
        
        switch opts.dataType
          
          case 'sisters'
            % show all sisters
            kitShowAllSisters(jobs{iJob},'channel',c,'contrast',opts.contrast,...
                 'zProject',opts.zProject);
            nData = length(dS.sisterList);
                
          case 'spots'
            % show all spots - defined using tracks
            kitShowAllSpots(jobs{iJob},'channel',c,'contrast',opts.contrast,...
                'zProject',opts.zProject, 'lineProfile',opts.lineProfile); %CCC edit 2022/10/18
            nData = length(dS.trackList);
        end
        
        switch opts.method
          % provide a message to request lists of spots/sisters  
          case 'select'
              kitLog('Please list all %s to be selected: ',opts.dataType);
              tempList = input('');
              % check for any sisters not in the list, provide warning and
              % remove if so
              incorrect = setdiff(tempList,1:nData);
              if ~isempty(incorrect)
                  warning('The following selected %s do not exist: %s. Will ignore.',opts.dataType,num2str(incorrect));
                  tempList = setdiff(tempList,incorrect);
              end
              % invert list to make ignored rather than selected (CCC 2022/10/20)
              ignoredList = setdiff(1:nData, tempList);
              % add to spotSel structure
              spotSelIDs = [string(num2str(jobs{iJob}.index)) string(num2str(tempList)) string(num2str(ignoredList))];
              %Cell, selected spots, ignored spots
              spotSel = cat(1, spotSel, spotSelIDs);

          case 'deselect'
              kitLog('Please list all %s to be ignored: ',opts.dataType);
              tempList = input('');
              % check for any sisters not in the list, provide warning and
              % remove if so
              incorrect = setdiff(tempList,1:nData);
              if ~isempty(incorrect)
                  warning('The following selected %s do not exist: %s. Will ignore.',opts.dataType,num2str(incorrect));
              end
              % invert the list to make it selected rather than ignored
              tempList = setdiff(1:nData,tempList);
              % invert list to make ignored rather than selected (CCC 2022/10/20)
              ignoredList = setdiff(1:nData, tempList);
              % add to spotSel structure
              spotSelIDs = [string(num2str(jobs{iJob}.index)) string(num2str(tempList)) string(num2str(ignoredList))];
              %Cell, selected spots, ignored spots
              spotSel = cat(1, spotSel, spotSelIDs);
        end
        
        % process the list
        if isempty(tempList)
            continue
        else
            nData = length(tempList);
            tempStruct = ones(nData,2)*iJob;
            tempStruct(:,2) = tempList;
            % collate selections
            allSels{iExpt+preprocExpts} = [allSels{iExpt+preprocExpts}; tempStruct];
        end
        
        % update progress
        prog = kitProgress(iJob/nJobs,prog);
        
    end
    
    % store final list in output, and save
    selectedData.selection = allSels;
    save(fullfile(filepath,filename),'selectedData');
    kitLog('Updated save file up to experiment %i of %i.',iExpt,nExpts);

end
if autoResize
    f1 = figure(1);
    f1.Resize = 'on';
end
kitLog('Data selection complete.');

end


    
