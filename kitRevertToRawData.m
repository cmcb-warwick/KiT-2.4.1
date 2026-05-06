function kitRevertToRawData(jobset)
% KITREVERTTORAWDATA(JOBSET) Restores original JOBSET data from JOBSET
% with filtered spots.
%
% C. C. Conway (2022)


% Get the data.
job = kitLoadAllJobs(jobset);

handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),jobset.options.spotMode));

% Predefine some handles required during looping
for iMov = 1:length(job)

    % get dataStruct
    for jChan = handles.chans
        dS = job{iMov}.dataStruct{jChan};
        if isfield(dS,'rawData')
            iC = dS.rawData.initCoord;
            sI = dS.rawData.spotInt;
            dS.initCoord = iC;
            dS.spotInt = sI;
            job{iMov}.dataStruct{jChan} = dS;
        else
            continue
        end
    end
       

    
    % save results
    job{iMov} = kitSaveJob(job{iMov});

end


% Re-run plane fitting for jobset with new filtered data.
kitLog('Re-fitting planes to filtered data.');
kitRunJob(jobset,'existing',1,'tasks',[2 6]);%kitRunJob
end

