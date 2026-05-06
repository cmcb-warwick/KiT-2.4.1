function selectedData = kitGetPairedSpotSelFromFilteredData(mS_unpairedFiltered, mS_pairedUnfiltered)
% KITGETPAIREDSPOTSELFROMFILTEREDDATA allows the creation of a complete
% spot selection for a paired movieStructure for the coordinate system
% channel using the filtered unpaired movieStructure (generally created 
% when performing intensity measurements only). You will not need to pass
% the selectedData from this function to kitUpdateSpotSelections, unlike
% most other KiT spot selection functions.
%
%   KITGETPAIREDSPOTSELFROMFILTEREDDATA(mS_unpairedFiltered, mS_pairedUnfiltered)
%   mS_unpairedFiltered should be an unpaired movieStructure where
%   coordinate system channel spots have been filtered using kitFilterSpots
%   in the command window or from the kitRun GUI (generally you will have
%   created this to measure intensities in 3 channels). mS_pairedUnfiltered
%   should be a paired movieStructure where no filtering has taken place.
%   You should ensure that both mS structures contain the same number of
%   cells in the same order, otherwise you will get an unusable output
%   and/or an error.
%
%   Edited from kitSelectData_sS C. C. Conway 2023

%automatically put dataType as 'spots', as normal
selectedData.dataType = 'spots';
selectedData.selection = {[]};
selectedData.rawSelection = {[]};

allSels{1} = [];
rawSels{1} = [];
nCells = length(mS_pairedUnfiltered);

%for each cell find out if spots have actually been selected. Skip if not
for iCell = 1:nCells
    if mS_unpairedFiltered{iCell}.dataStruct{2}.failed
        continue
    end
    %otherwise, find out which spots have passed filtering
    selectedIDs = find(~isnan(mS_unpairedFiltered{iCell}.dataStruct{2}.initCoord.allCoord(:,1)));
    nKTs = length(mS_pairedUnfiltered{iCell}.dataStruct{2}.tracks);
    tempList = [];
    rawList = [];
    nIC = length(selectedIDs);
    %
    for iKT = 1:nKTs
        iKTic = mS_pairedUnfiltered{iCell}.dataStruct{2}.tracks(iKT).tracksFeatIndxCG;
        for iIC = 1:nIC
            if iKTic == selectedIDs(iIC)
                tempList = [tempList iKT];
                rawList = [rawList iKTic];
            end
        end
    end
    % process the list
    if isempty(tempList)
        continue
    else
        nData = length(tempList);
        tempStruct = ones(nData,2)*iCell;
        rawStruct = tempStruct;
        tempStruct(:,2) = tempList;
        rawStruct(:,2) = rawList;
        % collate selections
        allSels{1} = [allSels{1}; tempStruct];
        rawSels{1} = [rawSels{1}; rawStruct];
    end
        

        
end
    
    % store final list in output
    selectedData.selection = allSels;
    selectedData.rawSelection = rawSels;

end


