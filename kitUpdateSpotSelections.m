function sS_new = kitUpdateSpotSelections(sS_old,mS_full,cellSel,force)
% Use this to update spot selections originally made using only a selection
% of cells, e.g. sS_..._sel_... = kitSelectData(mS_..._sel_...);
% Give it the FULL mS, as all information is required.

if nargin<4
    force = 0;
end
if nargin<3
    fprintf('Please provide a selection of cells. If the original spot selection\n');
    fprintf('was generated using the full movieStructure, type the following line:\n');
    fprintf('   sS_..._new = emanuele_updateSpotSels(sS_...,mS_...,[])\n');
    fprintf('i.e. with the [] provided where cellSelection would be.');
end

if isfield(sS_old,'rawSelection') && ~force
    fprintf('\nWARNING:\n')
    fprintf('This spot selection structure appears to be up to date.\n')
    fprintf('If you are certain that it was made using a selection of mS_...\n')
    fprintf('rather than the full dataset, then please type the following line:\n')
    fprintf('   sS_..._new = emanuele_updateSpotSels(sS_...,mS_...,cellSelection,1)\n')
    fprintf('i.e. with the ,1) at the end.\n\n')
    sS_new = sS_old;
    return
end

nExpts = length(sS_old.selection);
if nExpts > 1
    error('Please run this for only one experiment at a time. Combine them later.');
end
if isempty(cellSel)
    cellSel = 1:length(mS_full);
end

if max(sS_old.selection{1}(:,1))>max(cellSel)
    error('Maximum value in spot selection larger than largest value in cell selection. Not possible.');
end
if length(mS_full) < max(cellSel)
    error('Number of movies in movieStructure smaller than the larger value in cell selection. Ensure you provide the full movieStructure.');
end

if force
    fprintf('Forcing changes.');
end

chan = mS_full{1}.options.coordSystemChannel;
reps = 1+strcmp(sS_old.dataType,'sisters');
rawSelection = repmat(sS_old.selection{1},reps,1);
rawSelection(:) = NaN;

sS_new = sS_old;
nKTs = size(sS_old.selection{1},1);
for iKT = 1:nKTs
    iMov = sS_new.selection{1}(iKT,1);
    sS_new.selection{1}(iKT,1) = cellSel(iMov);
    
    iMov = cellSel(iMov);
    dS = mS_full{iMov}.dataStruct{chan};
    switch sS_new.dataType
        case 'spots'
            tids = sS_new.selection{1}(iKT,2);
            sids = dS.trackList(tids).featIndx;
            rawSelection(iKT,:) = sS_new.selection{1}(iKT,:);
            rawSelection(iKT,2) = sids;
        case 'sisters'
            tids = dS.sisterList(1).trackPairs(iKT,:);
            sids = [dS.trackList(tids(1)).featIndx dS.trackList(tids(2)).featIndx];
            rawSelection(2*iKT-1:2*iKT,:) = repmat(sS_new.selection{1}(iKT,:),2,1);
            rawSelection(2*iKT-1:2*iKT,2) = sids;
    end
        
end

sS_new.rawSelection{1} = rawSelection;

end








