function makeiMinStepsWithCellPhase(nBins, IntiM, iM, cellPhases, normalised, markerdist, intID, expnum)
%2026.01.16, making script so I can see whether cell phase influences
%intraKT delta. IntiM and iM should be produced from 
%kitNormCombineiMsIntiMs, and cellPhases should be provided in
%the same order that the experiments were given to the combining function.
%cellPhases should be provided as cell of tables. See makeiMinSteps for
%more info on the operation of this function!
%
% Copyright (c) 2026 C. C. Conway
    
cellPhaseClassificationData = [];

for iPair = 1:length(IntiM.label)
    currentExpt = str2num(IntiM.label(iPair, 1:2));
    currentCell = str2num(IntiM.label(iPair, 3:4));

    if strcmpi('Rosette', cellPhases{currentExpt}.Phase{currentCell})
        cellPhaseClassificationData = cat(1, cellPhaseClassificationData, [1 1]);

    elseif strcmpi('Congressing', cellPhases{currentExpt}.Phase{currentCell})
        cellPhaseClassificationData = cat(1, cellPhaseClassificationData, [2 2]);

    elseif strcmpi('Late prometa', cellPhases{currentExpt}.Phase{currentCell})
        cellPhaseClassificationData = cat(1, cellPhaseClassificationData, [3 3]);

    elseif strcmpi('Metaphase', cellPhases{currentExpt}.Phase{currentCell})
        cellPhaseClassificationData = cat(1, cellPhaseClassificationData, [4 4]);

    end

end




Total = size(find(~isnan(IntiM.intensity.mean.inner(:))),1);
stepsize = floor(Total/nBins);
bin = 1;
n = stepsize;
k = n;
if strcmpi(normalised,'norm') || strcmpi(normalised,'normalised')
    nStr = 'N';
    nv = IntiM.intensity.mean.outer./IntiM.intensity.mean.inner;
else
    nStr = '';
    nv = IntiM.intensity.mean.outer;
end

while n < Total+1
    cellClassCopy = cellPhaseClassificationData;
    iM_copy = iM;
    Subset = mink(maxk(nv(:), n), k); %20250430
    
    %no nan removal needed here bc logical comparison defaults comparisons
    %between nan and other values to false
    positionLow = nv<min(Subset);
    positionLowSisSep = max(positionLow, [], 2);
    positionHigh = nv>max(Subset);
    positionHighSisSep = max(positionHigh, [], 2);
    
    iM_copy.microscope.raw.delta.x.all(positionLow) = NaN;
    iM_copy.microscope.raw.delta.y.all(positionLow) = NaN;
    iM_copy.microscope.raw.delta.z.all(positionLow) = NaN;
    cellClassCopy(positionLow) = NaN;
    if isfield(iM_copy.microscope, 'sisSep')
        iM_copy.microscope.sisSep.threeD(positionLowSisSep) = NaN;
        iM_copy.microscope.raw.swivel.threeD.all(positionLow) = NaN;
    end
    
    iM_copy.microscope.raw.delta.x.all(positionHigh) = NaN;
    iM_copy.microscope.raw.delta.y.all(positionHigh) = NaN;
    iM_copy.microscope.raw.delta.z.all(positionHigh) = NaN;
    cellClassCopy(positionHigh) = NaN;
    if isfield(iM_copy.microscope, 'sisSep')
        iM_copy.microscope.sisSep.threeD(positionHighSisSep) = NaN;
        iM_copy.microscope.raw.swivel.threeD.all(positionHigh) = NaN;
    end
    %EPM = early prometa; IPM = intermediate prometa (ie mid prometa but
    %wanted things to be alphabetical when saving for my sanity); LPM =
    %late prometa; MMM = metaphase.
    phaseStrings = {'EPM', 'IPM', 'LPM', 'MMM'};
    for iPhase = 1:4
        iMSecondCopy = iM_copy;
        removePosition = cellClassCopy ~= iPhase;
        removeSisSep = max(removePosition, [], 2);

        iMSecondCopy.microscope.raw.delta.x.all(removePosition) = NaN;
        iMSecondCopy.microscope.raw.delta.y.all(removePosition) = NaN;
        iMSecondCopy.microscope.raw.delta.z.all(removePosition) = NaN;
        if isfield(iMSecondCopy.microscope, 'sisSep')
            iMSecondCopy.microscope.sisSep.threeD(removeSisSep) = NaN;
            iMSecondCopy.microscope.raw.swivel.threeD.all(removePosition) = NaN;
        end
        
        fieldname = sprintf('iM_exp%02d_dist%s_%sint%s_ss%d_%sBin%02d', expnum, markerdist, nStr, intID, stepsize, phaseStrings{iPhase}, bin);
        workingStruct = struct(fieldname, iMSecondCopy);
        f = fieldnames(workingStruct);

        if iPhase == 1
            for i = 1:length(f)
                masterEPMStruct.(f{i}) = workingStruct.(f{i});
            end

        elseif iPhase == 2
            for i = 1:length(f)
                masterIPMStruct.(f{i}) = workingStruct.(f{i});
            end

        elseif iPhase == 3
            for i = 1:length(f)
                masterLPMStruct.(f{i}) = workingStruct.(f{i});
            end

        elseif iPhase == 4
            for i = 1:length(f)
                masterMMMStruct.(f{i}) = workingStruct.(f{i});
            end

        end
    end
   

    n = n+k;
    bin = bin+1;

end

for iPhase = 1:4
    fieldname = sprintf('iM_exp%02d_dist%s_%sint%s_ss%d_%s.mat', expnum, markerdist, nStr, intID, stepsize, phaseStrings{iPhase});
    if iPhase == 1
        save(fieldname, '-struct', 'masterEPMStruct');
    elseif iPhase == 2
        save(fieldname, '-struct', 'masterIPMStruct');
    elseif iPhase == 3
        save(fieldname, '-struct', 'masterLPMStruct');
    elseif iPhase == 4
        save(fieldname, '-struct', 'masterMMMStruct');
    end
end


end