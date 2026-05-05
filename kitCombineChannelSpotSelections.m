function sS_comb = kitCombineChannelSpotSelections(sS1, sS2)
% KITCOMBINECHANNELSPOTSELECTIONS allows the user to combine two spot
% selections made on the same jobset in two different channels to a single
% spot selection where only the spots selected in both channels are
% featured. Both spot selections must be either a initial spot selection or
% a final spot selection (kitUpdateSpotSelections either has not or has
% been applied to both datasets, respectively).
%
% Please note that this has only been tested on data selections made with
% 'spots' dataType and only from single experiments with multiple cells.

% Copyright (c) 2024 C. C. Conway

final_sS1 = isfield(sS1, 'rawSelection');
final_sS2 = isfield(sS2, 'rawSelection'); %checking if both sS structures are initial or final
dT_sS1 = sS1.dataType;
dT_sS2 = sS2.dataType; %checking both sS structures share the same dataType
nExpts_sS1 = length(sS1.selection);
nExpts_sS2 = length(sS2.selection); %checking both sS structures have the same number of experiments
if ~strcmp(dT_sS1, dT_sS2)
    error('Provided sS structures have different dataType values. Please review.')
elseif nExpts_sS1 ~= nExpts_sS2
    error('Provided sS structures have different number of experiments. Please review.')
else
    warning off backtrace
    if final_sS1 ~= final_sS2
        warning('One provided sS structure is initial selection and th other is final selection (created with kitUpdateSpotSelections). The output of this function will need to be run through kitUpdateSpotSelections to use downstream.');
    elseif final_sS1
        warning('Both provided sS structures are final selection. You will not need to run the output of this function through kitUpdateSpotSelections.')
    else
        warning('Both provided sS structures are initial selection. You will need to run the output of this function through kitUpdateSpotSelections to use downstream.')
    end

    sS_comb = struct; %make new structure for combining
    sS_comb.dataType = dT_sS1; %reading in data type
    sS_comb.selection = cell(1, nExpts_sS1); %making new cell experiment
    for iSel = 1:nExpts_sS1
        sS_sel = [];
        sS1_sel = sS1.selection{iSel};
        sS2_sel = sS2.selection{iSel};
        for iKT1 = 1:length(sS1_sel)
            wkCell1 = sS1_sel(iKT1, 1);
            wkKT1 = sS1_sel(iKT1, 2);
            for iKT2 = 1:length(sS2_sel)
                wkCell2 = sS2_sel(iKT2, 1);
                wkKT2 = sS2_sel(iKT2, 2);
                if wkCell1 == wkCell2
                    if wkKT1 == wkKT2
                        bothCellKT = [wkCell1 wkKT1];
                        sS_sel = cat(1, sS_sel, bothCellKT);
                    end
                end
            end
        end
        sS_comb.selection{iSel} = sS_sel;
    end

    %check if both datasets are 'final' (created with kitUpdateSpotSelections)
    %if so, make new rawSelection that will not need to be combined to make final
    if (final_sS1 + final_sS2) == 2
        sS_comb.rawSelection = cell(1, nExpts_sS1);
        for iSel = 1:nExpts_sS1
            sS_rawSel = [];
            sS1_rawSel = sS1.rawSelection{iSel};
            sS2_rawSel = sS2.rawSelection{iSel};
            for iKT1 = 1:length(sS1_rawSel)
                wkCell1 = sS1_rawSel(iKT1, 1);
                wkKT1 = sS1_rawSel(iKT1, 2);
                for iKT2 = 1:length(sS2_rawSel)
                    wkCell2 = sS2_rawSel(iKT2, 1);
                    wkKT2 = sS2_rawSel(iKT2, 2);
                    if wkCell1 == wkCell2
                        if wkKT1 == wkKT2
                            bothCellKT = [wkCell1 wkKT1];
                            sS_rawSel = cat(1, sS_rawSel, bothCellKT);
                        end
                    end
                end
            end
            sS_comb.rawSelection{iSel} = sS_rawSel;
        end
    end
end

end