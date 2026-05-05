
%
% Place your data entries here. 
%
% ExptDataLoc is location relative to MCMC directory. Can be relative or absolute.
% Give file name and in ExptDats give type of fluorophore
%
% Add more as needed.
%

ExptDataLoc= '../ExptData/NormMeanZero/';
FileNames{1}='intraMeasurements_HeLaTwoColourCENPA_central.mat'; % Dual label
FileNames{2}='intraMeasurements_HeLaCENPANdc80_central.mat';
FileNames{3}='intraMeasurements_RPECENPCNdc809G3_central.mat';
FileNames{4}='intraMeasurements_RPECENPCBub1_central.mat';
FileNames{5}='intraMeasurements_RPECENPC488594_central.mat'; % Duallabel. Untreated only.
FileNames{6}='iM_161122_161206_170130_170213_170309_9G3C_dmso_sel_centr.mat';

% Create an associated ExptDat
ExptDats{1}.name='HeLaTwoColourCENPA'; ExptDats{1}.fluorophore='Ab';
ExptDats{2}.name='HeLaCENPANdc80'; ExptDats{2}.fluorophore='GFP';
ExptDats{3}.name='RPECENPCNdc809G3'; ExptDats{3}.fluorophore='Ab';
ExptDats{4}.name='RPECENPCBub1'; ExptDats{4}.fluorophore='Ab';
ExptDats{5}.name='RPECENPCTwoAb'; ExptDats{5}.fluorophore='Ab';
ExptDats{6}.name='RPECENPC9G3'; ExptDats{6}.fluorophore='Ab';
