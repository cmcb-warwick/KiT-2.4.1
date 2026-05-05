
%
% Place your data entries here. 
%
% ExptDataLoc is location relative to MCMC directory. Can be relative or absolute.
% Give file name and in ExptDats give type of fluorophore
%
% Add more as needed.
%

ExptDataLoc= '../ExptData/';
FileNames{1}='intraMeasurements_RPECENPCNdc809G3_centralGood.mat'; % Compilation of all CENPC to 9G3 data

% Create an associated ExptDat
ExptDats{1}.name='RPECENPC9G3'; ExptDats{1}.fluorophore='Ab';
