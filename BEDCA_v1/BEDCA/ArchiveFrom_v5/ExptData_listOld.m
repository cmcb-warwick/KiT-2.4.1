
%
% Place your data entries here. 
%
% ExptDataLoc is location relative to MCMC directory. Can be relative or absolute.
% Give file name and in ExptDats give type of fluorophore
%
% Add more as needed.
%

% 
% this is treatment name list, first is name used in MCMC files, 2nd entry is search string in your Experiment file names so your name must include it in your KiT file names. Any others can be used
   % in treatment specification in analysis code.
   % Ensure unique. Do not distinguish with capitals v small letters.
   % Note DMSO == untreated so appears twice
   % We used standard concentrations (Noc, Taxol) so allowed names that didnt specify conc.
  TreatmentLib={{'DMSO','DMSO'},{'DMSO','untreated'},{'Taxol15min','tax15min','tax'},{'Nocod2hr','noc2hr','noc'},{'Nocod45min','noc45min'},{'Nocod30min','noc30min'},{'Nocod15min','noc15min'},{'DMSOMG','DMSOMG'},{'3uMnocMG','3uMnocMG'},{'3uMnocRevMG','3uMnocRevMG'},{'1uMTaxolMG','1uMtaxMG','taxMG'},{'1uMTaxolRevMG','1uMtaxRevMG','taxRevMG'},{'CTRrnai','CTRi'},{'CenpCrnai','CenpCi'},{'CenpTrnai','CenpTi'}};

ExptDataLoc= '../ExptData_Roscioli/';
%FileNames{1}='intraMeasurements_RPECENPCNdc809G3_centralGood.mat'; % Compilation of all CENPC to 9G3 data

% Compilation sets over days (in ExptLabel)
FileNames{1}='iM_05_sel_07_08_13_18_Bub1C_DMSO_selectedSpotsBub1.mat';
FileNames{2}='iM_07_08_13_18_Bub1C_3uMnoc2hr_spotsExpanded.mat';
FileNames{3}='iM_161122_161206_170130_170213_170925_9G3C_DMSO_sel_centr.mat';
FileNames{4}='iM_161122_170929_9G3C_3uMnoc2hr_sel_centr.mat';
FileNames{5}='iMs_RPECENPCBub1(1_130_Green)_Noc_nearAllFullDepol_exp0507081318.mat';
FileNames{6}='iM_180903_Ndc80GFPC_DMSO_selectedSpots9G3_Ndc80GFP.mat';
FileNames{7}='iM_180903_Ndc80GFPC_3uMnoc2hr_spots9G3_Ndc80GFP.mat';
FileNames{8}='iM_180903_Ndc80GFPC_1uMtax15min_spots9G3_Ndc80GFP.mat';
FileNames{9}='iM_190423_e6_aGSka1C_2115upd.mat';
FileNames{10}='exp11Lina_iM_05032019_9G3CENPC_HeLa_NOpreExtract.mat';

FileNames{11}='iM_190930_C9G3_DMSOMG_9G3Ndc80GFPCenpCspots.mat';

FileNames{12}='iM_200307_CMad1N_3uMnoc2hr_i_spotspMad1Mad1N.mat';
FileNames{13}='iM_200307_CMad1N_3uMnoc2hr_ii_spotspMad1Mad1N.mat';
FileNames{14}='iM_200307_CMad1N_3uMnoc2hr_i_ii_spotspMad1Mad1N.mat';
FileNames{15}='iM_20030212_pMad1Mad1N_3uMnoc2hr_ii_spotspMad1Mad1N.mat';
FileNames{16}='exptLGiM_28062020exp28_NGSka19G3CenpC_DMSO_spotsSkaCenpC_9G3CenpC.mat';
FileNames{17}='exptLGiM_28062020exp28_NGSka19G3CenpC_DMSO_spotsSkaCenpC_Ska9G3.mat';
FileNames{18}='exptLGiM_28062020exp28_NGSka19G3CenpC_DMSO_spotsSkaCenpC_SkaCenpC.mat';


% Create an associated ExptDat. Marker concatenated name.
ExptDats{1}.name='RPECENPCBUB1';
ExptDats{2}.name='RPECENPCBUB1';
ExptDats{3}.name='RPECENPC9G3';
ExptDats{4}.name='RPECENPC9G3';
ExptDats{5}.name='RPECENPCBUB1';
ExptDats{6}.name='RPECENPCNdc80GFP';
ExptDats{7}.name='RPECENPCNdc80GFP';
ExptDats{8}.name='RPECENPCNdc80GFP';
ExptDats{9}.name='RPELinaSka';
ExptDats{10}.name='HeLaLinaSka';
ExptDats{11}.name='HeLaCENPC9G3';
ExptDats{12}.name='RPECENPCMad1N';
ExptDats{13}.name='RPECENPCMad1N';
ExptDats{14}.name='RPECENPCMad1N';
ExptDats{15}.name='RPEpMad1Mad1N';
ExptDats{16}.name='HeLaLina9G3CenpC';
ExptDats{17}.name='HeLaLinaSka9G3';
ExptDats{18}.name='HeLaLinaSkaCenpC';

% This is groupings by paired markers. .name should match across group
ExptDats{1}.group=[1 2];ExptDats{2}.group=[1 2];  % Bub1
ExptDats{3}.group=[3 4];ExptDats{4}.group=[3 4];  % 9G3C
ExptDats{6}.group=[6 7 8];
ExptDats{7}.group=[6 7 8];
ExptDats{8}.group=[6 7 8];
ExptDats{9}.group=[9];

ExptDats{1}.treatment='DMSO';ExptDats{3}.treatment='DMSO';
ExptDats{2}.treatment='noc2hr';ExptDats{4}.treatment='noc2hr';
ExptDats{5}.treatment='noc2hr';
ExptDats{6}.treatment='DMSO';
ExptDats{7}.treatment='noc2hr';
ExptDats{8}.treatment='tax15min';
ExptDats{11}.treatment='DMSOMG';


ExptDats{1}.fluorophore='Ab';
ExptDats{2}.fluorophore='Ab';
ExptDats{3}.fluorophore='Ab';
ExptDats{4}.fluorophore='Ab';
ExptDats{5}.fluorophore='Ab';
ExptDats{6}.fluorophore='GFP';
ExptDats{7}.fluorophore='GFP';
ExptDats{8}.fluorophore='GFP';
ExptDats{9}.fluorophore='GFP';
ExptDats{10}.fluorophore='Ab';
ExptDats{11}.fluorophore='Ab';
ExptDats{12}.fluorophore='Ab';
ExptDats{13}.fluorophore='Ab';
ExptDats{14}.fluorophore='Ab';
ExptDats{15}.fluorophore='Ab';
ExptDats{16}.fluorophore='Ab';
