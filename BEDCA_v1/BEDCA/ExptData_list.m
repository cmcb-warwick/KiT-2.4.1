
%
% Place your data entries here. 
%
% ExptDataLoc is location relative to MCMC directory. Can be relative or absolute.
% Give file name and in ExptDats give type of fluorophore
%
% Add more as needed.
%

% 
% This is a treatment name list, first is name used in MCMC files, 2nd entry is search string in your Experiment file names - your KiT file names
% must include it and uniquely identify the files (it is vital not to use search names that are contained within other search terms). 
% Ensure unique. Code does not distinguish between capitals v small letters.
% Note DMSO == untreated so appears twice as two different search terms.
  % Later entries are aka; can be used as runtreatment in EuclDistCorrection_Driver_fn
%
% We used standard concentrations (Noc, Taxol) so allowed names that didnt specify conc.

% Add new treatments to array
TreatmentLib = {{'DMSO','DMSO'},{'Bin01','Bin01'},{'Bin02','Bin02'},...
    {'Bin03','Bin03'},{'Bin04','Bin04'},{'Bin05','Bin05'},{'Bin06','Bin06'},...
    {'Bin07','Bin07'},{'Bin08','Bin08'},{'Bin09','Bin09'},{'Bin10','Bin10'},...
    {'Bin11','Bin11'},{'Bin12','Bin12'},{'Bin13','Bin13'},{'Bin14','Bin14'},...
    {'Bin15','Bin15'},{'Bin16','Bin16'},{'Bin17','Bin17'},{'Bin18','Bin18'},...
    {'Bin19','Bin19'},{'Bin20','Bin20'},{'Bin21','Bin21'},{'Bin22','Bin22'},...
    {'Bin23','Bin23'},{'Bin24','Bin24'},{'Bin25','Bin25'},{'DMSO','untreated'},...
    {'Taxol15min','tax15min','tax'},{'Nocod2hr','noc2hr','noc'},{'Nocod45min','noc45min'},...
    {'Nocod30min','noc30min'},{'Nocod15min','noc15min'},{'DMSOMG','DMSOMG'},...
    {'3uMnocMG','3uMnocMG'},{'3uMnocRevMG','3uMnocRevMG'},{'1uMTaxolMG','1uMtaxMG','taxMG'},...
    {'1uMTaxolRevMG','1uMtaxRevMG','taxRevMG'},{'CTRrnai','CTRi'},{'CenpCrnai','CenpCi'},...
    {'CenpTrnai','CenpTi'},{'NocOnly','NocOnly'},{'ZMOnly','ZMOnly'},{'NocZM','NocZM'}, ...
    {'lowKKOnly','lowKKOnly'}, {'highInt_KKlow','highInt_KKlow'}, {'lowInt_KKlow','lowInt_KKlow'}};
% TreatmentLib = {{'DMSO','DMSO'},{'Bin01_17','Bin01_17'},{'Bin18_25','Bin18_25'},...
%     {'DMSO','untreated'},...
%     {'Taxol15min','tax15min','tax'},{'Nocod2hr','noc2hr','noc'},{'Nocod45min','noc45min'},...
%     {'Nocod30min','noc30min'},{'Nocod15min','noc15min'},{'DMSOMG','DMSOMG'},...
%     {'3uMnocMG','3uMnocMG'},{'3uMnocRevMG','3uMnocRevMG'},{'1uMTaxolMG','1uMtaxMG','taxMG'},...
%     {'1uMTaxolRevMG','1uMtaxRevMG','taxRevMG'},{'CTRrnai','CTRi'},{'CenpCrnai','CenpCi'},...
%     {'CenpTrnai','CenpTi'},{'NocOnly','NocOnly'},{'ZMOnly','ZMOnly'},{'NocZM','NocZM'}};
% Directory with data (relative location to this file)
ExptDataLoc = 'ExptData/';

% FILE NAMES IN ABOVE DIRECTORY
% Compilation sets over days (in ExptLabel)
% the exp26/exp28/exp30 pool and exp27/29/30 pool are DMSO, exp28/exp30 pool are taxol)
FileNames{1} = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSO.mat'; %DMSO,example


% Create an associated ExptDat. Use concatenated markers name.

ExptDats{1}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin01.mat';
ExptDats{2}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin02.mat';
ExptDats{3}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin03.mat';
ExptDats{4}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin04.mat';
ExptDats{5}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin05.mat';
ExptDats{6}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin06.mat';
ExptDats{7}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin07.mat';
ExptDats{8}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin08.mat';
ExptDats{9}.name   = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin09.mat';
ExptDats{10}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin10.mat';
ExptDats{11}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin11.mat';
ExptDats{12}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin12.mat';
ExptDats{13}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin13.mat';
ExptDats{14}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin14.mat';
ExptDats{15}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin15.mat';
ExptDats{16}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin16.mat';
ExptDats{17}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin17.mat';
ExptDats{18}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin18.mat';
ExptDats{19}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin19.mat';
ExptDats{20}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin20.mat';
ExptDats{21}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin21.mat';
ExptDats{22}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin22.mat';
ExptDats{23}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin23.mat';
ExptDats{24}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin24.mat';
ExptDats{25}.name  = 'iM_exp03_distCenpCNdc80_NintMad2_ss168_DMSOBin25.mat';

% This is groupings by paired markers. .name should match across group
%ExptDats{1}.group=[1 2];ExptDats{2}.group=[1 2];  % Bub1

%
% Treatments

ExptDats{1}.treatment   = 'Bin01';
ExptDats{2}.treatment   = 'Bin02';
ExptDats{3}.treatment   = 'Bin03';
ExptDats{4}.treatment   = 'Bin04';
ExptDats{5}.treatment   = 'Bin05';
ExptDats{6}.treatment   = 'Bin06';
ExptDats{7}.treatment   = 'Bin07';
ExptDats{8}.treatment   = 'Bin08';
ExptDats{9}.treatment   = 'Bin09';
ExptDats{10}.treatment  = 'Bin10';
ExptDats{11}.treatment  = 'Bin11';
ExptDats{12}.treatment  = 'Bin12';
ExptDats{13}.treatment  = 'Bin13';
ExptDats{14}.treatment  = 'Bin14';
ExptDats{15}.treatment  = 'Bin15';
ExptDats{16}.treatment  = 'Bin16';
ExptDats{17}.treatment  = 'Bin17';
ExptDats{18}.treatment  = 'Bin18';
ExptDats{19}.treatment  = 'Bin19';
ExptDats{20}.treatment  = 'Bin20';
ExptDats{21}.treatment  = 'Bin21';
ExptDats{22}.treatment  = 'Bin22';
ExptDats{23}.treatment  = 'Bin23';
ExptDats{24}.treatment  = 'Bin24';
ExptDats{25}.treatment  = 'Bin25';

%Fluorophore

ExptDats{1}.fluorophore   = 'Ab';
ExptDats{2}.fluorophore   = 'Ab';
ExptDats{3}.fluorophore   = 'Ab';
ExptDats{4}.fluorophore   = 'Ab';
ExptDats{5}.fluorophore   = 'Ab';
ExptDats{6}.fluorophore   = 'Ab';
ExptDats{7}.fluorophore   = 'Ab';
ExptDats{8}.fluorophore   = 'Ab';
ExptDats{9}.fluorophore   = 'Ab';
ExptDats{10}.fluorophore  = 'Ab';
ExptDats{11}.fluorophore  = 'Ab';
ExptDats{12}.fluorophore  = 'Ab';
ExptDats{13}.fluorophore  = 'Ab';
ExptDats{14}.fluorophore  = 'Ab';
ExptDats{15}.fluorophore  = 'Ab';
ExptDats{16}.fluorophore  = 'Ab';
ExptDats{17}.fluorophore  = 'Ab';
ExptDats{18}.fluorophore  = 'Ab';
ExptDats{19}.fluorophore  = 'Ab';
ExptDats{20}.fluorophore  = 'Ab';
ExptDats{21}.fluorophore  = 'Ab';
ExptDats{22}.fluorophore  = 'Ab';
ExptDats{23}.fluorophore  = 'Ab';
ExptDats{24}.fluorophore  = 'Ab';
ExptDats{25}.fluorophore  = 'Ab';




