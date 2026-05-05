

%
% List of simulated data with metadata
%
%

SimDataLoc= '../SimData/';
%FileNames{1}='intraMeasurements_RPECENPCNdc809G3_centralGood.mat'; % Compilation of all CENPC to 9G3 data

% Compilation sets over days (in ExptLabel)
FileNames{1}='simInflation_sample1000_v6_lownoise.mat';
FileNames{2}='simInflation_sample1000_v6_normal.mat';
FileNames{3}='simInflation_sample1000_v6_muvar.mat';
FileNames{4}='simInflation_sample1000_v6_smallmuvar';

% 10 runs with variable mu, sample size 1000
for k=1:10
	FileNames{4+k}=['simInflation_sample1000_v6_muvar_' num2str(k) '.mat'];
ExptDats{4+k}.name='muvariation';ExptDats{4+k}.treatment='';
end

% 10 runs with variable mu, sample size 5000
for k=1:10
	FileNames{14+k}=['simInflation_sample5000_v6_muvar_' num2str(k) '.mat'];
ExptDats{14+k}.name='muvariation';ExptDats{14+k}.treatment='';
end

FileNames{25}='simInflation_N1000_mu40_25_45';
ExptDats{25}.name='';
FileNames{26}='simInflation_N1000_mu60_25_45';
ExptDats{26}.name='';
FileNames{27}='simInflation_N1000_mu40_25_45_Set2';
ExptDats{27}.name='';

FileNames{28}='PE-struct_12_Trianglestate 2019-8-13 (l=[1.400e+03,6.0e+01,5.0e+01,2.0e+01],merr=[2.5e+01,4.0e+_etc.mat';
FileNames{29}='simInflation__PeterEmTestValues.mat'

% Create an associated ExptDat. Marker concatenated name.
ExptDats{1}.name='lownoise';
ExptDats{2}.name='normal';
ExptDats{3}.name='muvariation';
ExptDats{4}.name='smmuvariation';
ExptDats{28}.name='PeterE_Test12';
ExptDats{29}.name='SimulatedPeterE_Test12';

% This is groupings by paired markers. .name should match across group
ExptDats{1}.group=[1 2];ExptDats{2}.group=[1 2];  % Bub1
ExptDats{3}.group=[3 4];ExptDats{4}.group=[3 4];

ExptDats{1}.treatment='';ExptDats{3}.treatment='';
ExptDats{2}.treatment='';ExptDats{4}.treatment='';




