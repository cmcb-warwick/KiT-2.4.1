

%
% Code to run MCMC on expts to call EuclDistCorrection_Driver_fn (cf EuclDistCorrection_Driver)
%

% Set renderer
set(0,'DefaultFigureRendererMode','manual')
set(0,'DefaultFigureRenderer','OpenGl')


  
%  
% Load default filters
%
% You can change these filters for specific runs
filterpars.filterdist=250; % max 3D distance allowed in sample data
filterpars.zfilterz=100; % zfilter value. Comparison only, not used in MCMC. 
options=[];


%
% FileNames{1}='intraMeasurements_RPECENPCNdc809G3_centralGood.mat'; % Compilation
%

ExptData_list  % loads the ExptData List
[datsm datsp nams]=plot_data(ExptDataLoc,FileNames{6}); % Plot the data

%options.a=4; % Default
%options.nsteps=25000;  % Default is 10000
%options.nchains=5; % default is 4

%options.mu=5;options.sigmeans=[12 30];
filterpars.filterdist=200; % max 3D distance allowed in sample data

%
% Emanuele USE THIS to run. Increase nsteps (x2 or x4) if convergence poor Rc> 1.1 in any of mu, sigx, sigz
%
  EuclDistCorrection_Driver_fn(6,'DMSO','a=4',filterpars,struct('sigs',[25 45],'nsteps',1600000,'a',4)); % use struct()


  return;
[algp conv]=EuclDistCorrection_Driver_fn(1,'DMSO','',filterpars,options); % use struct()

disp('Runs saved in')
algp.rundir


return
  % Dead
  Dead
  options.a=25;options.mu=5;options.sigmeans=[12 30];
filterpars.filterdist=100; % max 3D distance allowed in sample data

  EuclDistCorrection_Driver_fn(5,'DMSO','',filterpars,options)

%
% FileNames{4}='NormMeanZero/intraMeasurements_RPECENPCBub1_central.mat';
%

options.a=10;
%options.mu=5;options.sigmeans=[12 30];
filterpars.filterdist=200; % max 3D distance allowed in sample data

EuclDistCorrection_Driver_fn(4,'DMSO','',filterpars,options); % Gave mu=0.

options.a=100;
EuclDistCorrection_Driver_fn(4,'noc','',filterpars,options)


%
% FileNames{6}='NormMeanZero/iM_161122_161206_170130_170213_170309_9G3C_dmso_sel_centr.mat' 
%

options.a=10;
%options.mu=5;options.sigmeans=[12 30];
filterpars.filterdist=200; % max 3D distance allowed in sample data

[datsm datsp nams]=plot_data(ExptDataLoc,FileNames{6}); % Plot the data
EuclDistCorrection_Driver_fn(6,'DMSO','',filterpars,options); % Gave mu=0.

