
%
% driver for mcmc_ed_main.
% Loads data and sets up run.
%
% Distances measured in nm
%
% NJB May 2017

initialisetype='priors'; %'true'; % 'priors';
expansionfactor=1;    % For prior
datset='ExptFormat2'; %'Expt'; % 'Sim'; % ExptFormat2
algstr='ed';

% Set path
path(path,genpath(pwd))
path(path,'..')  % For simulate_3D
  
% Load data. dat is matrix (r,th,phi,indexOnZthreshold)
switch datset
case 'Sim'
FileName='simInflation_sample1000_v5';
load(['../SimData/' FileName '.mat']);
disp(['Data file = ' FileName]); 
SAVDIR=[]; %['../MCMCruns/Sim_' FileName];

%Filter
%dat=dat(find(dat(:,2)>pi/4 & dat(:,2)<3*pi/4),:); % filter on theta

unts={'nm','rad','rad'}; % True?

case 'Expt'
% This data has a biased distribution in theta. Why?

FileName='InitialDataSet/iM_161122_161206_170130_170213_170309_9G3C_dmso_sel';
A=load(['../ExptData/' FileName]);
dat=extract_sphericalcoords(A.iM_161122_161206_170130_170213_170309_9G3C_dmso_sel.plate.raw.delta);  % r, theta

% turn to nm
dat(:,1)=1000*dat(:,1);

% Filter out large distances.
J=find(dat(:,1)<=300);
dat=dat(J,:);

% Reflection randomisation in z-axis.
%dat(:,2)=dat(:,2)+(pi-2*dat(:,2)).*(rand(size(dat,1),1)>0.5);
%disp('Performed reflection symmetry in z plane randomisation');

unts={'nm','rad','rad'};
params.mu=80;
params.sig=expansionfactor*[25 75]; % Rough values

FileNameSh='iM_161122_161206_170130_170213_170309_9G3C_dmso_sel';

SAVDIR=['../MCMCruns/Expt_' FileNameSh];

case 'ExptFormat2'
% Has a top structure

% This data has a biased distribution in theta. Why?

%FileName='NormMeanZero/intraMeasurements_HeLaTwoColourCENPA_central.mat';
FileName='NormMeanZero/intraMeasurements_HeLaCENPANdc80_central.mat';
%FileName='NormMeanZero/intraMeasurements_RPECENPCNdc809G3_central.mat';


disp(['Loading ' FileName])

%FileName='intraMeasurements_RPECENPCNdc809G3.mat'; %'InitialDataSet/iM_161122_161206_170130_170213_170309_9G3C_dmso_sel';
A=load(['../ExptData/' FileName]);
%dat=extract_sphericalcoords(A.iM_161122_161206_170130_170213_170309_9G3C_dmso_sel.plate.raw.delta);  % r, theta
%dat=extract_sphericalcoords(A.iM_RCENPC9G3_3uMnoc2hr.plate.raw.delta);  % r, theta
  
%  dat=extract_sphericalcoords(A.iM_160714_HHCENPA_untreated_pairedSelectedCentralised.plate.raw.delta);  % r, theta
dat=extract_sphericalcoords(A.iM_160714_HeCENPAttNdc80C_untreated_pairedSelectedCentralised.plate.raw.delta);  % r, theta

%dat=extract_sphericalcoords(A.iM_RCENPC9G3_DMSO_central.plate.raw.delta);
%dat=extract_sphericalcoords(A.iM_RCENPC9G3_3uMnoc2hr_central.plate.raw.delta);
%dat=extract_sphericalcoords(A.iM_RCENPC9G3_1uMtax15min_central.plate.raw.delta);

% turn to nm
dat(:,1)=1000*dat(:,1);

% Filter out large distances.
J=find(dat(:,1)<=300);
dat=dat(J,:);

% Reflection randomisation in z-axis.
%dat(:,2)=dat(:,2)+(pi-2*dat(:,2)).*(rand(size(dat,1),1)>0.5);
%disp('Performed reflection symmetry in z plane randomisation');

unts={'nm','rad','rad'};
params.mu=80;
params.sig=expansionfactor*[25 75]; % Rough values

J=strfind(FileName,'/');
if ~isempty(J) & J(end)==length(FileName)
  J(end)=[];
end

if ~isempty(J)
FileNameSh=FileName((J(end)+1):end);
 else
   FileNameSh=FileName;
end
FileNmaeSh=[FileNameSh '_' algstr];

SAVDIR=['../MCMCruns/Expt_' FileNameSh];

end

if ~isempty(SAVDIR)
mkdir(SAVDIR);
end

% data processing
% ensure thetaX>0
J=find(dat(:,2)<0);
if ~isempty(J)
disp(['Warning: ' num2str(length(J)) ' data points have negative theta. Taking abs value.']);
dat(J,2)=-dat(J,2);
end

% MCMC Run length. marg (with twist) needs 200000. modsel 500000
algparams.n=40000;
algparams.subsample=1;
algparams.burnin=10000;
algparams.bloks=5;

%
% Priors
%
% Gamma G(a,b) has mean a/b, var a/b^2. Rel error is 1/sqrt(a).  
% sig has mean sqrt(b/a), and variance (a/b^2)/(4 a^3/b^3)=b/(4a^2)
  % Prior tau_x Need b/a = 25.5^2=650, for rel error 10% -> a=100, b=65000
%
 
algparams.priors.variable={'mu','taux','tauz'};
algparams.priors.types={'Gaussian','Gamma','Gamma'}; % Gaussian mean variance,...  Gamma power decaycoef not decaylength as used in matlab
a=1;
algparams.priors.params={[70 10000],[a params.sig(1)^2*a],[a params.sig(2)^2*a]}; % 10% error priors 
%algparams.priors.params={[70 10000000],[1. 1/100],[1 1/100]}; % Weak priors for Gamma are [1 0.01]

%algparams.priors.params={[80 10000000],[a 1/(params.sig(1)*gamma(a)/gamma(a-0.5))^2],[1 1/(params.sig(2)*gamma(a)/gamma(a-0.5))^2]}; % Strong priors

disp(['Priors on sdx, sdz are: rel error  ' num2str(100/sqrt(a)) '\% with means ' num2str(params.sig(1)) ', ' num2str(params.sig(2))]);  

if strcmp(datset,'Sim')
dphi=dat(:,6)-dat(:,3);
J=find(dphi>pi);if ~isempty(J) dphi(J)=dphi(J)-2*pi;end
J=find(dphi<-pi);if ~isempty(J) dphi(J)=dphi(J)+2*pi;end
 else
   dphi=[];
end

% Model priors
%algparams.modpriors=[0.5 0.5]; % mu prior (modsel alg) for mu=zero


switch initialisetype
case 'true'
% Load initial conditions.
algparams.init=[params.mu params.sig];   % True values
algparams.initthetas=dat(:,5);           % True thetas vec n
algparams.initphis=zeros(size(dat,1),1); % True values vec n. Angle relative to true vector

case 'priors'  % Draw from priors.

% mu. Positive.
mu=prior_draw(algparams.priors,'mu');
while (mu < 0)
   mu=prior_draw(algparams.priors,'mu');
end

algparams.init=[mu 1/sqrt(prior_draw(algparams.priors,'taux')) 1/sqrt(prior_draw(algparams.priors,'tauz'))];
algparams.initthetas=pi*rand(size(dat,1),1);
algparams.initphis=normrnd(0,params.sig(1)/(100*params.mu),size(dat,1),1); % Based on expected values

%algparams.initphis=2*pi*rand(size(dat,1),1)-pi; % This causes convergence problems. too wide
algparams.initphis=algparams.initphis-pi*floor(algparams.initphis/pi);

 
otherwise
 disp(['No such initialidsation type as ' initialisetype]);
 disp('ABORTING.');
return

end

algparams.initialise=initialisetype;

plot_inflatedat(dat,[]);

% Overriding. WARNING
%algparams.init(1)=1000;
%algparams.initthetas=dat(:,5);
%algparams.initphis=zeros(size(dat,1),1); %
%algparams.init(1)=0; % set mu=0


disp(['Data set ' FileName]);

% MCMC algorithm. Phi is difference
[mcmcparams, mcmcrun, mcmcruntheta, mcmcrunphi]=mcmc_ed_main3Dmarg([dat(:,[1 2]) dphi],algparams);
%[mcmcparams, mcmcrun, mcmcruntheta, mcmcrunphi]=mcmc_ed_main3Dmodsel([dat(:,[1 2]) dphi],algparams);

if strcmp(datset,'Sim')
mcmcparams.truevalues=[params.mu params.sig(1)^-2 params.sig(2)^-2];
 else
   mcmcparams.truevalues=[];
end
mcmc.params.units=unts;
mcmcparams.savedir=SAVDIR;
mcmcparams.filename=FileName;
mcmcparams.alg=algparams;
mcmcparams.burnin=algparams.burnin;

if ~isempty(SAVDIR)
save([SAVDIR '/mcmcrun.mat'],'mcmcparams','mcmcrun','mcmcruntheta','mcmcrunphi');
end

plot_mcmcRun_basicvariables(mcmcparams,mcmcrun,algparams.burnin);

% Look at phi, theta variables.
  options=mcmcparams; % SAVDIR and filename
  options.param='phi';
if ~isempty(mcmcrunphi)
[poormovers_phi mvs]=plot_mcmc_ang(mcmcrunphi,dphi,algparams.burnin,50,options);     % phi
end
options.param='theta';
  [poormovers_theta mvs]=plot_mcmc_ang(mcmcruntheta,dat(:,5),algparams.burnin,50,options); % theta

  create_outputs(dat,mcmcparams,mcmcrun,10000,0);
