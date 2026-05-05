


%
% Driver for RunMCMC_EuclDist
% Loads data and sets up run.
%
% Distances measured in nm
%
% Weak mu prior. Medium/strong sigx,sigz prior. a=3 has convergence problems
%
% NJB July 2017
%
% ExptData.name ExptDat.fluorophore (GFP, Ab)
%

%
% Choose one of datasets to run
%
  dataid={4,'DMSO'};  % First is id below in FileNames{}, Second category is 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(Not done yet).  Load extra categories associated with new treatments in loadDistData.m
   specnam='_a4'; % Add a specific nam, eg relating prior

   % You can change these filters
filterpars.filterdist=250; % max 3D distance allowed in sample data
filterpars.zfilterz=100; % zfilter value. Comparison only, not used in MCMC. 
   
% These are the data sets that are available. Files in ExptDataLoc
% Each has a structure with treatment categories (DMSO, tax15min, noc2hr)

ExptDataLoc= '../ExptData/';
FileNames{1}='NormMeanZero/intraMeasurements_HeLaTwoColourCENPA_central.mat'; % Dual label
FileNames{2}='NormMeanZero/intraMeasurements_HeLaCENPANdc80_central.mat';
FileNames{3}='NormMeanZero/intraMeasurements_RPECENPCNdc809G3_central.mat';
FileNames{4}='NormMeanZero/intraMeasurements_RPECENPCBub1_central.mat';
FileNames{5}='NormMeanZero/intraMeasurements_RPECENPC488594_central.mat'; % Dual label. Untreated only.


% Add more if needed.  


% Create an associated ExptDat
ExptDats{1}.name='HeLaTwoColourCENPA'; ExptDats{1}.fluorophore='Ab';
ExptDats{2}.name='HeLaCENPANdc80'; ExptDats{2}.fluorophore='GFP';
ExptDats{3}.name='RPECENPCNdc809G3'; ExptDats{3}.fluorophore='Ab';
ExptDats{4}.name='RPECENPCBub1'; ExptDats{4}.fluorophore='Abmedium';
ExptDats{5}.name='RPECENPCTwoAb'; ExptDats{5}.fluorophore='Ab';


%
% Run parameters. No need to edit below unless you want specific priors.
%

RunDIR='../MCMCruns';
initialisetype='priors'; %'true'; (For simulations only) % 'priors';
params.expansionfactor=1;    % For prior
algstr='ed';

params.mu=60;
switch ExptDats{dataid{1}}.fluorophore
case 'Ab'
params.sig=params.expansionfactor*[15 45]; % Rough values
a=100; % Parameter for the prior spread in the PSF.
 case 'GFP' % Fluorescent GFP/mCherry pair Live.
params.sig=params.expansionfactor*[25 75]; % Rough values
a=100; % Parameter for the prior spread in the PSF.
case 'Abmedium' % Medium strength prior
params.sig=params.expansionfactor*[15 45]; % Rough values
a=4; % Parameter for the prior spread in the PSF.
otherwise
disp('Need to select a fluorophore type to deetermine sigx, sigz priors');
end

%
% Shortened name
brkstr=max(strfind(FileNames{dataid{1}},'/'));
if ~isempty(brkstr) & brkstr(end)==length(FileNames{dataid{1}})
  brkstr(end)=[];
end

if ~isempty(brkstr)
FileNameSh=[FileNames{dataid{1}}((brkstr(end)+1):(end-4)) '_' algstr]; % Creates a Short file name (Last Part Of FileName), removing .mat
 else
FileNameSh=[FileNames{dataid{1}} '_' algstr]; 
   end
   ExptDat=ExptDats{dataid{1}};

   
% Set path
  path(path,genpath(pwd));
path(path,'..');  % For simulate_3D


%  
% Load data. dat is matrix (r,th,phi,indexOnZthreshold)
%
  
[dats treats nams]=loadDistData([ExptDataLoc FileNames{dataid{1}}],dataid{2},filterpars); % A single or multiple structure dats{k}.3Ddat .treatment

% More primitive
  % A=load([ExptDataLoc FileNames{dataid{1}}]);
  
unts={'nm','rad','rad'}; % Assumed original data is in microns

%						      
% Saving directory for runs
%
SAVDIR=['Expt_' FileNameSh];

%
% Priors 
%

%
% MCMC Run length. marg (with twist) needs 200000. modsel 500000
%


algparams.numchains=4;
algparams.n=100000;
algparams.subsample=1;
algparams.burnin=50000;
algparams.bloks=5;
algparams.savdir=SAVDIR;
algparams.scaleparams=params;
algparams.scaleparams.a=a;
algparams.initialisetype=initialisetype; % Initialisation of chain
algparams.paramnames={'mu','taux','tauz'};
algparams.unts=unts;
algparams.truevalues=[];
algparams.filters=filterpars;

ModNam='EuclDistMarg'; % Model name. Linked to RunMCMC_eucldist


%algparams.init_hMC % You may want to set this if simulated data

%
% Priors
%
% Gamma G(a,b) has mean a/b, var a/b^2. Rel error is 1/sqrt(a).  
% sig has mean sqrt(b/a), and variance (a/b^2)/(4 a^3/b^3)=b/(4a^2)
% Prior tau_x Need b/a = 25.5^2=650, for rel error 10% -> a=100, b=65000
%
 
algparams.priors.variable={'mu','taux','tauz'};
algparams.priors.types={'Gaussian','Gamma','Gamma'}; % Gaussian mean variance,...  Gamma power decaycoef not decaylength as used in matlab
algparams.priors.params={[params.mu 10000],[a params.sig(1)^2*a],[a params.sig(2)^2*a]}; % 10% error priors 
%algparams.priors.params={[70 10000000],[1. 1/100],[1 1/100]}; % Weak priors for Gamma are [1 0.01]

%algparams.priors.params={[80 10000000],[a 1/(params.sig(1)*gamma(a)/gamma(a-0.5))^2],[1 1/(params.sig(2)*gamma(a)/gamma(a-0.5))^2]}; % Strong priors

disp(['Priors on sdx, sdz are: rel error  ' num2str(100/sqrt(a)) '\% with means ' num2str(params.sig(1)) ', ' num2str(params.sig(2))]);  

%
% Main MCMC loop
%

for k=1:length(dats) % Over treatments

 dat=dats{k};

% data processing
% ensure thetaX>0
J=find(dat(:,2)<0);
if ~isempty(J)
disp(['Warning: ' num2str(length(J)) ' data points have negative theta. Taking abs value.']);
dat(J,2)=-dat(J,2);
end

disp(['Data set ' FileNames{dataid{1}} ' with treatment '  treats{k}]);

nam=[treats{k} specnam];

convergencediagnos=RunMCMC_eucldist(dat,ExptDat,nam,algparams,RunDIR,1,[]);

% This loads the runs and stats if only one run. You can comment out.
if length(dats)==1
mcmcdat=load_mcmc_runs([RunDIR '/' algparams.savdir '_' nam],['MCMC_' ModNam],[],'single');
[mcmcparams mcmcrun mcmcruntheta mcmcrunphi]=unpack_mcmcdat(mcmcdat);
load([RunDIR '/' algparams.savdir '_' nam '/MCMCStats.mat']);  % loads tab
end

end % k
