function [algparams convergencediagnos]=RunMCMC_ECsimulateddata(exptnum,runtreatment,specnam,filterpars,options)
  % [algparams convergencediagnos]=RunMCMC_ECsimulateddata(exptnum,runtreatment,specnam,filterpars,options)
  %
  % Driver (as function) for RunMCMC_EuclDist
%
%
% For simulated data runs. Analogous to EuclDistCorrection_Driver_fn
% Loads priors and data
%
% 
  % Loads data and sets up run.
  % From SimData_List give number of dataset. QUERY USED? runtreatment should be 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(assumes fields have tax, noc strings and anything else is assumed DMSO).  NOTE: Load extra categories associated with new treatments in loadDistData.m
  %
  % Distances measured in nm
  %
  % There are two filters. filterpars.filterdist only uses 3D data below this distance
  % filterpars.zfilterz used in comparisons (but not in MCMC itself)
  %
    % options (optional [] for none) include setting: a options.a (default 10), priors options.sigmeans=[20 40], length run options.nsteps, number chains options.nchains, to save thetas options.savethetas (=1 save, =0 dont save (default))
  %
  % options.toconvergence gives parameters of convergence criteria to run algorithm to convergence
    % See ReRunMCMC_eucldist for format
  %
  % Weak mu prior. Medium/strong sigx,sigz prior. a=3 has convergence problems
  %
  % Compare EuclDistCorrection_Driver_fn.m
  % NJB Oct 2018
  %

% Set path
path(path,genpath(pwd));
%path(path,'..');  % For simulate_3D
  
algstr='EDC';
ModNam='EuclDistMarg'; % Model name. Linked to RunMCMC_eucldist

SimData_List  % loads the SimData List
implot=1;

% ExptData.name ExptDat.fluorophore (GFP, Ab)
%
  dataid={exptnum,runtreatment};  % First is id below in FileNames{}, Second category is 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(Not done yet).  

%
% Run parameters. No need to edit below unless you want specific priors.
%

RunDIR='../MCMCruns';
initialisetype='priors'; %'true'; (For simulations only) % 'priors';
params.expansionfactor=1;    % For prior. Not being used, keep at 1.
nsamples=100000; % Target number of samples if chains long enough -- changes subsampling

if isempty(filterpars)
filterpars=struct('filterdist',200,'zfilterz',100);
end

if isfield(options,'a')
a=options.a;
specnam=[specnam '_a' num2str(a)]; % Adds a to specific name
 else
   a=10; % Parameter for the prior spread in the positional error/PSF.
end

if isfield(options,'mu')
   params.mu=options.mu;
else
params.mu=60;
end

if isfield(options,'sigmeans')
params.sig=options.sigmeans; % Rough values
else
params.sig=params.expansionfactor*[20 60]; % DEFAULT
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
   FileNameSh=[FileNames{dataid{1}}(1:(end-4)) '_' algstr]; 
   end

   ExptDat=ExptDats{dataid{1}}; % Metadata

%						      
% Saving directory for runs
%
SAVDIR=['Sim_' FileNameSh];

%  
% Load data. dat is matrix (r,th,phi,indexOnZthreshold)
%

%  [status msg]= mkdir([SimDataLoc 'Figures']);
%  [dats treats nams]=loadDistData([ExptDataLoc FileNames{dataid{1}}],dataid{2},filterpars,[]); % [ExptDataLoc 'Figures/' FileNameSh '_']); % A single or multiple structure dats{k}.3Ddat .treatment
  
% More primitive
dats={};simparams=[];
Srce=load([SimDataLoc FileNames{dataid{1}}]); % Loads dats treats simparams
if isfield(Srce,'simparams')
simparams=Srce.simparams;
 else
   if isfield(Srce,'params')
simparams=Srce.params;
end
end

if isfield(Srce,'dats')
dats=Srce.dats;
treats=Srce.treats;
else
dats{1}=Srce.dat;
treats{1}='NA';
end

unts={'nm','rad','rad'}; % Assumed original data is in microns

%
% MCMC Run length. marg (with twist) needs 200000. modsel 500000
%

if isfield(options,'nchains')
  algparams.numchains=options.nchains;
  else
algparams.numchains=4;
end

%length run options.nsteps, number chains options.nchains, to save thetas options.savethetas (=1 save (default), =0 dont save)

if isfield(options,'nsteps') &  options.nsteps>0
  algparams.n=ceil(options.nsteps);
else
algparams.n=100000;
end

algparams.subsample=ceil(algparams.n/nsamples);
algparams.burnin=ceil(algparams.n/2);
if isfield(options,'nblocks') &  options.nblocks>0
algparams.bloks=ceil(options.nblocks);
  else
algparams.bloks=5;
end
algparams.savdir=SAVDIR;
if isfield(options,'savethetas') & options.savethetas==1
algparams.savangles =1; % Save angles theta
disp('Saving inferred angles theta: this make files very large');
 else
   algparams.savangles =0; % Save angles theta
   end
algparams.scaleparams=params;
algparams.scaleparams.a=a;
algparams.initialisetype=initialisetype; % Initialisation of chain
algparams.paramnames={'mu','taux','tauz'};
algparams.unts=unts;
if ~isempty(simparams)
% mcmcparams.variables={'mu','taux','tauz'};
algparams.truevalues=[simparams.mu simparams.sig.^(-2)];
disp('Using true values \mu, \tau_x \tau_z')
disp(algparams.truevalues);
 else
algparams.truevalues=[];
end
algparams.filters=filterpars;

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

% FORCE PRIOR
% algparams.priors.param{2}=[42.6 1.9367e+04]; % Force prior in taux [a,k]


%algparams.priors.params={[80 10000000],[a 1/(params.sig(1)*gamma(a)/gamma(a-0.5))^2],[1 1/(params.sig(2)*gamma(a)/gamma(a-0.5))^2]}; % Strong priors

disp(['Priors on sdx, sdz are: rel error  ' num2str(100/sqrt(a)) '\% with means ' num2str(params.sig(1)) ', ' num2str(params.sig(2))]);  

%
% Main MCMC loop
%

rundirs={};
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
[convergencediagnos savdir]=RunMCMC_eucldist(dat,ExptDat,nam,algparams,RunDIR,1,[]);
rundirs{k}=savdir;

disp('Convergence statistics:');
disp(convergencediagnos.GR)

algparams.rundir=rundirs;

if isfield(options,'toconvergence')

%rerunparams=options.toconvergence;

convergencediagnos=ReRunMCMC_eucldist(algparams,RunDir,implot,modoptions,options.toconvergence);


end


end % k
