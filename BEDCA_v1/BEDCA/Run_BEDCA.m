function [algparams convergencediagnos]=Run_BEDCA(exptnum,runtreatment,specnam,filterpars,options)
  % [algparams convergencediagnos]=Run_BEDCA(exptnum,runtreatment,specnam,filterpars,options)
  %
  % Driver function for RunMCMC_EuclDist
  % Load data and sets up run.
  % From ExptData_List give number of dataset. runtreatment should be 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(assumes fields have tax, noc strings and anything else is assumed DMSO).
  % NOTE: To load extra categories associated with new treatments, edit loadDistData.m
  %
  % Distances measured in nm
  %
  % There are two (optional) filters. filterpars.filterdist only uses 3D data below this distance
  % filterpars.zfilterz used in comparison plots (but not in MCMC itself)
  %
  % options (optional [] for none) include: a options.a (default 4, ie 50% relative error in prior precision tau_x, tau_z), priors options.sigmeans=[20 40], length run options.nsteps, number chains options.nchains, to save angles thetas options.savethetas (=1 save (default), =0 dont save, often useful if large data sets), option.notpaired =1 if KTs not paired 
  % options.checkonly Outputs info to screen but no run is started.
  % options.plotsvisible % Turns off plotting figures to screen with options.plotsvisible=0. Doesnt WORK on Matlab 2016
  % options.verbose =1 For verbose output of mcmc algorithm Default is 0 (not verbose)
  %    
  % options.toconvergence gives parameters of convergence criteria to run algorithm to convergence
  % See ReRunMCMC_eucldist for format
  %
  %
  % algparams.savangles =0; % Save angles theta
  %
  % Weak mu prior. Medium/strong sigx,sigz prior. a=3 has convergence problems
  %
  % NJB July 2017. Last edit March 2020 to add 'checkonly' if statement and storing .div in RunMCMC_eucldist
  %

  % set renderer
set(0,'DefaultFigureRendererMode','manual')
set(0,'DefaultFigureRenderer','OpenGl')

    printheader();
    
if ~isempty(options) & isfield(options,'plotsvisible') & options.plotsvisible==0 
  disp('Figure visibility is OFF. In case of premature exit turn visibility on with set(0,''DefaultFigureVisible'',''on'')');
disp(' ');
set(0,'DefaultFigureVisible','off');
end

algstr='BEDCA';
ModNam='EuclDistMarg'; % Model name. Linked to RunMCMC_eucldist

ExptData_list  % loads your ExptData List and TreatmentLib

% ExptData.name ExptDat.fluorophore (GFP, Ab)
%
  dataid={exptnum,runtreatment};  % First is id below in FileNames{}, Second category is 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(Not done yet).  Load extra categories associated with new treatments in loadDistData.m

%
% Run parameters. Dont edit below unless you want specific priors outside those available.
%

RunDIR='../MCMCruns';
initialisetype='priors'; %'true'; (For simulations only) % 'priors';
params.expansionfactor=1;    % For prior. Not being used, keep at 1.
nsamples=100000; % Target number of samples if chains long enough -- changes subsampling

if isempty(filterpars)
filterpars=struct('filterdist',200,'zfilterz',100);
 else
   if ~isfield(filterpars,'filterdist') filterpars.filterdist=200;end
   if ~isfield(filterpars,'zfilterz') filterpars.zfilterz=100;end
end

% this is parameter 'a' that sets the tau variance/mean ratio of prior.
if isfield(options,'a')
a=options.a;
specnam=[specnam '_a' num2str(a)]; % Adds a to specific name
 else
   a=4; % Parameter for the prior spread in the PSF.
end

if isfield(options,'mu')
   params.mu=options.mu;
else
params.mu=60;
end

if isfield(ExptDats{dataid{1}},'fluorophore')
switch ExptDats{dataid{1}}.fluorophore
case 'Ab' % Medium strength prior
params.sig=params.expansionfactor*[20 40]; % Rough values
case 'GFP' % Fluorescent GFP/mCherry pair Live.
params.sig=params.expansionfactor*[25 75]; % Rough values
%case 'Ab' 
%params.sig=params.expansionfactor*[15 45]; % Rough values
otherwise
disp('Need to select a fluorophore type to determine sigx, sigz priors');
end
 else
   disp('No fluorophore type specified (Ab, GFP), using defaults');
params.sig=params.expansionfactor*[20 60];
end

if isfield(options,'sigmeans')
params.sig=options.sigmeans; % Rough values
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
   ExptDat=ExptDats{dataid{1}};

% Set path
path(path,genpath(pwd));
%path(path,'..');  % For simulate_3D

%						      
% Saving directory for runs
%
SAVDIR=['Expt_' FileNameSh];

%  
% Load data. dat is matrix (r,th,phi,indexOnZthreshold)
%

[status msg]= mkdir([ExptDataLoc 'Figures']);

[dats treats nams]=loadDistDatatreatLib([ExptDataLoc FileNames{dataid{1}}],TreatmentLib,dataid{2},filterpars,[]);

if 1==0
bob TO EDIT
if isfield (options,'notpairedKT') & options.notpairedKT==1
  [dats treats nams]=loadDistDataOrig([ExptDataLoc FileNames{dataid{1}}],dataid{2},filterpars,[]);
else
    [dats treats nams]=loadDistData([ExptDataLoc FileNames{dataid{1}}],dataid{2},filterpars,[]); % [ExptDataLoc 'Figures/' FileNameSh '_']); % A single or multiple structure dats{k}.3Ddat .treatment
end

%  [ExptDataLoc FileNames{dataid{1}}]
%  filterpars
%  dats
end % 0==1
    
unts={'nm','rad','rad'}; % Assumed original data is in microns

%
% MCMC Run length. marg (with twist) needs 200000. modsel 500000
%

  if isfield(options,'verbose')
algparams.verbose=options.verbose;
  else
    algparams.verbose=0;
  end
  
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
algparams.burnin=ceil(algparams.n/2);  % Default burnin in 50%
if isfield(options,'nblocks') &  options.nblocks>0
algparams.bloks=ceil(options.nblocks);
  else
algparams.bloks=5;
end
algparams.savdir=SAVDIR;
if isfield(options,'savethetas') & options.savethetas==0
algparams.savangles =0; % Dont Save angles theta
 else
   algparams.savangles =1; % Save angles theta (default)
   end
algparams.scaleparams=params;
algparams.scaleparams.a=a;
algparams.initialisetype=initialisetype; % Initialisation of chain
algparams.paramnames={'mu','taux','tauz'};
algparams.unts=unts;
algparams.truevalues=[];
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

if ~isempty(options) & isfield(options,'checkonly') & options.checkonly==1

% This just outputs set up of run and data info
  for k=1:length(dats) % Over treatments

dat=dats{k};
if size(dat,2)>6
diversitymeasure_v2(dat,1);
end
end % k

convergencediagnos=[];

  else

rundirs={};
for k=1:length(dats) % Over treatments

dat=dats{k};
if size(dat,2)>6
diversitymeasure_v2(dat,1);
end
%diversitymeasure(dat(:,7));

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

%disp('Convergence statistics:');
%disp(convergencediagnos.GR)

algparams.rundir=rundirs;

if isfield(options,'toconvergence') % Not working. bob
%rerunparams=options.toconvergence;
convergencediagnos=ReRunMCMC_eucldist(algparams,RunDir,1,[]);
end




end % k treatments

end % checkonly

set(0,'DefaultFigureVisible','on');
