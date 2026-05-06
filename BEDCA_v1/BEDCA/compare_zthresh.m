

%
% driver for mcmc_ed_main to compare with/without z threshold.
% Based on mcmc_ed
% Loads data and sets up run.
%
% Distances measured in nm
%
% NJB May 2017


initialisetype='priors';

% Set path
path(path,genpath(pwd))

% Load data. dat is matrix (r,th,phi,indexOnZthreshold)
FileName='simInflation_sample1000_v5';

load(['../SimData/' FileName '.mat']);
disp(['Data file = ' FileName]); 
SAVDIR=['../MCMCruns/Sim_' FileName '_zfilter'];

if ~isempty(SAVDIR)
mkdir(SAVDIR);
end

J=find(dat(:,4)==1); % Satisfies z threshold
% Filter data
dat=dat(J,:);


% data processing
% ensure thetaX>0
  J=find(dat(:,2)<0);
if ~isempty(J)
disp(['Warning: ' num2str(length(J)) ' data points have negative theta. Taking abs value.']);
dat(J,2)=-dat(J,2);
end
  
algparams.n=40000;
algparams.subsample=1;
algparams.burnin=10000;

%
% Priors
% For Tau dist Gamma(a,b) E[sigma]=b^1/2 Gamma(a-1/2)/Gamma(a)
  % Variance goes as var(sigma) \approx 1/(2 (ab)^.5)
  %
algparams.priors.variable={'mu','taux','tauz'};
algparams.priors.types={'Gaussian','Gamma','Gamma'}; % Gaussian mean variance,...  Gamma power decaycoef
  
algparams.priors.params={[80 10000000],[1 0.01],[1 0.01]}; % Weak priors

a=100;
algparams.priors.params={[80 10000000],[a 1/(params.sig(1)*gamma(a)/gamma(a-0.5))^2],[1 1/(params.sig(2)*gamma(a)/gamma(a-0.5))^2]}; % Strong priors


dphi=dat(:,6)-dat(:,3);
J=find(dphi>pi);if ~isempty(J) dphi(J)=dphi(J)-2*pi;end
J=find(dphi<-pi);if ~isempty(J) dphi(J)=dphi(J)+2*pi;end

switch initialisetype
case 'true'
% Load initial conditions.
algparams.init=[params.mu params.sig];   % True values
algparams.initthetas=dat(:,5);           % True thetas vec n
algparams.initphis=dphi;            % True values vec n

 case 'priors'  % Draw from priors.

 % mu. Positive.
mu=prior_draw(algparams.priors,'mu');
 while (mu < 0)
   mu=prior_draw(algparams.priors,'mu');
end

algparams.init=[mu 1/sqrt(prior_draw(algparams.priors,'taux')) 1/sqrt(prior_draw(algparams.priors,'tauz'))];
algparams.initthetas=pi*rand(size(dat,1),1);
algparams.initphis=2*pi*rand(size(dat,1),1)-pi;
 
 otherwise
 disp(['No such initialisation type as ' initialisetype]);
disp('ABORTING.');
return

end

algparams.initialise=initialisetype;

[mcmcparams, mcmcrun, mcmcruntheta, mcmcrunphi]=mcmc_ed_main3D(dat,algparams);

mcmcparams.truevalues=[params.mu params.sig(1)^-2 params.sig(2)^-2];
mcmcparams.savedir=SAVDIR;
mcmcparams.filename=[FileName '_zfilter'];

if ~isempty(SAVDIR)
save([SAVDIR '/mcmcrun_zfilter.mat'],'mcmcparams','mcmcrun','mcmcruntheta','mcmcrunphi');
end

plot_mcmcRun_basicvariables(mcmcparams,mcmcrun,algparams.burnin);


% Look at phi, theta variables.
  options=mcmcparams; % SAVDIR and filename
  options.param='phi';
  [poormovers_phi mvs]=plot_mcmc_ang(mcmcrunphi,dphi,algparams.burnin,50,options); %         phi
    options.param='theta';
  [poormovers_theta mvs]=plot_mcmc_ang(mcmcruntheta,dat(:,5),algparams.burnin,50,options); % theta

