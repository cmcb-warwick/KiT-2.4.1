function convergencediagnos=ReRunMCMC_eucldist(rundata,MCMCRunDir,implot,toconvergence)
% convergencediagnos=ReRunMCMC_eucldist(rundata,MCMCRunDir,implot,toconvergence)
%
% ReRuns the MCMC using last state in mcmcparams 
% Provide the MCMCRunDir location
%
% runparams is used to check number chains etc%
%
% toconvergence.trials (number of retrys), .Rc (criteria for convergence on given statistic .R, .Rc, .MPSRF, multiple allowed), .thinningrate (2-10, factor by which to extend run)
% 
% implot=1 To plot figures
%
% Returns a convergence diagnostic structure (test_convergence)
%
%
% Burroughs. Aug 2018. 

  if isempty(toconvergence)
  toconvergence.trials=4;
toconvergence.Rc=1.05; % (criteria for convergence on given statistic .R, .Rc, .MPSRF, multiple allowed),
  toconvergence.thinningrate=3;
end

  
% Loads convergence diagnostic
ffc=ls([MCMCRunDir '/*convergencediag.mat'])  % Lists all .mat files
  %  class(ffc), char 
ConvD=load(ffc)

%
  ff=ls([MCMCRunDir '/*output.mat']);  % Lists all .mat files

  % Check against rundata

  if ~isempty(rundata)

  if length(ff)~=rundata.numchains;
disp('Number of chains doesnt match rundata. ABORT');
return;
end

end


% Rerun MCMC
dstamp=[];
for k=1:length(ff)

	A=load(ff{k});

if isempty(dstamp)
dstamp=A.mcmcparams.datestamp;
 else
   if sum(dstamp ~= A.mcmcparams.datestamp)~=0
   disp(['Datestamp mismatch on chain ' numstr(k) '. ABORT. Rerun Expt from beginning, output corrupted.']);
return;
end
end

ncurr=size(A.mcmcrun,1);
subsamplecurr=A.mcmcrun.subsample;

rundata=mcmcparams.alg;
rundata.init=A.mcmcparams.laststate;
rundata.initthetas=A.mcmcparams.laststatethetas;

rundata.n=A.mcmcparams.runlength*A.mcmcparams.subsample*(toconvergence.thinningrate-1);
rundata.burnin=ceil(A.mcmcparams.runlength*A.mcmcparams.subsample*(toconvergence.thinningrate-1)/2);


disp(['Burnin ' num2str(burnin) ', subsampling rate ' num2str(subsample)]);

% Restart chain with same subsampling
[rrmcmcparams, rrmcmcrun, rrmcmcruntheta, rrmcmcrunphi]=mcmc_ed_main3Dmarg(dat(:,[1 2]),rundata);

% Glue parts together.
mcmcrun=[A.mcmcrun; rrmcmcrun];
mcmcruntheta=[A.mcmcruntheta; rrmcmcruntheta];

mcmcrun=mcmcrun(1:toconvergence.thinningrate:end,:); % Thin

mcmcruntheta=mcmcruntheta(1:toconvergence.thinningrate:end,:); % Thin if being saved

% Create mcmcparams
mcmcparams.toconvergence=toconvergence;
mcmcparams.toconvergence.subsamplesorig=subsamplecurr;
mcmcparams.toconvergence.norig=ncurr;

mcmcparams.subsample=mcmcparams.subsample*toconvergence.thinningrate;
mcmcparams.burnin=ceil((rundata.burnin+ncurr)/toconvergence.thinningrate);

end % k

convergencediagnos=[];

