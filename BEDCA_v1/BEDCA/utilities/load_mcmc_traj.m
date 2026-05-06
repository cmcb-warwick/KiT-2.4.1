function mcmcdat=load_mcmc_traj(RunDIR,ExptDat,modnam,celnum,trajnum)
% mcmcdat=load_mcmc_traj(RunDIR,modnam,celnum,trajnum)
%
% Short cut driver for load_mcmc_runs for a single chain load
%
% NJB 

  mcmcdat=load_mcmc_runs([RunDIR '/' ExptDat.name '_Cell' num2str(celnum)],['MCMC_' modnam '_Traj' num2str(trajnum)],[],'single');
