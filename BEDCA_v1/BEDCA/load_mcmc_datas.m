function D=load_mcmc_datas(exptnum,runtreatment,specnam)
%
  %
  %
  %
% NOTE.
% This loads the runs and stats if only one run. You can comment out.


Not DONe

  SAVDIR=['Expt_' FileNameSh];

brkstr=max(strfind(FileNames{dataid{1}},'/'));
if ~isempty(brkstr) & brkstr(end)==length(FileNames{dataid{1}})
  brkstr(end)=[];
end

  
  
mcmcdat=load_mcmc_runs([RunDIR '/' algparams.savdir '_' nam],['MCMC_' ModNam],[],'single');
[mcmcparams mcmcrun mcmcruntheta mcmcrunphi]=unpack_mcmcdat(mcmcdat);
load([RunDIR '/' algparams.savdir '_' nam '/MCMCStats.mat']);  % loads tab
end
