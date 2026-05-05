function Trajredo=report_convergence(TrajCells,RunDIR,ExptDat,modname,thresholds,recompute,saveflag)
% Trajredo=report_convergence(TrajCells,RunDIR,ExptDat,modname,thresholds,recompute,saveflag)
%
% This assessess convergence on available diagnostics
% reads convergence files, or recomputes if recompute=1 (saves if saveflag=1)
%
% NJB July 2013


if isempty(RunDIR)
  RunDIR='MCMCruns';
end

if isempty(recompute)
  recompute=0;
if isempty(saveflag)
  saveflag=0;
end

end

switch modname

case {'polewardhMC_SisterSwitcher','poleMShMC','poleMShMC_Sisters','polehMC_SisterSwitcher','Poleward MShMC:multiple poleward hMC sister switching algorithm','polewardhMC_SisterSwitcherproj','poleMShMCproj','poleMShMC_SistersProj','Proj Poleward MShMC:multiple poleward hMC sister switching algorithm','Projected SISTER poleward hMC switcher Model A: coherent/incoherent model'}


numfailed=[];
cnt=0;
for k=1:length(TrajCells) % Run over cells

Direc=[RunDIR '/' ExptDat.name '_Cell' num2str(TrajCells(k).cellIdx)];

FList=load_mcmc_runs(Direc,[],[],[]); % List trajectory names, over all models

for f=FList % Run over Trajectories

% You may need to reset modname to allow for user choice 
if length(f{1})> 5+length(modname) & strcmp(f{1}(1:(5+length(modname))),['MCMC_' modname]) % Correct model

cnt=cnt+1;

if recompute % Recompute diagnostics
close all;
mcmcdiagnos=test_convergencemcmcRuns(Direc,f{1},[],[],1,1); % test convergence of model

if saveflag & ~isempty(mcmcdiagnos)
  save([Direc '/' f{1} 'convergencediag.mat'],'mcmcdiagnos');
disp(['Saved run at ' Direc '/' f{1}  'convergencediag.mat']);
end

 else
[mcmcdat mcmcdiagnos]=load_mcmc_runs(Direc,f{1},[],'all');
end % recompute

disp(mcmcdiagnos.GR.Rc)  % Prints Rc output
failedconv=list_unconvergedruns(mcmcdiagnos,'GR',thresholds);
loc=strfind(f{1},'Traj');  % Unravel trajectory number from name
trajid=str2num(f{1}(loc(end)+4:end));
numfailed(cnt,:)=[TrajCells(k).cellIdx trajid  length(failedconv) mcmcdiagnos.GR.Rc];

end

end %f


close all;
disp(['Done Cell ' num2str(TrajCells(k).cellIdx)]);


end %k

% Plot failures by cell
cellid=unique(numfailed(:,1));
cnts=[];
for c=cellid'
  J=find(numfailed(:,1)==c);
cnts=[cnts; [length(J) sum(numfailed(J,3))]];
end %c

cnts
figure;bar(cnts)
legend('No. Traj','Not Conv');
xlabel('Cell');

% Save figure and convergence analysis
print('-depsc',[RunDIR '/convergencebycell_' modname '.eps']);
save([RunDIR '/convergence_' modname '.mat'],'numfailed');

disp(['Failure fraction: ' num2str(sum(numfailed(:,3))/size(numfailed,1)) ' at threshold ' num2str(thresholds(1))])


% numfailed is a matrix with CellIdx TrajId Failure

J=find(numfailed(:,3));
if ~isempty(J)
Trajredo=retrieve_TrajData(TrajCells,dataset({numfailed(J,1),'Cell'},{numfailed(J,1),'Traj'}));
end

case {'1D_harmonicwellMCMC','nocadozoletreated','MCMC1dHW'}

ModNam='MCMC1dHW';

for k=1:length(TrajCells)

disp('Not written');
bob

end %k


otherwise

disp('No such model');
end





