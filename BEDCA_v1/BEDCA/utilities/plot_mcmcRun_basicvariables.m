function plot_mcmcRun_basicvariables(mcmcparams,mcmcrun,burnin)
% plot_mcmcRun_sisterSwt(mcmcparams,mcmcrun,burnin)
%
% This plots the MCMC run basic variables in mcmcrun
%
% burnin can be specified or leave [] and uses default in mcmcparams
% typ='all'/[]DEFAULT, 'summary'
% 
% NJB 2014?

fontsize=20;

if isfield(mcmcparams,'varnames')
paramnams=mcmcparams.varnames;
 else
   paramnams=mcmcparams.variables;
end

truevalues=mcmcparams.truevalues;
SAVDIR=mcmcparams.savedir;
FileNam=mcmcparams.filename;
%dt=mcmcparams.timespacing;
%Dat=mcmcparams.datafile;

if ~isempty(SAVDIR)
FileNam=[FileNam 'Figs/'];
[status msg]=mkdir([SAVDIR '/' FileNam]);
end

if isempty(burnin)
burnin=mcmcparams.burnin;
end


%
% Plot histograms for each variable
%

for j=intersect(find(var(mcmcrun,0,1)>0),1:length(paramnams))

figure;
subplot(1,2,1)
plot(mcmcrun(:,j));
ylabel(paramnams{j});xlabel('sample number');
title('Full chain')
subplot(1,2,2)
plot(mcmcrun(burnin:end,j));
ylabel(paramnams{j});xlabel('sample number');
title('From burnin');


figure;
hist(mcmcrun(burnin:end,j),100)
hold on
str={'MCMC'};
if ~isempty(truevalues)
plot(truevalues(j),0,'r*','MarkerSize',15)
str{end+1}='TrueValues';
end
% Plot priors
xl=xlim;yl=ylim;
p=mcmcparams.priors.params{j};
[h x]=hist(mcmcrun(burnin:end,j),100);
hm=max(h);

switch mcmcparams.priors.types{j}
case 'Gaussian'
xp=min(xl(1),p(1)-4*min(p(2),20*(xl(2)-xl(1)))):(xl(2)-xl(1))/100:max(xl(2),p(1)+4*min(p(2),20*(xl(2)-xl(1))));
plot(xp,hm*normpdf(xp,p(1),p(2))/normpdf(p(1),p(1),p(2)),'r');
xlim(xl);
str{end+1}=['Prior (Mean ' num2str(p(1)) ')'];

 case 'Gamma' % x^(a-1) exp -bx
 mn=p(1)/p(2); sd=sqrt(p(1))/p(2);
xp=min(xl(1),mn-4*sd):(xl(2)-xl(1))/100:max(xl(2),mn+4*sd);
plot(xp,hm*gampdf(xp,p(1),1/p(2))/gampdf(p(1)/p(2),p(1),1/p(2)),'r');
xlim(xl)
str{end+1}=['Prior (Mean ' num2str(mn) ')'];

otherwise
disp(['No plotting for prior type ' mcmcparams.priors.types{j}])
end

  if length(str)==3
  legend(str{1},str{2},str{3});
  else
      legend(str{1},str{2});
end

title(['Posterior ' paramnams{j} ' mn, sd ' num2str(mean(mcmcrun(burnin:end,j))) ', ' num2str(std(mcmcrun(burnin:end,j))) ],'FontSize',fontsize)
  xlabel(paramnams{j},'FontSize',fontsize)
set(gca,'FontSize',fontsize);

if ~isempty(SAVDIR)
saveas(gcf,[SAVDIR '/' FileNam 'Posterior_' paramnams{j} '.fig'],'fig');
  print('-depsc',[SAVDIR '/' FileNam 'Posterior_' paramnams{j} '.eps']);
end

end %j


