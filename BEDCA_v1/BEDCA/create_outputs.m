function tab=create_outputs(dat,mcmcparams,mcmcrun,nsim,SAV)
  % tab=create_outputs(dat,mcmcparams,mcmcrun,nsim,SAV)
  %
  % this plots outputs and compares data to a reconstruction from a simulation.
  % Returns a table of statistics
  %
  % nsim is number of simulated 3D vectors to generate
  % Give zfilter if you want to compare to a filter 3D data set
  %
  % NJB May 2017

fontsize=20;
boxs=50;

paramnams=mcmcparams.variables;
truevalues=mcmcparams.truevalues;

if ~isempty(SAV) & SAV ~=0
SAVDIR=mcmcparams.savedir;
%FileNam=mcmcparams.filename;
 else
   SAVDIR=[];
end

if ~isempty(SAVDIR) 
FileNam='Figures/';
[status msg]=mkdir([SAVDIR '/' FileNam]);
end

burnin=mcmcparams.burnin;
mcmcrun=mcmcrun(burnin:end,:); % Delete burnin

%
% Plot data and inferred mu distribution
%
figure;hold on;
%J=find(dat(:,6));  % z filter flag
%hist(dat(J,1),boxs);
%
hist(dat(:,1),boxs);
[h x]=hist(dat(:,1),boxs,'Color','c');
hm=max(h);
% Stats
y=ylim;
plot(median(dat(:,1))*[1 1],y,'k--');
plot(mean(dat(:,1))*[1 1],y,'k:');

% mu posterior
j=1;
[h x]=hist(mcmcrun(:,j),ceil(boxs/5));
plot(x,h*hm/max(h),'r');
title('3D distance data and inferred distance histograms','FontSize',fontsize)
xlabel('3D distance','FontSize',fontsize)
set(gca,'FontSize',fontsize);

legend('Measured 3D','Median data','Mean data','MCMC');

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR '/' FileNam 'DataWithPostEsts.fig'],'fig');
  print('-depsc',[SAVDIR '/' FileNam 'DataWithPostEsts.eps']);
end

%
% Plot histograms for each variable
%

figure;
subplot(1,3,1)
j=1;
hist(mcmcrun(:,j),boxs);
hold on
if ~isempty(truevalues)
plot(truevalues(j),0,'r*','MarkerSize',15)
end
% Plot prior

title(['Posterior mu'],'FontSize',fontsize)
  xlabel('mu','FontSize',fontsize)
set(gca,'FontSize',fontsize);

subplot(1,3,2)  % taux/sigmax
j=2;
%[h x]=hist(1./sqrt(mcmcrun(:,j)),boxs);
%plot(x,h/(sum(h)*(x(2)-x(1))),'k');
hist(1./sqrt(mcmcrun(:,j)),boxs);
hold on
if ~isempty(truevalues)
plot(1/sqrt(truevalues(j)),0,'r*','MarkerSize',15)
end
x=xlim();
[params typ]=get_priorparams(mcmcparams.priors,'taux');
  xv=x(1):(x(2)-x(1))/1000:x(2);
plot(xv,gampdf(1./(xv.^2),params(1),params(2))./(xv.^3),'b')

title(['Posterior sigx'],'FontSize',fontsize)
xlabel('sigx','FontSize',fontsize)
set(gca,'FontSize',fontsize);

subplot(1,3,3)
j=3;
hist(1./sqrt(mcmcrun(:,j)),boxs)
hold on
if ~isempty(truevalues)
plot(1/sqrt(truevalues(j)),0,'r*','MarkerSize',15)
end
title(['Posterior sigz'],'FontSize',fontsize)
xlabel('sigz','FontSize',fontsize)
set(gca,'FontSize',fontsize);

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR '/' FileNam 'PosteriorParams.fig'],'fig');
  print('-depsc',[SAVDIR '/' FileNam 'PosteriorParams.eps']);
end


% Create table d3D (measurements)  mu sigx, sigz  (mean, meadian, sd of parameters)

disp(['Inflation in mu : inferred distance ' num2str(mean(mcmcrun(:,1))) ' shifts to observed mean 3D distance ' num2str(mean(dat(:,1))) ', (median ' num2str(median(dat(:,1))) ').']); 

if isfield(mcmcparams,'datestamp')
tab.datestamp=mcmcparams.datestamp;
end
if isfield(mcmcparams,'chain')
tab.chainnumber=mcmcparams.chain;
end
tab.samplen=[burnin size(mcmcrun,1)];
tab.params=[mcmcparams.variables {'sdx','sdz'}];
tab.means=[mean(mcmcrun) mean(sqrt(1./mcmcrun(:,2:3)))];
tab.medians=[median(mcmcrun) median(sqrt(1./mcmcrun(:,2:3)))];
tab.sd=[std(mcmcrun) std(sqrt(1./mcmcrun(:,2:3)))];
tab.datmean=mean(dat(:,1));
tab.datmedian=median(dat(:,1));
tab.datsd=std(dat(:,1));
tab.datskew=skewness(dat(:,1),0); %moment(dat(:,1),3)/tab.datsd^3;
tab.samplesize=size(dat,1);
tab.inflationmedian=tab.datmedian/tab.medians(1);

% Use dat column 6, zfilter satisfied
  J=find(dat(:,6)); 
tab.datzfiltermean=mean(dat(J,1));
tab.datzfiltermedian=median(dat(J,1));
tab.datzfiltersd=std(dat(J,1));
tab.filteredsamplesize=length(J);
tab.inflationzfiltermedian=tab.datzfiltermedian/tab.medians(1);

% Recreate distribution d from simulation with 95% confidence based on Poisson
figure;hold on;
[h x]=hist(dat(:,1),boxs);
plot(x,h/(sum(h)*(x(2)-x(1))),'k');
xlabel('3D distance d','FontSize',fontsize);

% Perform simulation with theta restricted to that observed.
%datsim=simulate_3Ddataset(mean(mcmcrun(:,1)),1/sqrt(mean(mcmcrun(:,2))),1/sqrt(mean(mcmcrun(:,3))),nsim,100,dat(:,2));

if isempty(nsim) | nsim==0
return;
end

datsim=simulate_3DdataFromMCMC_v2(mcmcrun,nsim,100,[]);
tab.simReconstruction=datsim;

[h x]=hist(datsim(:,1),250);
plot(x,h/(sum(h)*(x(2)-x(1))),'r');

legend('Expt','Sim');
set(gca,'FontSize',fontsize);
title('3D distance with theta matched simulation');
if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR '/' FileNam 'DataSimComp.fig'],'fig');
  print('-depsc',[SAVDIR '/' FileNam 'DataSimComp.eps']);
end

