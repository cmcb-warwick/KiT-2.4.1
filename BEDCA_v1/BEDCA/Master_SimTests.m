

%
%
%

% Set Renderer
set(0,'DefaultFigureRendererMode','manual')
set(0,'DefaultFigureRenderer','OpenGl')
  fontsize=20;

% Set path
path(path,genpath(pwd));
  
% Save location. [] for no sav
SAVDIR='../Figures/Figs_MCMCSimVary_mu/';

% Run MCMC
%for k=5:14
	 %[algparams convergencediagnos]=RunMCMC_ECsimulateddata(k,'','',[],[]);% struct('a',3));
% end


% Load data
filstem='Sim_simInflation_sample1000_v6_muvar_';
filmid='_EDC_mu=';
filend='_a3/';

muv=[0:2:8 10:10:70];
Hmns=figure;hold on;
Hmeds=figure;hold on;
Hsd=figure;hold on;

for k=1:length(muv)

mns=[];
meds=[];
sds=[];
for j=1:10  % Repeats
	A=load(['../MCMCruns/' filstem num2str(j) filmid num2str(muv(k)) filend 'MCMCStats_chain1.mat']);

% Extract statistics. {'mu'  'taux'  'tauz'  'sdx'  'sdz'}
% means
mns=[mns; A.tab.means];

% medians
meds=[meds; A.tab.medians];

% sd
sds=[sds; A.tab.sd];

end %j

figure(Hmns); %mu
errorbar(muv(k),mean(mns(:,1),1),std(mns(:,1),0,1),std(mns(:,1),0,1),'k');
vmns(k,:)=mean(mns,1);

figure(Hmeds); %mu
errorbar(muv(k),mean(meds(:,1),1),std(meds(:,1),0,1),std(meds(:,1),0,1),'k');
vmeds(k,:)=mean(meds,1);

figure(Hsd); %mu
errorbar(muv(k),mean(sds(:,1),1),std(sds(:,1),0,1),std(sds(:,1),0,1),'k');
vsds(k,:)=mean(sds,1);

end %k

figure(Hmns)
plot(muv,vmns(:,1),'k')
plot(muv,muv,'r--')
set(gca,'FontSize',fontsize)
xlabel('mu','FontSize',fontsize)
ylabel('Posterior mean mu','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Inferred_mumean_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Inferred_mumean_varymu.eps']);
end

figure(Hmeds)
plot(muv,vmeds(:,1),'k')
plot(muv,muv,'r--')
set(gca,'FontSize',fontsize)
xlabel('mu','FontSize',fontsize)
ylabel('Posterior median mu','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Inferred_mumed_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Inferred_mumed_varymu.eps']);
end

  
figure(Hsd)
set(gca,'FontSize',fontsize)
plot(muv,vsds(:,1),'k')
xlabel('mu','FontSize',fontsize)
ylabel('Posterior sd mu','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Inferred_musd_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Inferred_musd_varymu.eps']);
end

  
% Plot representative histograms
j=2;
%Select k
kv=[1 2 4 6 7 9 12];

for i=1:3
	H(i)=figure;hold on
end
	col={'r','m','g','c','b','k','r--','m--'}

i=1;
for k=kv
B=load(['../MCMCruns/' filstem num2str(j) filmid num2str(muv(k)) filend '/MCMC_EuclDistMargoutput1.mat']);
burnin=B.mcmcparams.burnin;

for s=1:3
       [h x]=hist(B.mcmcrun(burnin:end,s),200);
figure(H(s));
       plot(x,h/(sum(h)*(x(2)-x(1))),col{i});
end
       str{i}=['mu=' num2str(muv(k))];
i=i+1;
end % k


figure(H(1));
legend(str)
set(gca,'FontSize',fontsize)
xlabel('Posterior mu','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Posterior_mu_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Posterior_mu_varymu.eps']);
end


figure(H(2));
legend(str)
set(gca,'FontSize',fontsize)
xlabel('Posterior tau_x','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Posterior_taux_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Posterior_taux_varymu.eps']);
end



figure(H(3));
legend(str)
set(gca,'FontSize',fontsize)
xlabel('Posterior tau_z','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Posterior_tauz_varymu.fig'],'fig');
  print('-depsc',[SAVDIR 'Posterior_tauz_varymu.eps']);
end

% this case has mu tail mu=0
  k=1;
  B=load(['../MCMCruns/' filstem num2str(j) filmid num2str(muv(k)) filend '/MCMC_EuclDistMargoutput1.mat']);
twoDhisto(B.mcmcrun(:,1:2)',[200,200],[0 35 1.3/1000 3.5/1000],'',[])
xlabel('mu','FontSize',fontsize)
ylabel('taux','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'JointPosterior_mutaux_mu=' num2str(muv(k)) '.fig'],'fig');
  print('-depsc',[SAVDIR 'JointPosterior_mutaux_mu=' num2str(muv(k)) '.eps']);
end

close
xlabel('mu','FontSize',fontsize)
ylabel('taux','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'JointPosteriorlogscale_mutaux_mu=' num2str(muv(k)) '.fig'],'fig');
  print('-depsc',[SAVDIR 'JointPosteriorlogscale_mutaux_mu=' num2str(muv(k)) '.eps']);
end


k=6;
  B=load(['../MCMCruns/' filstem num2str(j) filmid num2str(muv(k)) filend '/MCMC_EuclDistMargoutput1.mat']);
twoDhisto(B.mcmcrun(:,1:2)',[200,200],[0 35 1.3/1000 3.5/1000],'',[])
xlabel('mu','FontSize',fontsize)
ylabel('taux','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'JointPosterior_mutaux_mu=' num2str(muv(k)) '.fig'],'fig');
  print('-depsc',[SAVDIR 'JointPosterior_mutaux_mu=' num2str(muv(k)) '.eps']);
end

close
xlabel('mu','FontSize',fontsize)
ylabel('taux','FontSize',fontsize)

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'JointPosteriorlogscale_mutaux_mu=' num2str(muv(k)) '.fig'],'fig');
  print('-depsc',[SAVDIR 'JointPosteriorlogscale_mutaux_mu=' num2str(muv(k)) '.eps']);
end


