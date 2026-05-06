function paper_figs(exptnum,runtreatment,nam,options)
% paper_figs(exptnum,runtreatment,nam,options)
%
% Give specific directory in MCMCruns or experiment number/treatment 
% eg paper_figs([],[],'Expt_iM_180903_Ndc80GFPC_DMSO_selectedSpots9G3_Ndc80GFP_EDC_DMSO',[])
  % or paper_figs(6,'DMSO',[],[])
%
%
% options to control plot
  % options.boxs options.linewidth
%  
% This plots figures for paper and saves in Figures in aboce directory
%
% NJB april 2019

% set renderer
set(0,'DefaultFigureRendererMode','manual')
set(0,'DefaultFigureRenderer','OpenGl')

  
% plot defaults
unts={'nm','rad','rad'}; % Assumed original data is in microns

fontsize=20;
linewidth=2;
boxs=25;
if ~isempty(options)
if isfield(options,'boxs')
boxs=options.boxs;
end
if isfield(options,'linewidth')
linewidth=options.linewidth;
end
end

% Load Expt list

% From EuclDist_Driver
ExptData_list
RunDIR='../MCMCruns';
FigDIR='Figures';

algstr='EDC';

if ~isempty(exptnum) && exptnum>0
dataid={exptnum,runtreatment};  % First is id below in FileNames{}, Second category is 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO'/'untreated', 'all'(Not done yet).  Load extra categories associated with new treatments in loadDistData.m

brkstr=max(strfind(FileNames{dataid{1}},'/'));
if ~isempty(brkstr) & brkstr(end)==length(FileNames{dataid{1}})
  brkstr(end)=[];
end

if ~isempty(brkstr)
FileNameSh=[FileNames{dataid{1}}((brkstr(end)+1):(end-4)) '_' algstr]; % Creates a Short file name (Last Part Of FileName), removing .mat
 else
   FileNameSh=[FileNames{dataid{1}}(1:(end-4)) '_' algstr];
end
     
% Build treatment name
if ~isempty(strfind(lower(nam),'tax'))
treatnam='Taxol';
 else
   if ~isempty(strfind(lower(nam),'noc'))
treatnam='Nocod';
 else
   treatnam='DMSO';
end
end

ExptDat=ExptDats{dataid{1}};


%						      
% Saving directory for mcmcruns
%
SAVDIR=['Expt_' FileNameSh '_' treatnam '/'];
SAVDIR=[RunDIR '/' SAVDIR];

 else
%						      
% Saving directory for mcmcruns, given directory in MCMCruns
%

   SAVDIR=[RunDIR '/' nam '/'];
end

disp(['Loading rundata for ' SAVDIR]);

%load([SAVDIR 'MCMCStats_chain1.mat']);      % Load stats of first chain
load([SAVDIR 'MCMC_EuclDistMargoutput1.mat']); % Load first chain 

paramnams=mcmcparams.variables;
truevalues=mcmcparams.truevalues;
burnin=mcmcparams.burnin;
mcmcrun=mcmcrun(burnin:end,:); % Delete burnin
dat=mcmcparams.datafile;

%
% Plot data and estimated mu
%
figure;hold on;
J=find(dat(:,6));
[h x]=hist(dat(J,1),boxs);
hm=max(h);
%[h x]=hist(dat(:,1),boxs);
plot(x,h*hm/max(h),'k','LineWidth',linewidth); % 3D raw data
% Stats
y=ylim;
%plot(median(dat(J,1))*[1 1],y,'c--');

% mu posterior
j=1;
[h x]=hist(mcmcrun(:,j),boxs);
plot(x,h*hm/max(h),'r','LineWidth',linewidth);

legend('Expt 3D (filtered)','\Delta_{ED}');

% Plot means
plot(mean(dat(J,1))*[1 1],y,'k:','LineWidth',linewidth);
%plot(tab.means(1)*[1 1],y,'r:','LineWidth',linewidth);

title('3D measurements and \Delta_{EC}','FontSize',fontsize)
xlabel('3D distance','FontSize',fontsize)
set(gca,'FontSize',fontsize);
xlim([0 150])


if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Figures/PFig_DataWithPostEsts.fig'],'fig');
  print('-depsc',[SAVDIR 'Figures/PFig_DataWithPostEsts.eps']);
end







