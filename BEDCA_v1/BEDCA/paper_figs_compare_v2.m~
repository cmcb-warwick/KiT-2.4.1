function paper_figs_compare(exptnum,runtreatments,server,options)
% paper_figs_compare(exptnum,runtreatment,server,options)
%
  % This plots the posteriors for the runtreatments
  %
% Give  experiment number/treatment from ExpData_List 
% eg  paper_figs_compare(6,'DMSO',[],[])
%
  % Give treatments to plot {'DMSO','taxol','Nocod'}
%
% options to control plot
% options.boxs options.linewidth options.xmax (maximum of x in plot,
% default 150, also allows options.xmax='auto' for autoscaling)
% options.legend == 0 turns off legends. On by default.
% options.raw  Plot raw delta3D
  %
% This plots figures for paper and saves in Figures in first treatment directory as PFig_CompareTreatmentsPost
%
  % Related paper_figs
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

if exist('MCMCruns','dir')==7
RunDIR='MCMCruns'; % Different NJB '../MCMCruns', Emanuele 'MCMCruns'
else
RunDIR='../MCMCruns';
end
  FigDIR='Figures';

if ~isempty(server)
    switch server
        case {'McAinshLabindiv','Ceresindiv'}
            RunDIR=McAinshLabDIR;
            case {'McAinshLabpool','Cerespool','McAinshLabpooled','Cerespooled'}
            RunDIR=McAinshLabDIRpooled;
        otherwise
            disp('Problem with server name:ABORTING')
            return;
    end
end
            
algstr='EDC';

% Run over treatments
figure;hold on; y=[0 1]; % range
cnt=1;
col={'k','r','b','g','m','c'};

for k=1:length(runtreatments)
	runtreatment=runtreatments{k};
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
if ~isempty(strfind(lower(runtreatment),'tax')) | ~isempty(strfind(lower(runtreatment),'taxol'))
treatnam='Taxol15min';
 else
   if ~isempty(strfind(lower(runtreatment),'noc')) | ~isempty(strfind(lower(runtreatment),'nocod')) | ~isempty(strfind(lower(runtreatment),'nocodozole'))
treatnam='Nocod2hr';
 else
   treatnam='DMSO';
end
end

ExptDat=ExptDats{dataid{1}};
disp(['Loading treatment ' treatnam]);


%						      
% Saving directory for mcmcruns
%
SAVDIR=['Expt_' FileNameSh '_' treatnam '/'];
SAVDIR=[RunDIR '/' SAVDIR];

disp(['Loading rundata for ' SAVDIR]);

load([SAVDIR 'MCMCStats_chain1.mat']);      % Load stats of first chain
load([SAVDIR 'MCMC_EuclDistMargoutput1.mat']); % Load first chain 

paramnams=mcmcparams.variables;
truevalues=mcmcparams.truevalues;
burnin=mcmcparams.burnin;
mcmcrun=mcmcrun(burnin:end,:); % Delete burnin
dat=mcmcparams.datafile;

%
% Plot data and estimated mu
%

% mu posterior
j=1;
[h x]=hist(mcmcrun(:,j),boxs);
plot(x,h/max(h),col{k},'LineWidth',linewidth);
str{cnt}=runtreatments{k};cnt=cnt+1;

% Plot means

plot(tab.means(1)*[1 1],y,[col{k} ':'],'LineWidth',linewidth); % 
str{cnt}='Mean';cnt=cnt+1;

if ~isempty(options)  && isfield(options,'raw') && options.raw==1
if isfield(options,'zfilter') && options.zfilter==0 
  J=1:size(dat,1);
  str='Expt 3D';cnt=cnt+1;
else
    J=find(dat(:,6));
str='Expt 3D (filtered)';cnt=cnt+1;
end
[h x]=hist(dat(J,1),boxs);
%[h x]=hist(dat(:,1),boxs);
plot(x,h/max(h),[col{k} '--'],'LineWidth',linewidth); % 3D raw data

plot(mean(dat(J,1))*[1 1],y,[col{k} '-.'],'LineWidth',linewidth);
str='Mean data';cnt=cnt+1;
% Stats
y=ylim;
%plot(median(dat(J,1))*[1 1],y,'c--');
end

if k==1
SAVDIR1=SAVDIR;
end


end % k

if isempty(options) || ~isfield(options,'legend') || options.legend==1
  legend(str);
end


title('3D measurements and \Delta_{EC}','FontSize',fontsize)
xlabel('3D distance','FontSize',fontsize)
set(gca,'FontSize',fontsize);
if ~isempty(options) && isfield(options,'xmax')
    if strcmp(options.xmax,'auto')==0
    xlim([0 options.xmax]);
    end
end

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR1 'Figures/PFig_CompareTreatmentsPost.fig'],'fig');
  print('-depsc',[SAVDIR1 'Figures/PFig_CompareTreatmentsPost.eps']);
  print('-painters','-dsvg',[SAVDIR1 'Figures/PFig_CompareTreatmentsPost.svg']);
  print('-painters','-dpdf',[SAVDIR1 'Figures/PFig_CompareTreatmentsPost.pdf']);
end


return;
% PLot sd in xy/z

figure;hold on;

% taux posterior
j=2;
[h x]=hist(1./sqrt(mcmcrun(:,j)),boxs);
plot(x,h/(max(h)*(x(2)-x(1))),'r','LineWidth',linewidth);
hm=max(h);

j=3;
[h x]=hist(1./sqrt(mcmcrun(:,j)),boxs);
     plot(x,h/(max(h)*(x(2)-x(1))),'b','LineWidth',linewidth);

y=ylim();
% Plot means
plot(tab.means(4)*[1 1],y,'r:','LineWidth',linewidth);
plot(tab.means(5)*[1 1],y,'b:','LineWidth',linewidth);

if isempty(options) || ~isfield(options,'legend') || options.legend==1
legend('sd_x','sd_z','mean sdxy','mean sdz');
end

title('Spot centre accuracy','FontSize',fontsize)
xlabel('standard deviation, nm','FontSize',fontsize)
set(gca,'FontSize',fontsize);

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Figures/PFig_SpotAccPostEsts.fig'],'fig');
print('-depsc',[SAVDIR 'Figures/PFig_SpotAccPostEsts.eps']);
print('-painters','-dsvg',[SAVDIR 'Figures/PFig_SpotAccPostEsts.svg']);
print('-painters','-dpdf',[SAVDIR 'Figures/PFig_SpotAccPostEsts.pdf']);
end
