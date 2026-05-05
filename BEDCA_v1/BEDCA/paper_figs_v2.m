function paper_figs_v2(exptnum,runtreatment,nam,server,options)
% paper_figs_v2(exptnum,runtreatment,nam,server,options)
%
% Give specific directory in MCMCruns or experiment number/treatment (untreated,DMSO,tax,nocod)
% eg paper_figs([],[],'Expt_iM_180903_Ndc80GFPC_DMSO_selectedSpots9G3_Ndc80GFP_EDC_DMSO',[])
% or eg paper_figs(6,'DMSO',[],[])
%
  % Note: Untreated relates to original Kit file. Relabelled DMSO for MCMC rundata
  %
% Plots 3D data and posterior with means, and short bars showing sd.
% option.add1D adds the 1D data.
% option.oneDfile specifies file to use (added for DMSO/untreated issue)
%
% options to control plot
% options.boxs options.linewidth options.xmax (maximum of x in plot,
% default 150, also allows options.xmax='auto' for autoscaling)
% options.legend == 0 turns off legends. On by default.
%
% This plots figures for paper and saves in Figures in aboce directory
%
% NJB april 2019


NOT FINISHED. NEEDS 1D data loaded.
  
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

if ~isempty(exptnum) && exptnum>0
dataid={exptnum,runtreatment};  % First is id below in FileNames{}, Second category is 'tax'/'tax15min', 'noc'/'noc2hr', 'DMSO', 'untreated', 'all'(Not done yet).  Load extra categories associated with new treatments in loadDistData.m

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
treatnam=gettreatname(runtreatment,TreatmentLib);

if isempty(treatnam)
disp(['Failed to find treatment ' runtreatment ' in TreatmentLib']);
disp('Check name against TreatmentLib in ExptData_List');
disp('ABORTING');
return;
end


ExptDat=ExptDats{dataid{1}};
disp(['Loading treatment ' treatnam]);

if ~isempty(options) & isfield(options,'add1D')

if ~isempty(options) & isfield(options,'oneDfile')
A=load([ExptDataLoc options.oneDfile]);
else
A=load([ExptDataLoc FileNames{dataid{1}}]);
end
S=fields(A);
idfield=[];

switch lower(runtreatment)
case 'untreated'

for j=1:length(S)
if ~isempty(findstr('untreated',S{j})) idfield=j; break;
end
end

case 'dmso'

for j=1:length(S)
if ~isempty(findstr('DMSO',S{j})) idfield=j; break;
end
end

case {'tax','taxol'}

for j=1:length(S)
if ~isempty(findstr('tax',S{j})) idfield=j; break;
end
end

case {'noc','nocodazol','nocod'}

for j=1:length(S)
if ~isempty(findstr('noc2hr',S{j})) idfield=j; break;
end
end

otherwise
disp(['Trouble locating 1D data for treatment ' runtreatment])

end

if isempty(idfield)
disp(['Not assigned: Trouble locating 1D data for treatment ' runtreatment '. ABORTING'])
return;
end

if ~isempty(options) & isfield(options,'oneDfile')
disp(['Loading 1D data from ' ExptDataLoc options.oneDfile ', field ' S{idfield}]);
else
disp(['Loading 1D data from ' ExptDataLoc FileNames{dataid{1}} ', field ' S{idfield}]);
end

kitdat=getfield(A,S{idfield});

end

%
% Saving directory for mcmcruns
%
SAVDIR=['Expt_' FileNameSh '_' treatnam '/'];
SAVDIR=[RunDIR '/' SAVDIR];

 else
%						      
% Saving directory for mcmcruns, given directory in MCMCruns and given filename nam
%

   SAVDIR=[RunDIR '/' nam '/'];
end

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
figure;hold on;
if ~isempty(options) && isfield(options,'zfilter') && options.zfilter==0
  J=1:size(dat,1);
str{1}='Expt 3D';
else
    J=find(dat(:,6));
str{1}='Expt 3D (filtered)';
end
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

% Plot 1D
if ~isempty(options) & isfield(options,'add1D')
[h x]=hist(1000*kitdat.plate.raw.delta.oneD,boxs);
plot(x,h*hm/max(h),'b','LineWidth',linewidth); % 1D raw in microns

plot(1000*mean(kitdat.plate.raw.delta.oneD,'omitnan')*[1 1],y,'b:','LineWidth',linewidth);
plot(1000*(mean(kitdat.plate.raw.delta.oneD,'omitnan')+std(kitdat.plate.raw.delta.oneD,'omitnan')*[1 -1]/sqrt(sum(isnan(kitdat.plate.raw.delta.oneD)))),0.975*y(2)*[1 1],'b','LineWidth',2*linewidth);

str=[str {'\Delta_{ED}','\Delta_{1D}'}]; %,'mean 1D'}];
 else
   str=[str {'\Delta_{ED}'}]; %,'mean measured','mean EC'}];
end

% Plot means
plot(mean(dat(J,1))*[1 1],y,'k:','LineWidth',linewidth);

plot(mean(dat(J,1))+std(dat(J,1))*[1 -1]/sqrt(length(J)),0.975*y(2)*[1 1],'k','LineWidth',2*linewidth);
plot(tab.means(1)*[1 1],y,'r:','LineWidth',linewidth);
plot(tab.means(1)+tab.sd(1)*[1 -1],0.975*y(2)*[1 1],'r','LineWidth',2*linewidth);


if isempty(options) || ~isfield(options,'legend') || options.legend==1
  legend(str)
end

title('3D measurements and \Delta_{EC}','FontSize',fontsize)
xlabel('3D distance','FontSize',fontsize)
set(gca,'FontSize',fontsize);
if ~isempty(options) && isfield(options,'xmax')
    if strcmp(options.xmax,'auto')==0
    xlim([0 options.xmax]);
    end
else
  xlim([0 150]);
end

if ~isempty(SAVDIR)
  saveas(gcf,[SAVDIR 'Figures/PFig_DataWithPostEsts.fig'],'fig');
  print('-depsc',[SAVDIR 'Figures/PFig_DataWithPostEsts.eps']);
  print('-painters','-dsvg',[SAVDIR 'Figures/PFig_DataWithPostEsts.svg']);
  print('-painters','-dpdf',[SAVDIR 'Figures/PFig_DataWithPostEsts.pdf']);
end

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
