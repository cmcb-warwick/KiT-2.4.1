function [datsm datsp nams]=plot_data(DIR,FileName)
  % [datsm datsp nams]=plot_data(DIR,FileName)
  %
  % This plots the 3D distance data in the KT experiments. Plots microscope and plate data and saves in Expt Directory under Figures
  % Returns the x,y,z coordinate data
  %
  % NJB 2016


dats=[];
nams={};
[status msg]=mkdir([DIR '/Figures']);
SAVFILESTEM=[DIR '/Figures/']; 

fontsize=20;

A=load([DIR '/' FileName]);
disp(['Loading data file ' FileName]);

fstr=fields(A);

nams=fieldnames(A);
disp(['Loading all ' num2str(length(nams)) ' treatments.']);
nams

for k=1:length(nams)

if ~isempty(strfind(lower(nams{k}),'tax'))
treatments{k}='Taxol';
 else
   if ~isempty(strfind(lower(nams{k}),'noc'))
treatments{k}='Nocod';
 else
   treatments{k}='DMSO';
end
end

ff=getfield(A,nams{k});

% Microscope coords
dset=ff.microscope.raw.delta;
x=dset.x.all(:);  y=dset.y.all(:);  z=dset.z.all(:);  % Places cols on top of each other

datsm{k}=[x y z];

% Plots
FileName=FileName(1:(end-4));
figure;
plot3(x,y,z,'.');
xlabel('x','FontSize',fontsize);
ylabel('y','FontSize',fontsize);
zlabel('z','FontSize',fontsize);
title(['Microscope coords ' treatments{k}]); 

saveas(gcf,[SAVFILESTEM FileName '_Data3DscatterMS_' treatments{k} '.fig'],'fig');
print('-depsc',[SAVFILESTEM FileName '_Data3DscatterMS_' treatments{k} '.eps']);


% Plate coords
dset=ff.plate.raw.delta;
x=dset.x.all(:);  y=dset.y.all(:);  z=dset.z.all(:);  % Places cols on top of each other

datsp{k}=[x y z];

figure;
plot3(x,y,z,'.');
xlabel('x','FontSize',fontsize);
ylabel('y','FontSize',fontsize);
zlabel('z','FontSize',fontsize);
title(['Plate coords ' treatments{k}]); 

saveas(gcf,[SAVFILESTEM FileName '_Data3DscatterPl_' treatments{k} '.fig'],'fig');
print('-depsc',[SAVFILESTEM FileName '_Data3DscatterPl_' treatments{k} '.eps']);


end %k


