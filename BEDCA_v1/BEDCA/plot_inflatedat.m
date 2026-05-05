function plot_inflatedat(dat,SAVFILESTEM)
  % plot_inflatedat(dat,SAVFILESTEM)
  % 
  % Plots the spherical coordinate and intrakinetochore data
  %
  % NJB May 2017

fontsize=20;

figure;
hist(dat(:,1),100);
xlabel('3D distance (nm)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Experimental 3D distance data'); 

if ~isempty(SAVFILESTEM)
  saveas(gcf,[SAVFILESTEM '_3DrData.fig'],'fig');
  print('-depsc',[SAVFILESTEM '_3DrData.eps']);
end

figure; subplot(1,3,1)
hist(dat(:,2),100);
xlim([0 pi])
xlabel('Angle theta','FontSize',fontsize);
set(gca,'FontSize',fontsize);

subplot(1,3,2)
hist(cos(dat(:,2)),100);
xlabel('cos(theta)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Spherical angle histograms');
subplot(1,3,3)
hist(dat(:,3),100);
xlabel('Angle phi','FontSize',fontsize);
set(gca,'FontSize',fontsize);

if ~isempty(SAVFILESTEM)
  saveas(gcf,[SAVFILESTEM '_angularData.fig'],'fig');
  print('-depsc',[SAVFILESTEM '_angularData.eps']);
end


