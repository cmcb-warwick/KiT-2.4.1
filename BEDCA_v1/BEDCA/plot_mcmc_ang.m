function [poormovers mvs]=plot_mcmc_ang(mcmcrunphi,truevalues,burnin,subsamp,options)
  % plot_mcmc_ang(mcmcrunphi,truevalues,burnin,subsamp,options)
  %
  % Analyses the angular coordinates
  % subsamp used to plot a subsample of angular variables
  %
  % Give options with 
  % options.savedir options.filename to save figures.
  % options.figarrangement = 'singlefig' (default) produces one figure, 6 panels
  % options.figarrangement = 'multfig' generates a fig per plot
  %
  % NJB

  if isempty(options) | ~isfield(options,'figarrangement')
  singlefig=1; % Plot single figures.
  else
    singlefig=options.figarrangement;
end

if isempty(truevalues)
nrows=1;
 else
   nrows=2;
end

figure;
if singlefig
subplot(nrows,3,1); end
plot(burnin:size(mcmcrunphi,1),mcmcrunphi(burnin:end,1:subsamp:end));
xlabel('MCMC step');

% Number time it moves
J=abs(diff(mcmcrunphi(burnin:end,:),1,1))>10^-10;
mvs=sum(J,1)/(size(mcmcrunphi,1)-burnin+1);

poormovers=find(mvs<0.05);

means=mean(mcmcrunphi(burnin:end,:),1);
stds=std(mcmcrunphi(burnin:end,:),1,1);
[Y I]=sort(means);

if singlefig subplot(nrows,3,2); else figure; end
plot(mvs(I),'k*');
xlabel('Sorted angles');
ylabel('Perc steps moved');

if singlefig subplot(nrows,3,3); else figure; end
hold on;
plot(means(I),'k*');
J=find(stds(I)>0 & isreal(stds(I)) & isreal(means(I)));
if ~isempty(J)
%plot(J,means(I(J)),'*g')
errorbar(J,means(I(J)),stds(I(J)),'*k');
end

if ~isempty(truevalues)
plot(truevalues(I),'r.')
end
xlabel('Sorted angles');

if ~isempty(truevalues)
if singlefig subplot(2,3,4); else figure; end
hist(means'-truevalues,100); %'
xlabel('Difference mean from true value');

if singlefig subplot(2,3,5); else figure; end
hist((means'-truevalues)./stds',100);
xlabel('Difference z');
disp(['Proportion outside z=1.5 (sds) = ' num2str(sum(abs((means'-truevalues)./stds')>1.5)/length(means))]);

if singlefig subplot(2,3,6); else figure; end
z=(means'-truevalues)./stds';[Y I]=sort(z);
plot(z(I),(1:length(z))/length(z));
xlabel('Difference z');
title('Cumulative sum')
disp(['Proportion outside z=1.5 (sds) = ' num2str(sum(abs((means'-truevalues)./stds')>1.5)/length(means))]);
end

if ~isempty(options) & isfield(options,'savedir') & ~isempty(options.savedir) 

FileNam=[options.filename 'Figs/'];
SAVDIR=[options.savedir '/' FileNam];
[status msg]=mkdir(SAVDIR);

saveas(gcf,[SAVDIR '/HiddenVar_' options.param '.fig'],'fig');
  print('-depsc',[SAVDIR '/HiddenVar_' options.param '.eps']);
end

