function str=warningsconvergence_BEDCA(gelrubin)
  % str=warningsconvergence_BEDCA(gelrubin)
  % examines the convergence diagnostics and issues warnings
  % BEDCA MCMC
  % Give gelman rubin stataistics structure
  %
  % NJB 2020

  
  GRthresh=1.02;
  GRthreshvec=1.02;
str=[];

J=find(gelrubin.Rc>GRthresh);

if ~isempty(J)
disp(' ');

for j=1:length(J)
	disp(['CONVERGENCE WARNING: variable ' gelrubin.varnames{J(j)} ' has corrected GR statistic ' num2str(gelrubin.Rc(J(j))) '. Should be below ' num2str(GRthresh)]);
str=[str [' CONVERGENCE WARNING: variable ' gelrubin.varnames{J(j)} ' has corrected GR statistic ' num2str(gelrubin.Rc(J(j))) '. Should be below ' num2str(GRthresh)]];
disp(['Effective sample dimension is ' num2str(gelrubin.d(J(j)))])
end

disp(['Redo runs with more nsteps']);
disp(' ');
end

  if gelrubin.MPSRF>GRthreshvec
  disp(' ');
  disp(['Multi dimension GR stat is ' num2str(gelrubin.MPSRF) '. Advise rerunning with more nsteps.']);
disp(' ');
str=[str ['Multi dimension GR stat is ' num2str(gelrubin.MPSRF) '. Advise rerunning with more nsteps.']];
end
