function unpack_dataGrp(DIR,filenam,prestem)
  %
  %
  % This unpacks a .mat file with multuiple data sets
  %
  % NJ 2020



  hom=pwd;
cd(DIR)

A=load(filenam);
load(filenam); % load into workspace for save
fn=fieldnames(A);

for k=1:length(fn);
save([prestem fn{k}],fn{k});
end
  



cd(hom)
