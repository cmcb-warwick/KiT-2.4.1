function [mcmcdat mcmcdiagnos]=load_mcmc_runs(dirnam,filenam,N,outformat)
% mcmcdat=load_mcmc_runs(dirnam,filenam,N,outformat)
% 
% General routine to load runs/list available. Returns mcmcdat str with .mcmcparams, .mcmcrun etc
% or an array of such data over the multiple runs. Returns convergence diagnostic file if available.
%
% outformat= 'single' (takes first run), 'all' (returns all mcmcdat(k)), TODO 'poolsamples' (returns all samples concatenated). NOTE: On latter burnin must be reset in mcmcparams or some samples will be lost.
%
% N (optional []) is expected number of runs.
%
% For list set filenam=[] and dirnam='MCMCruns'/ or more specific, eg dirnam='MCMCruns/SpDisc2sGSp_May12_Cell1'
%
% mcmcdat is then a cell array of names
%
% TO ADD a search facility
%
% This loads the specified file. If filenam=[] lists the available runs.
%
% Author: Burroughs
% Date: Revised June 2012. V1 July 2013.


if isempty(dirnam)

  mcmcdat=[];mcmcdiagnos=[];
  L=dir();

disp('Available MCMC run dirs are:')
  namset={};

for i=1:length(L)

if strmatch('MCMC',L(i).name) %~isempty(pos)
  disp(L(i).name)
end

end % i

  return;
end


if ~isempty(filenam)
    disp(['Loading: ' dirnam '/' filenam 'output files' ])

if isempty(outformat)
    outformat='single';
end


switch outformat
case 'single'  % loads first run
mcmcdat=load ([dirnam '/' filenam 'output1.mat']);
mcmcdat.format='single';

case 'all'

L=dir(dirnam);
n=0;
for i=1:length(L)
str=L(i).name(1:min(length(filenam),length(L(i).name)));
numb=L(i).name(length(filenam)+7:strfind(L(i).name,'.mat')-1); % This is run number at string.

if strcmp(filenam,str) & ~L(i).isdir  & ~isempty(numb) & ~isempty(str2num(numb)) 
  n=n+1;
end

end %i

if ~isempty(N) & N~=n % Check there are N
  disp(['Found ' num2str(n) ' runs. Given N is incorrect. Loading available files.'])
end
N=n;

for r=1:N
  mr=load ([dirnam '/' filenam 'output' num2str(r) '.mat']);
  mr.format=['run' num2str(r)];
  mcmcdat(r)=mr;
end %r

end % switch

% Load convergence diagnostic if available

try
%[dirnam '/' filenam 'mcmcdiagnos.mat']
load ([dirnam '/' filenam 'convergencediag.mat']);
catch
mcmcdiagnos=[];
end

 else % No filenam. Output a list

f=strfind(dirnam,'/');

if isempty(f) | f(end) <length(dirnam)
  dirnam=[dirnam '/'];
end

f=strfind(dirnam,'/');
if strcmp(dirnam(1:8),'MCMCruns') & length(f) <2 % Listing cells

disp('Available Cell Runs are:');
ls(dirnam)

mcmcdat=ls(dirnam);

  mcmcdiagnos=[];

else

  L=dir(dirnam);cnt=1;nams={};

disp('Available runs are:')
namset={};

for i=1:length(L)

pos=strfind(L(i).name,'output'); 

if ~isempty(pos)
  str=L(i).name(1:(pos(1)-1));

if ~ismember(str,namset)  % This removes redundancy in chains
  namset=union(namset,str);
  nams{cnt}=str;cnt=cnt+1;
%disp(str)
end

end % pos

end %i

mcmcdat=nams;
mcmcdiagnos=[];

% Find number
for j=1:length(nams)
  n=0;
for i=1:length(L)
  str=L(i).name(1:min(length(nams{j}),length(L(i).name)));
numb=L(i).name(length(nams{j})+7:strfind(L(i).name,'.mat')-1); % This is run number at string.
if strcmp(str,nams{j}) & ~L(i).isdir  & ~isempty(numb) & ~isempty(str2num(numb)) 
  n=n+1;
end
end %i

  disp([nams{j} ', ' num2str(n) ' runs']);

end %j


end
end



