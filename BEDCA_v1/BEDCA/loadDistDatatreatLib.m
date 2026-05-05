function [dats treatments nams]=loadDistDatatreatLib(FileName,TreatLib,typ,filterpars,SAVFILESTEM)
  % [dats treatments nams]=loadDistDatatreatLib(FileName,TreatLib,typ,filterpars,SAVFILESTEM)
  % 
  % This loads 3D distance coord data for 2 fluorophore vectors (in microns)
  % Scaling microns -> nm is at the end of the function
  %
  % Give the TreatmentLib (as defined in ExptDataList)
  % Uses 2nd string, TreatLIb{k}{2}, in searches
  % 
  % typ is Treatment as given in TreatmentLib (either serach term or aka) 
  %
  % Add any other categories
  %
  % filterpars is zfilterz and  max 3D distance filter, .zfilterz .filterdist
  % If not wanted leave off field or make []
  %
  % Returns a multiple structure dats{k}.3Ddat .treatment
  % .3Ddat is (r theta phi sis-pair-id sis zfilterFlag cell-id KTpair-idKiT)
  % This assumes KTs are paired. Otherwise use loadDistDataOrig which doesnt have these extra columns.
  %
  % Distance in nms
  %
  % NJB July 2017. Revised 2018. Additional added 2019. NOTE: typ choices could be simplified with an array
% Needs revising to lift treatments from TreatmentLib
  
dats=[];
treatments=[];
nams={};
savfile=[];

% Extract filters
if isfield(filterpars,'filterdist') & ~isempty(filterpars.filterdist)
filterdist=filterpars.filterdist;
filteron=1;
 else
   filteron=0;
end

%
% This is not used in algorithm/comparsison only to Smith et al. eLife 2016
%
if isfield(filterpars,'zfilterz') & ~isempty(filterpars.zfilterz)
zfilterz=filterpars.zfilterz;
 else
   zfilterz=Inf;
end  

A=load(FileName);
disp(['Loading data file ' FileName]);

fstr=fields(A);

%
%
% Find search term. Easier to use contains, version 2017 or later
%
  scre=zeros(1,length(TreatLib));
  for k=1:length(TreatLib)
	  for j=2:length(TreatLib{k})
		  scre(k)=scre(k)+strcmp(typ,TreatLib{k}{j});
	end
end
	
H=find(scre>0);	

if length(H)>1
disp('Identity of treatment unclear. Respecify runtreatment ABORT.');
return;
end

if isempty(H)
disp('Cannot find runtreatment in TreatmentLib');
return;
end

searchstr=lower(TreatLib{H}{2});

%
% Find treatment in (data) fields of A
%
J=strfind(lower(fstr),searchstr);
fndstr=zeros(1,length(fstr));
for k=1:length(fstr) fndstr(k)=~isempty(J{k}); end

if sum(fndstr)~=1
  disp(['Error in finding treatment ' typ ', search string ' searchstr ', in ' FileName '. ABORT.']);
disp('Fields are:')
fstr

return;
end

treatments{1}=TreatLib{H}{1};
ff=getfield(A,fstr{find(fndstr)});
if ~isempty(SAVFILESTEM) savfile=[SAVFILESTEM treatments{1}];
 else savfile=[];
end
%below: ccc edit 2026.01.16, hopefully will prevent errors downstream at
%diversitymeasure_v2. Changed str2num(ff.label(:,3:4)) to
%str2num(ff.label(:,1:4)). This will hopefully account for experiment
%pooling in terms of the cell counting (I don't think this is used anywhere
%else in BEDCA so am pretty confident I can change this here without too
%much drama.
dats{1}=[extract_sphericalcoords(ff.microscope.raw.delta,zfilterz,savfile) [str2num(ff.label(:,1:4));str2num(ff.label(:,1:4))] [str2num(ff.label(:,5:7));str2num(ff.label(:,5:7))]];

%
% Changes microns to nm
%
for k=1:length(dats)
	dats{k}(:,1)=dats{k}(:,1)*1000;

if filteron
J=find(dats{k}(:,1)<=filterdist);
dats{k}=dats{k}(J,:);
end

end



