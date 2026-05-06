function treatnam=gettreatname(runtreatment,TreatmentLib)
  % treatnam=gettreatname(runtreatment,TreatmentLib)
  %
  % Give runtreatment (analysis treatments, one of (AKA) treatments in TreatmentLib, see ExptData_List)
  %
  % NJB Jan 2020

  fnd=[];
  for k=1:length(TreatmentLib)
	  aka=TreatmentLib{k};

for j=1:length(aka)
	if ~isempty(strfind(lower(runtreatment),lower(aka{j})))

	if isempty(fnd)
	fnd=[k j];
 else

   if length(TreatmentLib{k}{j})> length(TreatmentLib{fnd(1)}{fnd(2)})
fnd=[k j];
end
   
   end

break;
end


end %j
end %k

if ~isempty(fnd)
treatnam=TreatmentLib{fnd(1)}{1};
 else
   treatnam=[];
end
