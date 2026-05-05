function [params typ] =get_priorparams(priors,nam)
% [params typ]=get_priorparams(priors,nam)
%
% From priors structure, get parameters
% RTN=[params dist_type]
%
% eg
% priors.variable={'L','pbs'};
% priors.types={'Gaussian','Beta'};  Gaussian mean variance,...  Gamma power decaycoef 
% priors.params={[0.1 100],[4/5 1]};
%
% NJB. Sep 2012
%

J=find(ismember(priors.variable,nam));

if isempty(J)
  disp(['Failed to find parameter ' nam ' in prior distribution structure'])
priors
params=[];typ='Failed';
return
end


params=priors.params{J};
typ=priors.types{J};

end % function


