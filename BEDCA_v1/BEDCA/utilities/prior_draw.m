function X=prior_draw(priors,nam)
% X=prior_draw(priors,nam)
%
% This returns a draw from teh prior named nam
%
% eg
% priors.variable={'L','pbs'};
% priors.types={'Gaussian','Beta'};
% priors.params={[0.1 100],[4/5 1]};
%
% Gaussian mean, variance
% Gamma  index, rate (not scale)
%
% NJB. Sep 2012
%


J=find(ismember(priors.variable,nam));

if isempty(J)
  disp(['Failed to find parameter ' nam ' in prior distribution structure'])
priors
X=[];
return
end

params=priors.params{J(1)};

switch priors.types{J(1)}

case 'Gaussian'

X=sqrt(params(2))*randn()+params(1);

case 'Gamma'

X=gamrnd(params(1),1)/params(2);

case 'Beta'

X=betarnd(params(1),params(2));

case 'Exp'

X=gamrnd(1,1)/params(1);

otherwise

disp(['No such distribution name ' priors.types{J(1)}])

  X=[];
end

end % function
