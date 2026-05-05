function [mu sd n]=combinesample_stats(mus,sds,ns,typ)
  % [mu sd n]=combinesample_stats(mus,sds,ns)
  %
  % typ= 'SEM'
  % This takes the mean, SEM data (mu_i, sd_i) for sample size n_i
  % and returns mean, SEM data for combined sample
  %
  % typ='SD' or [] (default
%
% This takes the mean, (sample) SD data (mu_i, sd_i) for sample size n_i
  % and returns mean, SD data for combined sample
  %
  % Give a vector of means, SEMs (or SD) and sample sizes.
  %
  % NOTE: SEM assumed to be mu_i mean of sample i, unbiased variance sd_i^2 = (sample variance)/n_i, unbiased sample variance
  %
% NJB June 2019					

n=sum(ns);
mu=sum(ns.*mus)/n;

if isempty(typ) | strcmp(typ,'SD')

   vars=sds.^2; % Sample variances

var=(n/(n-1))*(sum((ns-1).*vars+ns.*mus.^2)/n - mu^2);
sd=sqrt(var);

		    else
vars=sds.^2.*ns; % Sample variances

var=(n/(n-1))*(sum((ns-1).*vars+ns.*mus.^2)/n - mu^2);
sd=sqrt(var/n);
end
