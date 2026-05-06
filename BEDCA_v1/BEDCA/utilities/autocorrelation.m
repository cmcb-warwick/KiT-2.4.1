function [R n]=autocorrelation(x,m)
% [R n]=autocorrelation(x,m)
%
  % Computes the autocorrelation R(k)=E[x x_{+k}]/var(x) for k=0 (variance),1...m
  % Returns n, number of samples for each k
% 
  % If NaN entries, skips
%

  if isempty(m)
    m=ceil(length(x)/4);
end

  for k=0:m

  xk=x((1+k):end);
  x0=x(1:length(xk));

J=intersect(find(~isnan(x0)),find(~isnan(xk)));


R(k+1)=mean(xk(J).*x0(J))-mean(x0(J))*mean(xk(J));
n(k+1)=length(J);

end %k

R=R/R(1);

end % function




