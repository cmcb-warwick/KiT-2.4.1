

% ktaux+0.5*(sumSqrpxy+sum((mu*si).^2-2*mu*rpxy.*si.*cphi))

si = sin(dat(:,5));  % True values
sXi = sin(dat(:,2));
ri=dat(:,1);

k=0.5*sum((params.mu*si).^2+(ri.*sXi).^2-2*params.mu*si.*ri.*sXi.*cos(dphi));
kapp=0.5*sum((params.mu*si-ri.*sXi).^2);

a=N-3/2;

a/k
