function dat=simulate_3DdataFromMCMC(mcmcrun,N,deltaz,mcmctheta)
  % dat=simulate_3DdataFromMCMC(mcmcrun,N,deltaz,mcmctheta)
%
  % DONT USE.
% This simulates a 3D data set of spot centres
  % with separation mu, accuracy sigx, sigz
  % that is drawn from the mcmcrun (give post-burnin only)
  %
  % If mcmctheta given also co-samples theta
  %
  % Returns data dat=(r,thetaX,phiX,withinzthreshold,theta,phi)
  % ie both spot position and vec n
  %
  % NJB May 2017

NOT working. NEEDS edits
  
dat=zeros(N,6);

% Determine blocks
blks=ceil(size(mcmcrun,1)/N);
sigxv=1./sqrt(mcmcrun(:,3));sigzv=1./sqrt(mcmcrun(:,4));

for b=1:blks
	for k=1:size(mcmcrun,1)

		mu=mcmcrun(k,1);
		theta=mcmcrun(k,2);
sigx=sigxv(k);sigz=sigzv(k);
		
% random direction. Distribution is \propto sin(theta), so cos(theta) is uniform

	if isempty(thetavalues)
theta=acos(2*rand()-1);  % \vec n 3D distribution.
	else
	  theta=mcmctheta(k);
end

phi=2*pi*rand();

% Simulated position with Gaussian noise. Different in z directions
x=sin(theta)*cos(phi)*mu+normrnd(0,sigx);
y=sin(theta)*sin(phi)*mu+normrnd(0,sigx);
z=cos(theta)*mu+normrnd(0,sigz);

% Compute shift theta, phi of observed position
r=sqrt(x^2+y^2+z^2);
thetaX=acos(z/r);

phiX=atan(y/x)+(x<0)*pi;
%
%Put PhiX in [0 2pi]
if phiX<0 phiX=phiX+2*pi; end

dat(k,:)=[r; thetaX; phiX;  abs(z)<deltaz; theta; phi];
end %k
end

end % function

  
