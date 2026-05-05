function dat=simulate_3DdataFromMCMC_v2(mcmcrun,N,deltaz,mcmctheta)
% dat=simulate_3DdataFromMCMC_v2(mcmcrun,N,deltaz,mcmctheta)
%
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

  dat=[];

% Determine blocks
blks=ceil(N/size(mcmcrun,1));
%disp(['Simulating ' num2str(blks) ' blocks.']);

if ~isempty(mcmctheta)
theta=mcmctheta;
M=size(theta,2);
disp('Using inferred theta from MCMC in simulation of 3D dataset');
 else
   M=100;
disp('Using uniform cos(theta) in simulation of 3D dataset');
end

sigxv=1./sqrt(mcmcrun(:,2));sigzv=1./sqrt(mcmcrun(:,3));
sigxm=sigxv*ones(1,M);sigzm=sigzv*ones(1,M);
muv=mcmcrun(:,1)*ones(1,M);

dat=zeros(blks*size(mcmcrun,1)*M,6);
cnt=1;

for b=1:blks
		
% random direction. Distribution is \propto sin(theta), so cos(theta) is uniform

if isempty(mcmctheta)
  theta=acos(2*rand(size(mcmcrun,1),M)-1);
end
	
phi=2*pi*rand(size(theta));

% Simulated position with Gaussian noise. Different in z directions
x=sin(theta).*cos(phi).*muv+normrnd(0,sigxm);
y=sin(theta).*sin(phi).*muv+normrnd(0,sigxm);
z=cos(theta).*muv+normrnd(0,sigzm);

% Compute shift theta, phi of observed position
r=sqrt(x.^2+y.^2+z.^2);
thetaX=acos(z./r);

phiX=atan(y./x)+(x<0).*pi;
%
%Put PhiX in [0 2pi]
phiX=phiX+2*pi.*(phiX<0);

%size(r),size(thetaX), size(phiX), size(theta), size(phi), size(abs(z(:))<deltaz)
%size([r(:) thetaX(:) phiX(:)  abs(z(:))<deltaz theta(:) phi(:)])

%dat=[dat; [r(:) thetaX(:) phiX(:)  abs(z(:))<deltaz theta(:) phi(:)]];
dat(cnt:(cnt+size(mcmcrun,1)*M-1),:)=[r(:) thetaX(:) phiX(:)  abs(z(:))<deltaz theta(:) phi(:)];
cnt=cnt+size(mcmcrun,1)*M;

end

end % function

  
