function dat=extract_sphericalcoords(dset,zfilterz,SAVFILESTEM)
  %  dat=extract_sphericalcoords(dset,zfilterz,SAVFILESTEM)
  %
  % This takes data dset with paired sister data and returns spherical coordinates for each sister
  %
  % dset.x .y .z  and each has field .all(N x 2), .P, .AP, .N .S
  % Give zfilterz for determining samples satisfying filter criterion 
  %
  % Returns array with cols r, theta, phi (in -pi/2 to pi/2), KTpair, sister, filter
  %
  % Can optionally draw the 3D plots. But better to use plot_data
  %
  % NJB May 2017

zfilterz=zfilterz/1000; % data in microns

x=dset.x.all(:);  y=dset.y.all(:);  z=dset.z.all(:);  % Places cols on top of each other

r=sqrt(x.^2 + y.^2 +z.^2);
N=size(dset.x.all,1);

dat=[r acos(z./r) atan(y./x) [(1:N)';(1:N)'] [ones(N,1);2*ones(N,1)] abs(z)<zfilterz];

%return

if ~isempty(SAVFILESTEM)
fontsize=20;
figure;plot3(r.*cos(dat(:,3)),r.*sin(dat(:,3)),z,'.')
xlabel('x','FontSize',fontsize);
ylabel('y','FontSize',fontsize);
zlabel('z','FontSize',fontsize);

  saveas(gcf,[SAVFILESTEM 'Data3Dscatter.fig'],'fig');
  print('-depsc',[SAVFILESTEM 'Data3Dscatter.eps']);
end

