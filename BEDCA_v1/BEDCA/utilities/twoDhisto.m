function [Fhs1,Fhs2,vcounts]=twoDhisto(X,boxes,ranges,str,Fgs)%,fitflag,fitthresh)
% [Fhs1,Fhs2,vcounts]=twoDhisto(X,boxes,ranges,str,Fgs)
% Draws a 2 D histogram
% X is array of data (x,y), and boxes, ranges is [boxx,boxy], [lowx highx lowy highy]
% X(1,:) is x, X(2,:) is y. If X(3,:) given computes mean image.
   %
   % str is printed on figures. Can take an array of figure handles Fgs.
%
   %
% returns figure handles


counts=[]; % Counts in 2D lattice

lowx=ranges(1);
lowy=ranges(3);

dx=(ranges(2)-lowx)/boxes(1);
dy=(ranges(4)-lowy)/boxes(2);


[n,m]=size(X);
counts=zeros(boxes(1),boxes(2));

vcounts=zeros(boxes(1),boxes(2));

meanx=zeros(1,boxes(1));
meany=zeros(1,boxes(2));
varx=zeros(1,boxes(1));;
vary=zeros(1,boxes(2));;


for k=1:m
nx=ceil((X(1,k)-lowx)/dx);
ny=ceil((X(2,k)-lowy)/dy);

if nx>0 & nx<=boxes(1) 
if ny>0 & ny<=boxes(2)
counts(nx,ny)=counts(nx,ny)+1;
if size(X,1)>2
vcounts(nx,ny)=vcounts(nx,ny)+X(3:end,k);
end

meanx(nx)=meanx(nx)+X(2,k);
meany(ny)=meany(ny)+X(1,k);
varx(nx)=varx(nx)+X(2,k)*X(2,k);
vary(ny)=vary(ny)+X(1,k)*X(1,k);

end
end
end

meanx=meanx./sum(counts,2)'; %'
meany=meany./sum(counts,1);

varx=(varx./sum(counts,2)'-meanx.*meanx)./(sum(counts,2)'-1);
vary=(vary./sum(counts,1)-meany.*meany)./(sum(counts,1)-1);

if isempty(Fgs)
Fhs1=figure;subplot(1,2,1)
 else
   Fhs1=figure(Fgs(1));
end
%Fhs1=figure;%subplot(2,1,1);
%ploterrors(lowx+dx*(1:boxes(1)),meanx,sqrt(varx),1.65,'k');

plot(lowx+dx*(1:boxes(1)),meanx,'k',lowx+dx*(1:boxes(1)),meanx+1.65*sqrt(varx),'s',lowx+dx*(1:boxes(1)),meanx-1.65*sqrt(varx),'s')
title([str ': mean y of sampling position']);
ylabel('Mean');
xlabel('x coord');

if isempty(Fgs)
%Fhs2=figure;
Fhs2=[];
subplot(1,2,2)
 else
   Fhs2=figure(Fgs(2));
end
%Fhs2=figure;
%subplot(2,1,2);

%ploterrors(lowy+dy*(1:boxes(1)),meany,sqrt(vary),1.65,'k');

plot(lowy+dy*(1:boxes(2)),meany,'k',lowy+dy*(1:boxes(2)),meany+1.65*sqrt(vary),'s',lowy+dy*(1:boxes(2)),meany-1.65*sqrt(vary),'s')
     title([str ': mean x of sampling position']);
ylabel('Mean');
xlabel('y coord');

if 0
{figure;
plot(lowy+dy*(1:boxes(2)),meany,'k',lowy+dy*(1:boxes(2)),meany+1.65*sqrt(vary),'s',lowy+dy*(1:boxes(2)),meany-1.65*sqrt(vary),'s');
xlabel('Estimated separation (arbitrary units)');ylabel('Test statistic');
}
end

					   
%figure;imshow(counts,[],'InitialMagnification','fit');colormap(jet);colorbar;
% This flips the axes
[n,m]=size(counts);
counts=counts';%'
figure;
%subplot(1,2,2);
imshow(max(-1,log(counts(m:-1:1,1:n))),[],'InitialMagnification','fit');
colormap(jet);colorbar;
title(['Log counts']);
xlabel('x');ylabel('y')

figure
%subplot(1,2,1);
imshow(counts(m:-1:1,1:n),[],'InitialMagnification','fit');
colormap(jet);colorbar;
title(['Position counts'])
xlabel('x');ylabel('y')

if size(X,1)>2
J=find(counts>0);
vcounts(J)=vcounts(J)./counts(J);
zs=ones(size(vcounts));zs(J)=0;
I(:,:,1)=vcounts/max(vcounts(:));
I(:,:,3)=vcounts.*(1-zs)/max(vcounts(:))+zs;
I(:,:,2)=vcounts/max(vcounts(:));

figure;

subplot(2,1,1)
imshow(I,[],'InitialMagnification','fit'); %colorbar();colormap(jet)
title(['Average ' str])

subplot(2,1,2)
imshow(vcounts,[],'InitialMagnification','fit'); colorbar();colormap(jet)

 else
   vcounts=[];
end
