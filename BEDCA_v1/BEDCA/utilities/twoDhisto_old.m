function twoDhisto(X,boxes,ranges)
% Draws a 2 D histogram
% X is array of data, and ranges,boxes is [boxx,boxy], [lowx highx lowy highy]
%



counts=[];

dx=(highx-lowx)/boxes(1);
dy=(highy-lowy)/boxes(2);


[n,m]=size(X);
counts=zero(boxes(1),boxes(2));

for k=1:m
nx=floor((X(1,k)-lowx)/dx);
ny=floor((X(2,k)-lowy)/dy);

if nx>=0 & nx<=boxes(1) 
if ny>=0 & ny<=boxes(2)
counts(nx,ny)=counts(nx,ny)+1;
end
end

end


figure;imshow(counts,[]);colormap(jet);colorbar;



