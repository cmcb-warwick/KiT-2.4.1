function [imgCrpd,coords] = indivSpot(img,coords,opts)

% get pixel resolution
imageSize = size(img);
imgw = opts.imgHalfWidth;
bgcol = opts.bgcol;

% predefine cropped image
imgCrpd = ones(2*imgw+1)*bgcol;

% calculate centre pixel
centrePxl = round(coords);

%CCC test 2025.03.26
maxImg = max(img, [], 3);
irange=stretchlim(maxImg,[0.1 0.9995]);
%test end

% max project over various z-slices around point - CCC edit 2023.10.18
if opts.zProject == 1 %for full image
    img = max(img(:,:,max(1,centrePxl(3)-4):min(centrePxl(3)+4,opts.imageSize(3))), [], 3);
elseif opts.zProject == -1 %for zoomed image when creating sS structure
    img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);
else
    img = max(img(:,:,max(1,centrePxl(3)-1):min(centrePxl(3)+1,opts.imageSize(3))), [], 3);
end
%img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);

% produce cropped image around track centre
% %CCC test 2023.07.17
% if centrePxl(1)-imgw+1 < 1
%     centrePxl(1) = imgw;
% elseif centrePxl(1)+imgw+1 > imageSize(2)
%     centrePxl(1) = imageSize(2)-imgw-1;
% end
% if centrePxl(2)-imgw+1 < 1
%     centrePxl(2) = imgw;
% elseif centrePxl(2)+imgw+1 > imageSize(1)
%     centrePxl(2) = imageSize(1)-imgw-1;
% end
% %CCC test end 2023.07.17 %this works how I want it to
% xReg = [max(centrePxl(1)-imgw+1,1) min(centrePxl(1)+imgw+1,imageSize(2))];
% yReg = [max(centrePxl(2)-imgw+1,1) min(centrePxl(2)+imgw+1,imageSize(1))];
% imgCrpd(1:diff(yReg)+1,1:diff(xReg)+1) = img(yReg(1):yReg(2),xReg(1):xReg(2));


%CCC test 2023.12.05
if centrePxl(1)-imgw < 1
    centrePxl(1) = imgw+1;
elseif centrePxl(1)+imgw > imageSize(2)
    centrePxl(1) = imageSize(2)-imgw-1;
end
if centrePxl(2)-imgw < 1
    centrePxl(2) = imgw+1;
elseif centrePxl(2)+imgw > imageSize(1)
    centrePxl(2) = imageSize(1)-imgw-1;
end

xReg = [max(centrePxl(1)-imgw,1) min(centrePxl(1)+imgw,imageSize(2))];
yReg = [max(centrePxl(2)-imgw,1) min(centrePxl(2)+imgw,imageSize(1))];
imgCrpd(1:diff(yReg)+1,1:diff(xReg)+1) = img(yReg(1):yReg(2),xReg(1):xReg(2));
%CCC test end 2023.12.05 %works how I want it to
%CCC test 2025.03.26

% define contrast stretch and apply
if opts.zProject == -1
    irange=stretchlim(imgCrpd,[0.1 0.9995]); %CCC 2025.03.26
end
imgCrpd = imadjust(imgCrpd, irange, []);

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - [xReg(1) yReg(1)];

end

