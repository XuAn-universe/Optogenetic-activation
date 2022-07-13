%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
% add a scale bar to motor mapping ROI
x = 497;
y = 1204;
n = round(1000/2.63);
figure;
im = imrotate(ROI, 90);
im(y:y+8, x:x+n) = 1;
imshow(im, [])

%%
% interpolate a ROI of an image
h = imfreehand('Closed', 'True');
wait(h);
BW = createMask(h);
im = roifill(x__2(:, :, 1), BW);
figure;imshow(im, [])

%%
% add a scale bar to intrinsic imaging ROI
xstart = 32;
ystart = 17;
height = 1024;
width = 1280;
ROI = GreenImage(ystart+1:ystart+height, xstart+1:width+xstart);
% add a scale bar to an image
x = 800;
y = 916;
n = round(1000/2.63);
figure;
im = ROI;
im(y:y+16, x:x+n) = min(im(:));
imshow(im, [])