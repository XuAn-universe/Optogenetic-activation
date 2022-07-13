%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
maxvalue = 0.6;
d_absmax = 0.8;
im = get(gco, 'CData');
im = rgb2hsv(im);
im(:, :, 1) = im(:, :, 1)/0.67;
im(:, :, 1) = im(:, :, 1)*2*maxvalue-maxvalue;
im(:, :, 1) = mat2gray(im(:, :, 1), [-d_absmax d_absmax])*0.67;
im = hsv2rgb(im);
set(gco, 'CData', im);