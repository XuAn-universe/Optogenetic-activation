function im = boost_IM_brightness(im, thr)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
if thr ~= 1
    vidFrame = rgb2hsv(im);
    intensity = vidFrame(:, :, 3);
    intensity = mat2gray(intensity, [0 thr]);
    vidFrame(:, :, 3) = intensity;
    im = hsv2rgb(vidFrame);
end