%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
CP = get(gca, 'CameraPosition');
CT = get(gca, 'CameraTarget');
CUV = get(gca, 'CameraUpVector');
CVA = get(gca, 'CameraViewAngle');
P = get(gca, 'Position');

%%
set(gca, 'CameraPosition', CP, 'CameraTarget', CT, 'CameraUpVector', CUV, 'CameraViewAngle', CVA, 'Position', P);