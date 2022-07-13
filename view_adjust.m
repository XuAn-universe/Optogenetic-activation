function view_adjust(hObject, eventdata, axis)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
view = get(gca, 'View');
set(axis, 'View', view);