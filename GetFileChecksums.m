%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
[file, path] = uigetfile('','Load Data', 'MultiSelect', 'on');
if path == 0
    return;
end
if ~iscell(file)
    if file ~= 0
        data = [path file];
        filechecksum = Simulink.getFileChecksum(data)
    else
        return;
    end
else
    filechecksum = cell(numel(file));
    for i = 1:numel(file)
        data = [path file{i}];
        filechecksum{i} = Simulink.getFileChecksum(data)
    end
end