function [LastExpFolder] = FindLastExpFolder(AnimalFolder)
% This function is looking for the last experiment's folder on the f: drive
% If no any ratxxxx folder this comes back with rat0000 otherwise this
% come back that folder name which contains the highest number
	curpwd = pwd;

    % The DefDir is 0 if the def_objects folder is missing, 1 otherwise
    LastExpFolder = 'exp000';
    % Change the current path to AnimalFolder
    cd(AnimalFolder);
    % The DirList structure contains the all files' and folders' data of
    % the f: drive
    DirList = dir;
    % Search that folder name which contains the biggest number
    for n=3:size(DirList)
        % Select that names which contains the rat string on the first
        % three positions
        if strncmp(DirList(n).name, 'exp', 3)
            % The MATLAB gives back a vector of zeros and ones when
            % comparing the value of two string. If there is minimum
            % one 1 in the vector the relation is true. The max() gives
            % back the maximum element of the result vector and if it 1
            % the comparison is true
            if max(DirList(n).name > LastExpFolder) == 1
                LastExpFolder = DirList(n).name;
            end
        end
    end
    cd (curpwd);
end
