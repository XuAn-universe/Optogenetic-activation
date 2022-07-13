function [LastTrial] = FindLastTrial(DataFolder)
% This function is looking for the last trial
% If no any trial this comes back with trial000 otherwise this
% come back that trial name which contains the highest number
	curpwd = pwd;

    LastTrial = 'trial000.mat';
    % Change the current path to DataFolder
    cd(DataFolder);
    % The DirList structure contains all files' and folders' data
    DirList = dir;
    % Search that trial name which contains the biggest number
    for n = 3:size(DirList)
        % Select that names which contains the trial string on the first
        % five positions
        if strncmp(DirList(n).name, 'trial', 5)
            % The MATLAB gives back a vector of zeros and ones when
            % comparing the value of two string. If there is minimum
            % one 1 in the vector the relation is true. The max() gives
            % back the maximum element of the result vector and if it 1
            % the comparison is true
            if max(DirList(n).name > LastTrial) == 1
                LastTrial = DirList(n).name;
            end
        end
    end
    cd (curpwd);
end