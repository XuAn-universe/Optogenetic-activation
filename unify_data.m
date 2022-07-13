function unify_data(datadir)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
if ~exist(datadir, 'dir')
    disp('Directory doesn''t exist');
    return;
end
LastExpFolder = FindLastExpFolder(datadir);
ExpNumber = str2double(LastExpFolder(4:6));
n = 0;
for i = numel(datadir):-1:1
    if datadir(i) == '\'
        n = n+1;
        loc(n) = i;
        if n == 2
            break;
        end
    end
end
stidir = [datadir(1:loc(2)) 'Stimulus Parameter' datadir(loc(1):end)];
for i = 1:ExpNumber
    cam1dir = [datadir '\exp' num2str(i, '%03d')];
    cam2dir = [datadir(1:end-1) ' PG2)\exp' num2str(i, '%03d')];
    if exist(cam1dir, 'dir') && exist(cam2dir, 'dir')
        movefile([cam2dir '\*'], cam1dir);
    else
        disp(['''' cam1dir ''' or ''' cam2dir ''' doesn''t exist']);
    end
    if exist(stidir, 'dir')
        copyfile([stidir '\exp' num2str(i, '%03d') '.mat'], cam1dir);
    else
        disp(['''' stidir ''' doesn''t exist']);
    end
end
disp('Done !');