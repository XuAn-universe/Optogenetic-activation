%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
foldername = {'I:\Videos for head fixed stimulation DLC tracking\20160716_Laser(Tcerg1Ai32 TSP4)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160716_Laser(Tcerg1Ai32 TSP4)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160717_Laser(Tcerg1Ai32 TSP5)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160717_Laser(Tcerg1Ai32 TSP5)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160717_Laser(Tcerg1Ai32 TSP5)\exp002\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160717_Laser(Tcerg1Ai32 TSP5)\exp002\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160523_Laser(Sema3EAi32 TSP1)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160523_Laser(Sema3EAi32 TSP1)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160523_Laser(Sema3EAi32 TSP2)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160523_Laser(Sema3EAi32 TSP2)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160831_Laser(Sema3EAi32 TSP3)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20160831_Laser(Sema3EAi32 TSP3)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20161006_Laser(Sema3EAi32 TSP4)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20161006_Laser(Sema3EAi32 TSP4)\exp001\PG2_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20161006_Laser(Sema3EAi32 TSP5)\exp001\PG1_videos',...
    'I:\Videos for head fixed stimulation DLC tracking\20161006_Laser(Sema3EAi32 TSP5)\exp001\PG2_videos'};
curpwd = pwd;
workbar(0, 'Computing Ongoing...', 'Progress');
for ii = 1:numel(foldername)
    cd(foldername{ii});
    DirList = dir;
    totalfile = numel(DirList);
    for i = 3:totalfile
        if strcmp(DirList(i).name(end-3:end), '.csv') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            [num, txt, raw] = xlsread(DirList(i).name);
            leftpaw.columnx = [];
            leftpaw.columny = [];
            leftpaw.columnlikelihood = [];
            rightpaw.columnx = [];
            rightpaw.columny = [];
            rightpaw.columnlikelihood = [];
            jaw.columnx = [];
            jaw.columny = [];
            jaw.columnlikelihood = [];
            for j = 1:size(txt, 1)
                for k = 1:size(txt, 2)
                    switch txt{j, k}
                        case 'LeftPaw'
                            if strcmp(txt{j+1, k}, 'x')
                                leftpaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                leftpaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                leftpaw.columnlikelihood = k;
                            end
                        case 'RightPaw'
                            if strcmp(txt{j+1, k}, 'x')
                                rightpaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                rightpaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                rightpaw.columnlikelihood = k;
                            end
                        case 'Jaw'
                            if strcmp(txt{j+1, k}, 'x')
                                jaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                jaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                jaw.columnlikelihood = k;
                            end
                    end
                end
            end
            if ~isempty(leftpaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_trajectory.mat'])
                trajectory = num(:, [leftpaw.columnx leftpaw.columny leftpaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_trajectory.mat'], 'trajectory');
            end
            if ~isempty(rightpaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_rightpaw_trajectory.mat'])
                trajectory = num(:, [rightpaw.columnx rightpaw.columny rightpaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_rightpaw_trajectory.mat'], 'trajectory');
            end
            if ~isempty(jaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_jaw_trajectory.mat'])
                trajectory = num(:, [jaw.columnx jaw.columny jaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_jaw_trajectory.mat'], 'trajectory');
            end
        end
        if strcmp(DirList(i).name(end-2:end), '.h5') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            delete(DirList(i).name);
        end
        if strcmp(DirList(i).name(end-6:end), '.pickle') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            delete(DirList(i).name);
        end
        if strcmp(DirList(i).name(end-3:end), '.mp4') && ~isempty(strfind(DirList(i).name, 'DeepCut')) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_annotated.mp4'])
            movefile(DirList(i).name, [DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_annotated.mp4']);
        end
    end
    workbar(ii/numel(foldername), ['folder' num2str(ii) '/' num2str(numel(foldername))], 'Progress');
    drawnow;
end
cd(curpwd);
msgbox('Done !');

%%
curpwd = pwd;
try
    cd(foldername);
end
foldername = uigetdir('', 'Select a folder to process');
if foldername ~= 0
    cd(foldername);
    DirList = dir;
    workbar(0, 'Computing Ongoing...', 'Progress'); 
    totalfile = numel(DirList);
    for i = 3:totalfile
        if strcmp(DirList(i).name(end-3:end), '.csv') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            [num, txt, raw] = xlsread(DirList(i).name);
            leftpaw.columnx = [];
            leftpaw.columny = [];
            leftpaw.columnlikelihood = [];
            rightpaw.columnx = [];
            rightpaw.columny = [];
            rightpaw.columnlikelihood = [];
            jaw.columnx = [];
            jaw.columny = [];
            jaw.columnlikelihood = [];
            for j = 1:size(txt, 1)
                for k = 1:size(txt, 2)
                    switch txt{j, k}
                        case 'LeftPaw'
                            if strcmp(txt{j+1, k}, 'x')
                                leftpaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                leftpaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                leftpaw.columnlikelihood = k;
                            end
                        case 'RightPaw'
                            if strcmp(txt{j+1, k}, 'x')
                                rightpaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                rightpaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                rightpaw.columnlikelihood = k;
                            end
                        case 'Jaw'
                            if strcmp(txt{j+1, k}, 'x')
                                jaw.columnx = k;
                            end
                            if strcmp(txt{j+1, k}, 'y')
                                jaw.columny = k;
                            end
                            if strcmp(txt{j+1, k}, 'likelihood')
                                jaw.columnlikelihood = k;
                            end
                    end
                end
            end
            if ~isempty(leftpaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_trajectory.mat'])
                trajectory = num(:, [leftpaw.columnx leftpaw.columny leftpaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_trajectory.mat'], 'trajectory');
            end
            if ~isempty(rightpaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_rightpaw_trajectory.mat'])
                trajectory = num(:, [rightpaw.columnx rightpaw.columny rightpaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_rightpaw_trajectory.mat'], 'trajectory');
            end
            if ~isempty(jaw.columnx) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_jaw_trajectory.mat'])
                trajectory = num(:, [jaw.columnx jaw.columny jaw.columnlikelihood]);
                save([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_jaw_trajectory.mat'], 'trajectory');
            end
        end
        if strcmp(DirList(i).name(end-2:end), '.h5') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            delete(DirList(i).name);
        end
        if strcmp(DirList(i).name(end-6:end), '.pickle') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
            delete(DirList(i).name);
        end
        if strcmp(DirList(i).name(end-3:end), '.mp4') && ~isempty(strfind(DirList(i).name, 'DeepCut')) %&& ~isfile([DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_annotated.mp4'])
            movefile(DirList(i).name, [DirList(i).name(1:strfind(DirList(i).name, 'DeepCut')-1) '_annotated.mp4']);
        end
        workbar(i/totalfile, ['file' num2str(i) '/' num2str(totalfile)], 'Progress');
    end
    drawnow;
end
cd(curpwd);
msgbox('Done !');