%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
exp_path = {'I:\Videos for head fixed stimulation DLC tracking\20161006_Laser(Sema3EAi32 TSP5)\exp001'};
pre = 0.5;
exc_path = 0;
excID = [0:381 1891];
FrameRate = 100;
nbins = 2000;
body_check = 3; % 1 left forelimb, 2 right forelimb, 3 jaw
data_path_front = cell(0);
data_path_side = cell(0);
speedf = [];
speeds = [];
accelerationf = [];
accelerations = [];
n = 1;
curpwd = pwd;
for i = 1:numel(exp_path)
    cd(exp_path{i});
    DirList = dir;
    totalfile = numel(DirList);
    workbar(0, 'Computing Ongoing...', 'Progress'); 
    for j = 3:totalfile
        if body_check == 3
            if numel(DirList(j).name) <= 17
                continue;
            end
            if strcmp(DirList(j).name(1:9), 'PG2camera') && strcmp(DirList(j).name(end-17:end), 'jaw_trajectory.mat')
                temp = load(DirList(j).name);
                trajectory = temp.trajectory;
                trialID = DirList(j).name(10:end-19);
                if any(excID == str2double(trialID)) && exc_path == i
                    if exist('jaw_exclude', 'dir') ~= 7
                        mkdir jaw_exclude;
                    end
                    movefile(DirList(j).name, 'jaw_exclude');
                    continue;
                end
                csvdata = load(['PG2gpio' trialID '.csv']);
                if any(csvdata(:, 1))
                    frame_onsite = find(csvdata(:, 1), 1);
                    trajectory_pre = trajectory(frame_onsite-pre*FrameRate:frame_onsite-1, 1:2);
                    speed = sqrt(diff(trajectory_pre(:, 1)).^2+diff(trajectory_pre(:, 2)).^2);
                    acceleration = abs(diff(speed));
                    speedf(n) = sum(speed);
                    accelerationf(n) = sum(acceleration);
                else
                    continue;
                end
                
                data_path_front{n} = [exp_path{i} '\PG2camera' trialID '_jaw_trajectory.mat'];
                n = n+1;
            end
        elseif body_check == 1 || body_check == 2
            if numel(DirList(j).name) <= 13
                continue;
            end
            if strcmp(DirList(j).name(1:9), 'PG1camera') && strcmp(DirList(j).name(end-13:end), 'trajectory.mat')
                if body_check == 1 && ~isempty(strfind(DirList(j).name, 'rightpaw'))
                    continue;
                end
                if body_check == 2 && isempty(strfind(DirList(j).name, 'rightpaw'))
                    continue;
                end   
                temp = load(DirList(j).name);
                trajectory = temp.trajectory;
                if body_check == 1
                    trialID = DirList(j).name(10:end-15);
                else
                    trialID = DirList(j).name(10:end-24);
                end
                if any(excID == str2double(trialID))  && exc_path == i
                    if exist('forelimb_exclude', 'dir') ~= 7
                        mkdir forelimb_exclude;
                    end
                    movefile(DirList(j).name, 'forelimb_exclude');
                    continue;
                end
                csvdata = load(['PG1gpio' trialID '.csv']);
                if any(csvdata(:, 1))
                    frame_onsite = find(csvdata(:, 1), 1);
                    try
                    trajectory_pre = trajectory(frame_onsite-pre*FrameRate:frame_onsite-1, 1:2);
                    catch
                    end
                    speed = sqrt(diff(trajectory_pre(:, 1)).^2+diff(trajectory_pre(:, 2)).^2);
                    acceleration = abs(diff(speed));
                    speeds(n) = sum(speed);
                    accelerations(n) = sum(acceleration);
                else
                    continue;
                end
                
                if body_check == 1
                    temp = load(['PG2camera' trialID '_trajectory.mat']);
                else
                    temp = load(['PG2camera' trialID '_rightpaw_trajectory.mat']);
                end
                trajectory = temp.trajectory;
                csvdata = load(['PG2gpio' trialID '.csv']);
                if any(csvdata(:, 1))
                    frame_onsite = find(csvdata(:, 1), 1);
                    try
                    trajectory_pre = trajectory(frame_onsite-pre*FrameRate:frame_onsite-1, 1:2);
                    catch
                    end
                    speed = sqrt(diff(trajectory_pre(:, 1)).^2+diff(trajectory_pre(:, 2)).^2);
                    acceleration = abs(diff(speed));
                    speedf(n) = sum(speed);
                    accelerationf(n) = sum(acceleration);
                else
                    continue;
                end
                
                if body_check == 1
                    data_path_front{n} = [exp_path{i} '\PG2camera' trialID '_trajectory.mat'];
                    data_path_side{n} = [exp_path{i} '\PG1camera' trialID '_trajectory.mat'];
                else
                    data_path_front{n} = [exp_path{i} '\PG2camera' trialID '_rightpaw_trajectory.mat'];
                    data_path_side{n} = [exp_path{i} '\PG1camera' trialID '_rightpaw_trajectory.mat'];
                end
                n = n+1;
            end
        end
        workbar(j/totalfile, ['file' num2str(j) '/' num2str(totalfile)], 'Progress');
        drawnow;
    end
end
cd(curpwd);
msgbox('Done !');

%%
curpwd = pwd;
if body_check == 3 || body_check == 2
    thr_speed = fit_gaussian2dist(speedf, nbins, 'Speed Front (pixel)');
    thr_acceleration = fit_gaussian2dist(accelerationf, nbins, 'Acceleration Front (pixel)');
    
    figure;
    scatter(speedf, accelerationf, 5, [0 0 0]);
    hold on;
    yl = ylim;
    xl = xlim;
    line([thr_speed; thr_speed], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
    line([xl' xl'], [thr_acceleration; thr_acceleration], 'Color', [1 0 0], 'LineStyle', '--');
    xlabel('Speed Front (pixel)');
    ylabel('Acceleration Front (pixel)');
    
%     IDexc = find(speedf > thr_speed(2) & accelerationf > thr_acceleration(2));
    IDexc = find(speedf > thr_speed(2) | accelerationf > thr_acceleration(2));
%     workbar(0, 'Moving Files...', 'Progress');
    parfor i = 1:numel(IDexc)
        [filepath, ~, ~] = fileparts(data_path_front{IDexc(i)});
        if ~strcmp(filepath, pwd)
            cd(filepath);
        end
        if body_check == 3
            if exist('jaw_exclude', 'dir') ~= 7
                mkdir jaw_exclude;
            end
            movefile(data_path_front{IDexc(i)}, 'jaw_exclude');
        else
            if exist('forelimb_exclude', 'dir') ~= 7
                mkdir forelimb_exclude;
            end
            movefile(data_path_front{IDexc(i)}, 'forelimb_exclude');
            movefile(data_path_side{IDexc(i)}, 'forelimb_exclude');
        end
%         workbar(i/numel(IDexc), ['file' num2str(i) '/' num2str(numel(IDexc))], 'Progress');
%         drawnow;
    end
elseif body_check == 1
    thr_speedf = fit_gaussian2dist(speedf, nbins, 'Speed Front (pixel)');
    thr_accelerationf = fit_gaussian2dist(accelerationf, nbins, 'Acceleration Front (pixel)');
    thr_speeds = fit_gaussian2dist(speeds, nbins, 'Speed Side (pixel)');
    thr_accelerations = fit_gaussian2dist(accelerations, nbins, 'Acceleration Side (pixel)');
    
    figure;
    subplot(2, 2, 1);
    scatter(speeds, speedf, 5, [0 0 0]);
    hold on;
    yl = ylim;
    xl = xlim;
    line([thr_speeds; thr_speeds], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
    line([xl' xl'], [thr_speedf; thr_speedf], 'Color', [1 0 0], 'LineStyle', '--');
    xlabel('Speed Side (pixel)');
    ylabel('Speed Front (pixel)');
    
    subplot(2, 2, 2);
    scatter(accelerations, accelerationf, 5, [0 0 0]);
    hold on;
    yl = ylim;
    xl = xlim;
    line([thr_accelerations; thr_accelerations], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
    line([xl' xl'], [thr_accelerationf; thr_accelerationf], 'Color', [1 0 0], 'LineStyle', '--');
    xlabel('Acceleration Side (pixel)');
    ylabel('Acceleration Front (pixel)');
    
    subplot(2, 2, 3);
    scatter(speeds, accelerations, 5, [0 0 0]);
    hold on;
    yl = ylim;
    xl = xlim;
    line([thr_speeds; thr_speeds], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
    line([xl' xl'], [thr_accelerations; thr_accelerations], 'Color', [1 0 0], 'LineStyle', '--');
    xlabel('Speed Side (pixel)');
    ylabel('Acceleration Side (pixel)');
    
    subplot(2, 2, 4);
    scatter(speedf, accelerationf, 5, [0 0 0]);
    hold on;
    yl = ylim;
    xl = xlim;
    line([thr_speedf; thr_speedf], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
    line([xl' xl'], [thr_accelerationf; thr_accelerationf], 'Color', [1 0 0], 'LineStyle', '--');
    xlabel('Speed Front (pixel)');
    ylabel('Acceleration Front (pixel)');
    
%     IDexc = find(speedf > thr_speedf(2) & accelerationf > thr_accelerationf(2) & speeds > thr_speeds(2) & accelerations > thr_accelerations(2));
%     IDexc = find((speedf > thr_speedf(2) & accelerationf > thr_accelerationf(2)) | (speeds > thr_speeds(2) & accelerations > thr_accelerations(2)));
    IDexc = find(speedf > thr_speedf(2) | accelerationf > thr_accelerationf(2) | speeds > thr_speeds(2) | accelerations > thr_accelerations(2));
%     workbar(0, 'Moving Files...', 'Progress');
    parfor i = 1:numel(IDexc)
        [filepath, ~, ~] = fileparts(data_path_front{IDexc(i)});
        if ~strcmp(filepath, pwd)
            cd(filepath);
        end
        if exist('forelimb_exclude', 'dir') ~= 7
            mkdir forelimb_exclude;
        end
        movefile(data_path_front{IDexc(i)}, 'forelimb_exclude');
        movefile(data_path_side{IDexc(i)}, 'forelimb_exclude');
%         workbar(i/numel(IDexc), ['file' num2str(i) '/' num2str(numel(IDexc))], 'Progress');
%         drawnow;
    end
end
cd(curpwd);
msgbox('Done !');