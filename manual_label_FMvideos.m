%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
trials = 1:20;
ntrials = numel(trials);

thr = 0.5;

curpwd = pwd;
try
    cd(foldername);
end
foldername = uigetdir('', 'Select a folder to process');
if foldername == 0
    return;
else
    cd(foldername);
end

nose_start1 = nan(ntrials, 3);
nose_stop1 = nan(ntrials, 3);
leftpaw_start1 = nan(ntrials, 3);
leftpaw_stop1 = nan(ntrials, 3);

nose_start2 = nan(ntrials, 3);
nose_stop2 = nan(ntrials, 3);
leftpaw_start2 = nan(ntrials, 3);
leftpaw_stop2 = nan(ntrials, 3);

hf = figure;
for i = 1:ntrials
    triggerstate = load(['PG1gpio' num2str(trials(i)) '.csv']);
    lightframeid = find(triggerstate == 1);
    vreader = VideoReader(['PG1camera' num2str(trials(i)) '.avi']);
    framestart = read(vreader, lightframeid(1)-1);
    framestop = read(vreader, lightframeid(end));
    if thr ~= 1
        vidFrame = rgb2hsv(framestart);
        intensity = vidFrame(:, :, 3);
        intensity = mat2gray(intensity, [0 thr]);
        vidFrame(:, :, 3) = intensity;
        framestart = hsv2rgb(vidFrame);
        
        vidFrame = rgb2hsv(framestop);
        intensity = vidFrame(:, :, 3);
        intensity = mat2gray(intensity, [0 thr]);
        vidFrame(:, :, 3) = intensity;
        framestop = hsv2rgb(vidFrame);
    end
    if i == 1
        subplot(1, 2, 1);
        him1 = imshow(framestart);
        title('Start');
        subplot(1, 2, 2);
        him2 = imshow(framestop);
        title('Stop');
    else
        him1.CData = framestart;
        him2.CData = framestop;
    end
    
    hf.Name = ['Trial' num2str(trials(i)) ' Nose Start'];
    [x, y, button] = ginput(1);
    nose_start1(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' LeftPaw Start'];
    [x, y, button] = ginput(1);
    leftpaw_start1(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' Nose Stop'];
    [x, y, button] = ginput(1);
    nose_stop1(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' LeftPaw Stop'];
    [x, y, button] = ginput(1);
    leftpaw_stop1(i, :) = [x y button];
    
    triggerstate = load(['PG2gpio' num2str(trials(i)) '.csv']);
    lightframeid = find(triggerstate == 1);
    vreader = VideoReader(['PG2camera' num2str(trials(i)) '.avi']);
    framestart = read(vreader, lightframeid(1)-1);
    framestop = read(vreader, lightframeid(end));
    him1.CData = framestart;
    him2.CData = framestop;
    
    hf.Name = ['Trial' num2str(trials(i)) ' Nose Start'];
    [x, y, button] = ginput(1);
    nose_start2(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' LeftPaw Start'];
    [x, y, button] = ginput(1);
    leftpaw_start2(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' Nose Stop'];
    [x, y, button] = ginput(1);
    nose_stop2(i, :) = [x y button];
    hf.Name = ['Trial' num2str(trials(i)) ' LeftPaw Stop'];
    [x, y, button] = ginput(1);
    leftpaw_stop2(i, :) = [x y button];
end
close(hf);
result.nose_start1 = nose_start1;
result.nose_stop1 = nose_stop1;
result.leftpaw_start1 = leftpaw_start1;
result.leftpaw_stop1 = leftpaw_stop1;
result.nose_start2 = nose_start2;
result.nose_stop2 = nose_stop2;
result.leftpaw_start2 = leftpaw_start2;
result.leftpaw_stop2 = leftpaw_stop2;
if exist('result.mat', 'file') ~= 2
    save('result.mat', 'result');
else
    user_response = modaldlg;
    switch lower(user_response)
        case 'no'
            % take no action
        case 'yes'
            save('result.mat', 'result');
    end
end
cd(curpwd);
msgbox('Done !');

%% combine results from different exp folders
n = 2;
result.nose_start1 = [];
result.nose_stop1 = [];
result.leftpaw_start1 = [];
result.leftpaw_stop1 = [];
result.nose_start2 = [];
result.nose_stop2 = [];
result.leftpaw_start2 = [];
result.leftpaw_stop2 = [];
for i = 1:n
    eval(['result.nose_start1 = [result.nose_start1; result' num2str(i) '.nose_start1];']);
    eval(['result.nose_stop1 = [result.nose_stop1; result' num2str(i) '.nose_stop1];']);
    eval(['result.leftpaw_start1 = [result.leftpaw_start1; result' num2str(i) '.leftpaw_start1];']);
    eval(['result.leftpaw_stop1 = [result.leftpaw_stop1; result' num2str(i) '.leftpaw_stop1];']);
    eval(['result.nose_start2 = [result.nose_start2; result' num2str(i) '.nose_start2];']);
    eval(['result.nose_stop2 = [result.nose_stop2; result' num2str(i) '.nose_stop2];']);
    eval(['result.leftpaw_start2 = [result.leftpaw_start2; result' num2str(i) '.leftpaw_start2];']);
    eval(['result.leftpaw_stop2 = [result.leftpaw_stop2; result' num2str(i) '.leftpaw_stop2];']);
end