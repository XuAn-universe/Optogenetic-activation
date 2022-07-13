function varargout = Video_Editor(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Video_Editor_OpeningFcn, ...
                   'gui_OutputFcn',  @Video_Editor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function Video_Editor_OpeningFcn(hObject, eventdata, handles, varargin)
image(rand(720, 1280, 3), 'Parent', handles.Movie_axes);
axis(handles.Movie_axes, 'image');
axis(handles.Movie_axes, 'off');
handles.output = hObject;
guidata(hObject, handles);

function varargout = Video_Editor_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Play_button_Callback(hObject, eventdata, handles)
readerobj = handles.readerobj;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = readerobj.FrameRate;
    FrameRate = round(FrameRate);
end
interval = 1/FrameRate;
Frames = readerobj.NumberOfFrames;
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex+1;
rotate_times = handles.rotate_times;
set(handles.Stop_button, 'UserData', 0);
while get(handles.Stop_button, 'UserData') == 0 && index <= Frames
    time0 = tic;
    im = read(readerobj, index);
    if rotate_times ~= 0
        im = imrotate(im, -90*rotate_times);
    end
    if ~isempty(handles.csvdata)
        if handles.csvdata(index, 1) == 1 && isempty(handles.hrect)
            handles.hrect = rectangle(handles.Movie_axes, 'Position', [0 0 handles.sqsize handles.sqsize], 'FaceColor', [1 1 1], 'EdgeColor', 'none');
        end
        if handles.csvdata(index, 1) ~= 1 && ~isempty(handles.hrect)
            delete(handles.hrect);
            handles.hrect = [];
        end
    end
    set(handles.him, 'CData', im);
    set(handles.frame_slider, 'Value', index);
    drawnow;
    index = index+1;
    while toc(time0) < interval
    end
%     pause(0.001);
    disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
end
currentindex = index-1;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);
disp(['Current Frame ' num2str(currentindex)]);

function LastFrame_button_Callback(hObject, eventdata, handles)
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-1;
set(handles.frame_slider, 'Value', index);
frame_slider_Callback(hObject, eventdata, handles)

function NextFrame_button_Callback(hObject, eventdata, handles)
readerobj = handles.readerobj;
Frames = readerobj.NumberOfFrames;
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex+1;
set(handles.frame_slider, 'Value', index);
frame_slider_Callback(hObject, eventdata, handles)

function FrameRate_text_Callback(hObject, eventdata, handles)

function FrameRate_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rotate_button_Callback(hObject, eventdata, handles)
handles.rotate_times = handles.rotate_times+1;
if handles.rotate_times == 4
    handles.rotate_times = 0;
end
im = imrotate(handles.currentim, -90);
handles.currentim = im;
handles.him = image(im, 'Parent', handles.Movie_axes);
axis(handles.Movie_axes, 'image');
handles.Movie_axes.Visible = 'off';
if ~isempty(handles.csvdata)
    if handles.csvdata(handles.currentindex, 1) == 1 && isempty(handles.hrect)
        handles.hrect = rectangle(handles.Movie_axes, 'Position', [0 0 handles.sqsize handles.sqsize], 'FaceColor', [1 1 1], 'EdgeColor', 'none');
    end
    if handles.csvdata(handles.currentindex, 1) ~= 1 && ~isemptye(handles.hrect)
        delete(handles.hrect);
        handles.hrect = [];
    end
end
guidata(hObject, handles);

function SetEnd_button_Callback(hObject, eventdata, handles)
currentindex = handles.currentindex;
set(handles.SetEnd_button, 'UserData', currentindex);

function SetStart_button_Callback(hObject, eventdata, handles)
currentindex = handles.currentindex;
set(handles.SetStart_button, 'UserData', currentindex);

function Stop_button_Callback(hObject, eventdata, handles)
set(handles.Stop_button, 'UserData', 1);

function SelectVideo_button_Callback(hObject, eventdata, handles)
handles.csvdata = [];
try
    delete(handles.hrect);
end
handles.hrect = [];
curpwd = pwd;
try
    cd(handles.pathname);
end
[filename, pathname] = uigetfile({'*.*'; '*.Mov'; '*.avi'; '*.mp4'; '*.mkv'; '*.rmvb'}, 'Select a video');
cd(curpwd);
if filename ~= 0
    set(handles.VideoName_text, 'String', [pathname filename]);
    [~, fname, ~] = fileparts([pathname filename]);
    try
        handles.csvdata = load([pathname fname(1:3) 'gpio' fname(10:end) '.csv']);
    end
    readerobj = VideoReader([pathname filename]);
    im = read(readerobj, 1);
    handles.him = image(im, 'Parent', handles.Movie_axes);
    axis(handles.Movie_axes, 'image');
    handles.Movie_axes.Visible = 'off';
    handles.sqsize = 90;
    if ~isempty(handles.csvdata)
        if handles.csvdata(1, 1) == 1
            handles.hrect = rectangle(handles.Movie_axes, 'Position', [0 0 handles.sqsize handles.sqsize], 'FaceColor', [1 1 1], 'EdgeColor', 'none');
        end
    end
    handles.readerobj = readerobj;
    handles.rotate_times = 0;
    handles.currentim = im;
    handles.currentindex = 1;
    handles.pathname = pathname;
    set(handles.Stop_button, 'UserData', 0);
    set(handles.trial_text, 'String', 'Trial 1');
    set(handles.trial_text, 'UserData', []);
    set(handles.success_radiobutton, 'Value', 0);
    set(handles.fail_radiobutton, 'Value', 0);
    set(handles.reachnumber_text, 'String', '0');
    set(handles.startframe_text, 'String', []);
    set(handles.endframe_text, 'String', []);
    set(handles.frame_slider, 'Max', readerobj.NumberOfFrames);
    set(handles.frame_slider, 'Min', 1);
    set(handles.frame_slider, 'Value', 1);
    set(handles.frame_slider, 'SliderStep', [1/readerobj.NumberOfFrames, 10/readerobj.NumberOfFrames]);
    set(handles.level1_radiobutton, 'Value', 0);
    set(handles.level2_radiobutton, 'Value', 0);
    set(handles.level3_radiobutton, 'Value', 0);
    set(handles.missreachnumber_text, 'String', '0');
    guidata(hObject, handles);
    disp('Ready !');
end

function VideoName_text_Callback(hObject, eventdata, handles)

function VideoName_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LastFrame_button_ButtonDownFcn(hObject, eventdata, handles)

function Playback_button_Callback(hObject, eventdata, handles)
readerobj = handles.readerobj;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = readerobj.FrameRate;
    FrameRate = round(FrameRate);
end
interval = 1/FrameRate;
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-1;
rotate_times = handles.rotate_times;
set(handles.Stop_button, 'UserData', 0);
while get(handles.Stop_button, 'UserData') == 0 && index >= 1
    time0 = tic;
    im = read(readerobj, index);
    if rotate_times ~= 0
        im = imrotate(im, -90*rotate_times);
    end
    if ~isempty(handles.csvdata)
        if handles.csvdata(index, 1) == 1 && isempty(handles.hrect)
            handles.hrect = rectangle(handles.Movie_axes, 'Position', [0 0 handles.sqsize handles.sqsize], 'FaceColor', [1 1 1], 'EdgeColor', 'none');
        end
        if handles.csvdata(index, 1) ~= 1 && ~isempty(handles.hrect)
            delete(handles.hrect);
            handles.hrect = [];
        end
    end
    set(handles.him, 'CData', im);
    set(handles.frame_slider, 'Value', index);
    drawnow;
    index = index-1;
    while toc(time0) < interval
    end
    disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
end
currentindex = index+1;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);
disp(['Current Frame ' num2str(currentindex)]);

function FirstFrame_button_Callback(hObject, eventdata, handles)
set(handles.frame_slider, 'Value', 1);
frame_slider_Callback(hObject, eventdata, handles)

function SaveVideo_button_Callback(hObject, eventdata, handles)
readerobj = handles.readerobj;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = readerobj.FrameRate;
    FrameRate = round(FrameRate);
end
startindex = get(handles.SetStart_button, 'UserData');
endindex = get(handles.SetEnd_button, 'UserData');
rotate_times = handles.rotate_times;
if startindex == endindex
    errordlg('Start frame = End frame!', 'ERROR');
else
    [filename, pathname] = uiputfile('*.mp4', 'Save as');
    if ~filename
        return;
    end
    ims = read(readerobj, [min(startindex, endindex) max(startindex, endindex)]);
    imvideo = imrotate(ims, -90*rotate_times);
    for i = 1:endindex-startindex+1
        id = i+startindex-1;
        if ~isempty(handles.csvdata)
            if handles.csvdata(id, 1)
                imvideo(:, :, :, i) = insertShape(imvideo(:, :, :, i), 'FilledRectangle', [0 0 handles.sqsize handles.sqsize], 'Color', 'w', 'Opacity', 1);
            end
        end
    end
    clear ims;
    vidObj = VideoWriter([pathname filename], 'MPEG-4');
    set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
    open(vidObj);
    writeVideo(vidObj, imvideo);
    close(vidObj);
    clear imvideo
    disp('Done!');
end

function StepBack1_button_Callback(hObject, eventdata, handles)
readerobj = handles.readerobj;
FrameRate = readerobj.FrameRate;
FrameRate = round(FrameRate);
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-FrameRate;
if index <= 0
    index = 1;
end
set(handles.frame_slider, 'Value', index);
frame_slider_Callback(hObject, eventdata, handles)

function Timekeep_text_Callback(hObject, eventdata, handles)

function Timekeep_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NextTrial_Callback(hObject, eventdata, handles)
trialid = get(handles.trial_text, 'String');
trialid = str2double(trialid(7:end));
trial = get(handles.trial_text, 'UserData');

if get(handles.success_radiobutton, 'Value')
    result.outcome = 1;
else
    if get(handles.fail_radiobutton, 'Value')
        result.outcome = 0;
        if get(handles.MissbutEat_checkbox, 'Value')
            result.missbuteat = 1;
        else
            result.missbuteat = 0;
        end
        result.MissatReach = str2double(get(handles.missreachnumber_text, 'String'));
        if result.MissatReach == 0
            errordlg('The number in the box of ''Miss @ Reach'' cannot be 0', 'Error');
            return;
        end
    else
        errordlg('Select whether it was a success or fail!', 'Error');
        return;
    end
end
result.reaches = str2double(get(handles.reachnumber_text, 'String'));
if result.reaches == 0
    errordlg('# of reaches should not be 0', 'Error');
    return;
end
result.startframe = str2double(get(handles.startframe_text, 'String'));
if isnan(result.startframe)
    errordlg('Please decide the start frame of the trial', 'Error');
    return;
end
result.endframe = str2double(get(handles.endframe_text, 'String'));
if isnan(result.endframe)
    errordlg('Please decide the end frame of the trial', 'Error');
    return;
end
if get(handles.level1_radiobutton, 'Value')
    result.level = 1;
elseif get(handles.level2_radiobutton, 'Value')
    result.level = 2;
elseif get(handles.level3_radiobutton, 'Value')
    result.level = 3;
else
    errordlg('Choose the level of the training', 'Error');
    return;
end
trial{trialid} = result;

set(handles.trial_text, 'UserData', trial);
set(handles.trial_text, 'String', ['Trial ' num2str(trialid+1)]);
if numel(trial) >= trialid+1
    result = trial{trialid+1};
    if result.outcome == 1
        set(handles.success_radiobutton, 'Value', 1);
        set(handles.fail_radiobutton, 'Value', 0);
        set(handles.MissbutEat_checkbox, 'Value', 0);
        set(handles.missreachnumber_text, 'String', '0');
    else
        set(handles.success_radiobutton, 'Value', 0);
        set(handles.fail_radiobutton, 'Value', 1);
        if result.missbuteat == 1
            set(handles.MissbutEat_checkbox, 'Value', 1);
        else
            set(handles.MissbutEat_checkbox, 'Value', 0);
        end
        set(handles.missreachnumber_text, 'String', num2str(result.MissatReach));
    end
    set(handles.reachnumber_text, 'String', num2str(result.reaches));
    set(handles.startframe_text, 'String', num2str(result.startframe));
    set(handles.endframe_text, 'String', num2str(result.endframe));
    switch result.level
        case 1
            set(handles.level1_radiobutton, 'Value', 1);
            set(handles.level2_radiobutton, 'Value', 0);
            set(handles.level3_radiobutton, 'Value', 0);
        case 2
            set(handles.level1_radiobutton, 'Value', 0);
            set(handles.level2_radiobutton, 'Value', 1);
            set(handles.level3_radiobutton, 'Value', 0);
        case 3
            set(handles.level1_radiobutton, 'Value', 0);
            set(handles.level2_radiobutton, 'Value', 0);
            set(handles.level3_radiobutton, 'Value', 1);
    end
    
    set(handles.frame_slider, 'Value', result.startframe);
    frame_slider_Callback(hObject, eventdata, handles)
else
    set(handles.success_radiobutton, 'Value', 0);
    set(handles.fail_radiobutton, 'Value', 0);
    set(handles.reachnumber_text, 'String', '0');
    set(handles.startframe_text, 'String', []);
    set(handles.endframe_text, 'String', []);
    set(handles.MissbutEat_checkbox, 'Value', 0);
    set(handles.missreachnumber_text, 'String', '0');
end

function LastTrial_Callback(hObject, eventdata, handles)
trialid = get(handles.trial_text, 'String');
trialid = str2double(trialid(7:end));
trial = get(handles.trial_text, 'UserData');
if trialid ~= 1
    result = trial{trialid-1};
    if result.outcome == 1
        set(handles.success_radiobutton, 'Value', 1);
        set(handles.fail_radiobutton, 'Value', 0);
        set(handles.MissbutEat_checkbox, 'Value', 0);
        set(handles.missreachnumber_text, 'String', '0');
    else
        set(handles.success_radiobutton, 'Value', 0);
        set(handles.fail_radiobutton, 'Value', 1);
        if result.missbuteat == 1
            set(handles.MissbutEat_checkbox, 'Value', 1);
        else
            set(handles.MissbutEat_checkbox, 'Value', 0);
        end
        set(handles.missreachnumber_text, 'String', num2str(result.MissatReach));
    end
    set(handles.reachnumber_text, 'String', num2str(result.reaches));
    set(handles.startframe_text, 'String', num2str(result.startframe));
    set(handles.endframe_text, 'String', num2str(result.endframe));
    switch result.level
        case 1
            set(handles.level1_radiobutton, 'Value', 1);
            set(handles.level2_radiobutton, 'Value', 0);
            set(handles.level3_radiobutton, 'Value', 0);
        case 2
            set(handles.level1_radiobutton, 'Value', 0);
            set(handles.level2_radiobutton, 'Value', 1);
            set(handles.level3_radiobutton, 'Value', 0);
        case 3
            set(handles.level1_radiobutton, 'Value', 0);
            set(handles.level2_radiobutton, 'Value', 0);
            set(handles.level3_radiobutton, 'Value', 1);
    end
    set(handles.trial_text, 'String', ['Trial ' num2str(trialid-1)]);
    
    set(handles.frame_slider, 'Value', result.startframe);
    frame_slider_Callback(hObject, eventdata, handles)
end

function triallength_text_Callback(hObject, eventdata, handles)

function triallength_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function frame_slider_Callback(hObject, eventdata, handles)
index = get(handles.frame_slider, 'Value');
index = fix(index);
readerobj = handles.readerobj;
rotate_times = handles.rotate_times;
im = read(readerobj, index);
im = imrotate(im, -90*rotate_times);
if ~isempty(handles.csvdata)
    if handles.csvdata(index, 1) == 1 && isempty(handles.hrect)
        handles.hrect = rectangle(handles.Movie_axes, 'Position', [0 0 handles.sqsize handles.sqsize], 'FaceColor', [1 1 1], 'EdgeColor', 'none');
    end
    if handles.csvdata(index, 1) ~= 1 && ~isempty(handles.hrect)
        delete(handles.hrect);
        handles.hrect = [];
    end
end
set(handles.him, 'CData', im);
currentindex = index;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);
disp(['Current Frame ' num2str(currentindex)]);

function frame_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function SaveResult_Callback(hObject, eventdata, handles)
curpwd = pwd;
try
    cd(handles.pathname);
end
data = get(handles.trial_text, 'UserData');
if ~isempty(data)
    [file, path] = uiputfile('*.mat','Save');
    if file ~= 0
        save([path file], 'data');
    end
end
cd(curpwd);

function startframe_text_Callback(hObject, eventdata, handles)

function startframe_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StartFrame_Callback(hObject, eventdata, handles)
set(handles.startframe_text, 'String', num2str(handles.currentindex));

function endframe_text_Callback(hObject, eventdata, handles)

function endframe_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EndFrame_Callback(hObject, eventdata, handles)
set(handles.endframe_text, 'String', num2str(handles.currentindex));

function reachnumber_text_Callback(hObject, eventdata, handles)

function reachnumber_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reachup_Callback(hObject, eventdata, handles)
reaches = str2double(get(handles.reachnumber_text, 'String'));
reaches = reaches+1;
set(handles.reachnumber_text, 'String', num2str(reaches));

function reachdown_Callback(hObject, eventdata, handles)
reaches = str2double(get(handles.reachnumber_text, 'String'));
reaches = reaches-1;
set(handles.reachnumber_text, 'String', num2str(max(reaches, 0)));

function LoadResult_Callback(hObject, eventdata, handles)
curpwd = pwd;
try
    cd(handles.pathname);
end
[file, path] = uigetfile('*.mat','Load');
if file ~= 0
    load([path file]);
    set(handles.trial_text, 'UserData', data);
    set(handles.trial_text, 'String', ['Trial ' num2str(numel(data)+1)]);
    set(handles.success_radiobutton, 'Value', 0);
    set(handles.fail_radiobutton, 'Value', 0);
    set(handles.reachnumber_text, 'String', '0');
    set(handles.startframe_text, 'String', []);
    set(handles.endframe_text, 'String', []);
    set(handles.level1_radiobutton, 'Value', 0);
    set(handles.level2_radiobutton, 'Value', 0);
    set(handles.level3_radiobutton, 'Value', 0);
    set(handles.MissbutEat_checkbox, 'Value', 0);
    set(handles.missreachnumber_text, 'String', '0');
    
    set(handles.frame_slider, 'Value', data{end}.endframe);
    frame_slider_Callback(hObject, eventdata, handles)
end
cd(curpwd);

function MissbutEat_checkbox_Callback(hObject, eventdata, handles)

function missreachnumber_text_Callback(hObject, eventdata, handles)

function missreachnumber_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MReachup_Callback(hObject, eventdata, handles)
missatreach = str2double(get(handles.missreachnumber_text, 'String'));
missatreach = missatreach+1;
set(handles.missreachnumber_text, 'String', num2str(missatreach));

function MReachdown_Callback(hObject, eventdata, handles)
missatreach = str2double(get(handles.missreachnumber_text, 'String'));
missatreach = missatreach-1;
set(handles.missreachnumber_text, 'String', num2str(max(missatreach, 0)));
