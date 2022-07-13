function varargout = Video_Editor_Fast(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Video_Editor_Fast_OpeningFcn, ...
                   'gui_OutputFcn',  @Video_Editor_Fast_OutputFcn, ...
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

function Video_Editor_Fast_OpeningFcn(hObject, eventdata, handles, varargin)
axis(handles.Movie_axes, 'off');
handles.output = hObject;
guidata(hObject, handles);

function varargout = Video_Editor_Fast_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Play_button_Callback(hObject, eventdata, handles)
global im_all
readerobj = handles.readerobj;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = readerobj.FrameRate;
    FrameRate = round(FrameRate);
end
interval = 1/FrameRate;
Frames = size(im_all, 3);
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex+1;
rotate_times = handles.rotate_times;
set(handles.Stop_button, 'UserData', 0);
stop_status = get(handles.Stop_button, 'UserData');
while stop_status == 0 && index <= Frames
    time0 = tic;
    stop_status = get(handles.Stop_button, 'UserData');
    im = im_all(:, :, index);
    im = imrotate(im, -90*rotate_times);
    set(handles.him, 'CData', im);
    drawnow;
    index = index+1;
    while toc(time0) < interval
    end
    disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
end
currentindex = index-1;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);

function LastFrame_button_Callback(hObject, eventdata, handles)
global im_all
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-1;
rotate_times = handles.rotate_times;
im = im_all(:, :, index);
im = imrotate(im, -90*rotate_times);
set(handles.him, 'CData', im);
currentindex = index;
disp(['Current Frame: ' num2str(index)]);
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);

function NextFrame_button_Callback(hObject, eventdata, handles)
global im_all
Frames = size(im_all, 3);
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex+1;
rotate_times = handles.rotate_times;
im = im_all(:, :, index);
im = imrotate(im, -90*rotate_times);
set(handles.him, 'CData', im);
currentindex = index;
disp(['Current Frame: ' num2str(index)]);
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);

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
him = imshow(im, [], 'Parent', handles.Movie_axes);
handles.him = him;
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
global im_all
[filename, pathname] = uigetfile({'*.*'; '*.Mov'; '*.avi'; '*.mp4'; '*.mkv'; '*.rmvb'}, 'Select a video');
if filename ~= 0
    set(handles.VideoName_text, 'String', [pathname filename]);
    readerobj = VideoReader([pathname filename]);
    Duration = readerobj.Duration;
    FrameRate = readerobj.FrameRate;
    Frames = round(Duration*FrameRate);
    Height = readerobj.Height;
    Width = readerobj.Width;
    bin = str2double(get(handles.bin_text, 'String'));
    if mod(Height, bin) || mod(Width, bin)
        errordlg('Check your Bin !', 'ERROR');
        return;
    end
    im_all = uint8(zeros(Height/bin, Width/bin, Frames));
    n = 0;
    while hasFrame(readerobj)
        n = n+1;
        im = readFrame(readerobj);
        im = rgb2gray(im);
        temp = zeros(size(im)/bin);
        for i = 1:bin
            for j = 1:bin
                temp = temp+double(im(i:bin:end, j:bin:end));
            end
        end
        temp = temp/bin^2;
        im_all(:, :, n) = uint8(temp);
    end
    him = imshow(im_all(:, :, 1), [], 'Parent', handles.Movie_axes);
    handles.readerobj = readerobj;
    handles.rotate_times = 0;
    handles.currentim = im;
    handles.currentindex = 1;
    handles.him = him;
    set(handles.Stop_button, 'UserData', 0);
    set(handles.SetStart_button, 'UserData', 1);
    set(handles.SetEnd_button, 'UserData', Frames);
    guidata(hObject, handles);
    disp('Done!');
end

function VideoName_text_Callback(hObject, eventdata, handles)

function VideoName_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LastFrame_button_ButtonDownFcn(hObject, eventdata, handles)

function Playback_button_Callback(hObject, eventdata, handles)
global im_all
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
stop_status = get(handles.Stop_button, 'UserData');
while stop_status == 0 && index >= 1
    time0 = tic;
    stop_status = get(handles.Stop_button, 'UserData');
    im = im_all(:, :, index);
    im = imrotate(im, -90*rotate_times);
    set(handles.him, 'CData', im);
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

function FirstFrame_button_Callback(hObject, eventdata, handles)
global im_all
rotate_times = handles.rotate_times;
im = im_all(:, :, 1);
im = imrotate(im, -90*rotate_times);
set(handles.him, 'CData', im);
currentindex = 1;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);

function SaveVideo_button_Callback(hObject, eventdata, handles)
global im_all
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
    [filename, pathname] = uiputfile('*.avi', 'Save as');
    if ~filename
        return;
    end
    vidObj = VideoWriter([pathname filename], 'Motion JPEG AVI');
    set(vidObj, 'FrameRate', FrameRate, 'Quality', 75);
    open(vidObj);
    if get(handles.bin_check, 'Value')
        ims = im_all(:, :, min(startindex, endindex):max(startindex, endindex));
        imvideo = imrotate(ims, -90*rotate_times);
        clear ims;
        imvideo = reshape(imvideo, size(imvideo, 1), size(imvideo, 2), 1, size(imvideo, 3));
        writeVideo(vidObj, imvideo);
    else
        readerobj = VideoReader(get(handles.VideoName_text, 'String'));
        ims = read(readerobj, [min(startindex, endindex) max(startindex, endindex)]);
        for i = 1:size(ims, 4)
            imvideo = imrotate(ims(:, :, :, i), -90*rotate_times);
            if size(imvideo, 1) > 1080
                bin = 2;
                temp = zeros(size(imvideo, 1)/bin, size(imvideo, 2)/bin, size(imvideo, 3));
                for m = 1:bin
                    for n = 1:bin
                        temp = temp+double(imvideo(m:bin:end, n:bin:end, :));
                    end
                end
                temp = temp/bin^2;
                imvideo = uint8(temp);
            end
%             if i >= 66-startindex+1 && i <=93-startindex+1
%                 imvideo(end-140:end, 1:140, 1:2) = 0;
%                 imvideo(end-140:end, 1:140, 3) = 255;
%             end
            writeVideo(vidObj, imvideo);
        end
        clear ims;
    end
    close(vidObj);
    clear imvideo;
    disp('Done!');
end

function StepBack1_button_Callback(hObject, eventdata, handles)
global im_all
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
rotate_times = handles.rotate_times;
im = im_all(:, :, index);
im = imrotate(im, -90*rotate_times);
set(handles.him, 'CData', im);
currentindex = index;
handles.currentindex = currentindex;
handles.currentim = im;
guidata(hObject, handles);

function bin_text_Callback(hObject, eventdata, handles)

function bin_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bin_check_Callback(hObject, eventdata, handles)
