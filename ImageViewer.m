function varargout = ImageViewer(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2017
% xan@cshl.edu
% Version: 1.0
%
% Update 1: Xu An, July 2019
% xan@cshl.edu
% Enable color rendering
%
% Update 2: Xu An, Sep 2019
% xan@cshl.edu
% Enable manual cell marking

% Update 3: Xu An, Oct 2021
% xu.an@duke.edu
% Make it compatible with STP data set
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @untitled_OpeningFcn, ...
    'gui_OutputFcn',  @untitled_OutputFcn, ...
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

function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

handles.radiobutton_selected = 1;
handles.RGBorder = [1 2 3];
set(handles.figure1, 'KeyPressFcn', @quit_marking);
handles.hpoint = [];
guidata(hObject, handles);

function varargout = untitled_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function ImageFolder_Callback(hObject, eventdata, handles)
if ~isfield(handles, 'DataDir')
    DataDir = uigetdir('C:\', 'Select an image folder to analyze');
else
    DataDir = uigetdir(handles.DataDir, 'Select an image folder to analyze');
end
if DataDir ~= 0
    handles.DataDir = DataDir;
    set(handles.ImageDir, 'String', DataDir);
    set(handles.image_list, 'String', []);
    Update_Callback(hObject, eventdata, handles);
end

function scalefactor_text_Callback(hObject, eventdata, handles)

function scalefactor_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ImageResize_Callback(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
DirList = dir;
n = 0;
for i = 3:size(DirList)
    if strncmp(DirList(i).name(end-3:end), '.jp2', 4)
        n = n+1;
        image_name{n} = DirList(i).name;
    end
end
scalefactor = str2double(get(handles.scalefactor_text, 'String'));
if isnan(scalefactor)
    scalefactor = 0.1;
end
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(image_name)
    im = imread(image_name{i});
    if scalefactor < 1 && scalefactor > 0
        im = imresize(im, scalefactor);
    end
%     save([image_name{i}(1:end-4) '.mat'], 'im');
    for j = 1:size(im, 3)
        if j == 1
            imwrite(im(:,:,j), [image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none');
        else
            imwrite(im(:,:,j), [image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none', 'WriteMode', 'append');
        end
    end
    workbar(i/numel(image_name), [num2str(i) '/' num2str(numel(image_name))], 'Progress'); 
end
cd (curpwd);

function Update_Callback(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
DirList = dir;
n = 0;
image_name = [];
for i = 3:size(DirList)
    if strncmp(DirList(i).name(end-4:end), '.tiff', 5)
        n = n+1;
        image_name{n} = DirList(i).name(1:end-5);
        STP = 0;
    end
    
    if strncmp(DirList(i).name(end-3:end), '.jp2', 4) && ~isempty(strfind(DirList(i).name, 'StitchedImage_Z'))
        n = n+1;
        image_name{n} = DirList(i).name(1:end-4);
        STP = 1;
    end
end
if n == 0
    guidata(hObject, handles);
    cd(curpwd);
    return;
end
number_sort = zeros(1, n);
for i = 1:n
    if ~STP
        try
            name = image_name{i}(11:end);
            slideid = str2double(name(1:find(name == ' ')-1));
            if any(name == '-')
                imageid = str2double(name(find(name == '-')+1:end));
            elseif any(name == '_')
                imageid = str2double(name(find(name == '_')+1:end));
            end
            number_sort(i) = slideid*100+imageid;
        catch
            number_sort(i) = NaN;
        end
    else
        number_sort(i) = str2double(image_name{i}(16:18));
    end
end
[~, index] = sort(number_sort);
image_name = image_name(index);
set(handles.image_list, 'String', image_name);
set(handles.image_list, 'Value', []);
set(handles.colorchannel_list, 'String', []);
set(handles.colorchannel_list, 'Value', []);
try
    load('Adjustment_Parameters.mat');
    if ~exist('area_ratios', 'var')
        area_ratios = ones(1, n)*0.2;
    end
    if ~exist('thrs', 'var')
        thrs = zeros(n, 3)/0;
    end
    if ~exist('cell_detection', 'var')
        cell_detection = cell(3, n);
    end
    save('Adjustment_Parameters.mat', 'rotation_angles', 'mirror_image', 'display_ranges', 'area_ratios', 'thrs', 'cell_detection');
catch
    rotation_angles = zeros(1, n);
    mirror_image = zeros(1, n);
    display_ranges = zeros(3, 2, n);
    display_ranges(:, 2, :) = 1;
    area_ratios = ones(1, n)*0.2;
    thrs = zeros(n, 3)/0;
    cell_detection = cell(3, n);
    save('Adjustment_Parameters.mat', 'rotation_angles', 'mirror_image', 'display_ranges', 'area_ratios', 'thrs', 'cell_detection');
end
handles.rotation_angles = rotation_angles;
handles.mirror_image = mirror_image;
handles.display_ranges = display_ranges;
handles.area_ratios = area_ratios;
handles.thrs = thrs;
handles.cell_detection = cell_detection;
handles.STP = STP;
guidata(hObject, handles);
cd(curpwd);

function image_list_Callback(hObject, eventdata, handles)
if numel(get(handles.image_list, 'Value')) ~= 1
    return;
end
image_name = get(handles.image_list, 'String');
image_name = image_name{get(handles.image_list, 'Value')};
curpwd = pwd;
cd(handles.DataDir);
if get(handles.highres_check, 'Value')
    im = imread([image_name '.jp2']);
    msgbox('Image has been loaded!', '')
else
    try
        iminfo = imfinfo([image_name '.tiff'], 'tiff');
        im = zeros(iminfo(1).Height, iminfo(1).Width, numel(iminfo), 'uint16');
        for i = 1:numel(iminfo)
            im(:,:,i) = imread([image_name '.tiff'], 'tiff', 'Index', i);
        end
    catch
        load([image_name '.mat']);
    end
end
handles.im = im;
guidata(hObject, handles);
channel = get(handles.colorchannel_list, 'Value');
if isempty(channel)
    set(handles.colorchannel_list, 'Value', 1);
end
color_channel = cell(1, size(im, 3));
for i = 1:size(im, 3)
    color_channel{i} = ['Channel ' num2str(i)];
end
if max(channel) > size(im, 3)
    set(handles.colorchannel_list, 'Value', []);
    set(handles.colorchannel_list, 'String', color_channel);
    cd(curpwd);
    return;
end
set(handles.colorchannel_list, 'String', color_channel);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);
rotation_angle = handles.rotation_angles(get(handles.image_list, 'Value'));
set(handles.rotation_text, 'String', num2str(rotation_angle));
set(handles.MirrorImage_check, 'Value', handles.mirror_image(get(handles.image_list, 'Value')));
area_ratio = handles.area_ratios(get(handles.image_list, 'Value'));
set(handles.arearatio_text, 'String', num2str(area_ratio));
cd(curpwd);

function intensity_histogram(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
cla(handles.imhist_axes);
im = handles.im;
hcolors(1:3) = {'r' 'g' 'b'};
if get(handles.radiobutton1, 'Value')
    if ~isempty(intersect(1, channel))
        im = im(:, :, 1);
        channel = 1;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton2, 'Value')
    if ~isempty(intersect(2, channel))
        im = im(:, :, 2);
        channel = 2;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton3, 'Value')
    if ~isempty(intersect(3, channel))
        im = im(:, :, 3);
        channel = 3;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
hcolor = hcolors{handles.RGBorder(channel)};
im = double(im)/(2^16-1);
if handles.rotation_angles(imageid) ~= 0
    im = imrotate(im, handles.rotation_angles(imageid), 'bicubic');
end
if handles.mirror_image(imageid) ~= 0
    im = fliplr(im);
end
[count, x] = imhist(im);
x = x(min(find(count ~= 0)):max(find(count ~= 0)));
count = count(min(find(count ~= 0)):max(find(count ~= 0)));
count = count.*x;
stem(handles.imhist_axes, x, count, hcolor, 'Marker', 'none');
ylim(handles.imhist_axes, [0 max(count)]);
xlim(handles.imhist_axes, [x(1) x(end)]);
hold(handles.imhist_axes, 'on');
hlow = plot(handles.imhist_axes, [max(handles.display_ranges(channel, 1, imageid), x(1)) max(handles.display_ranges(channel, 1, imageid), x(1))], [0 max(count)], [':' hcolor]);
hhigh = plot(handles.imhist_axes, [min(handles.display_ranges(channel, 2, imageid), x(end)) min(handles.display_ranges(channel, 2, imageid), x(end))], [0 max(count)], [':' hcolor]);
set(handles.imhist_axes, 'UserData', [hlow hhigh]);
title(handles.imhist_axes, ['Image ' num2str(imageid)]);
set(handles.MinValue_slider, 'Min', x(1));
set(handles.MinValue_slider, 'Max', x(end));
set(handles.MinValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MinValue_slider, 'Value', max(handles.display_ranges(channel, 1, imageid), x(1)));
set(handles.minvalue_text, 'String', num2str(max(handles.display_ranges(channel, 1, imageid), x(1))));
set(handles.MaxValue_slider, 'Min', x(1));
set(handles.MaxValue_slider, 'Max', x(end));
set(handles.MaxValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MaxValue_slider, 'Value', min(handles.display_ranges(channel, 2, imageid), x(end)));
set(handles.maxvalue_text, 'String', num2str(min(handles.display_ranges(channel, 2, imageid), x(end))));
clear im;
    
function image_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function colorchannel_list_Callback(hObject, eventdata, handles)
persistent hrect
try
    delete(hrect);
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
im = handles.im;
for i = 1:numel(channel)
    temp = im(:, :, channel(i));
    temp = double(temp)/(2^16-1);
    if handles.rotation_angles(imageid) ~= 0
        temp = imrotate(temp, handles.rotation_angles(imageid), 'bicubic');
    end
    if handles.mirror_image(imageid) ~= 0
        temp = fliplr(temp);
    end
    if get(handles.BBox_check, 'Value')
        mask = temp > handles.display_ranges(channel(i), 2, imageid);
    end
    temp = imadjust(temp, [handles.display_ranges(channel(i), 1, imageid) handles.display_ranges(channel(i), 2, imageid)], [0 1]);
    if i == 1
        im_display = zeros([size(temp), 3]);
    end
    im_display(:, :, handles.RGBorder(channel(i))) = temp;
    eval(['set(handles.radiobutton' num2str(channel(i)) ', ''Visible'', ''on'');']);
end
try
    handles.image_axes.Children.CData = im_display;
catch
    cla(handles.image_axes);
    image(im_display, 'Parent', handles.image_axes);
    axis(handles.image_axes, 'image');
    axis(handles.image_axes, 'off');
end
if numel(channel) == 1
    if channel == 1
        set(handles.radiobutton1, 'Value', 1);
    end
    if channel == 2
        set(handles.radiobutton2, 'Value', 1);
    end
    if channel == 3
        set(handles.radiobutton3, 'Value', 1);
    end
    intensity_histogram(hObject, eventdata, handles);
    
    if get(handles.BBox_check, 'Value')
        se = strel('disk', 3);
        mask = imerode(mask, se);
        mask = imdilate(mask, se);
        stats = regionprops(mask, 'Area', 'BoundingBox');
        if numel(stats) ~= 0
            area = zeros(1, numel(stats));
            boundingbox = zeros(numel(stats), 4);
            for j = 1:numel(stats)
                area(j) = stats(j).Area;
                boundingbox(j, :) = stats(j).BoundingBox;
                boundingbox(j, 3:4) = boundingbox(j, 1:2)+boundingbox(j, 3:4);
            end
            maxarea = max(area);
            boundingbox(area < maxarea*handles.area_ratios(imageid), :) = [];
            if size(boundingbox, 1) ~= 1
                boundingbox_all(1:2) = min(boundingbox(:, 1:2));
                boundingbox_all(3:4) = max(boundingbox(:, 3:4));
            else
                boundingbox_all = boundingbox;
            end
            boundingbox_all(3:4) = boundingbox_all(3:4)-boundingbox_all(1:2);
            hold(handles.image_axes, 'on');
            hrect = rectangle('Position', boundingbox_all, 'EdgeColor', [1 1 1], 'LineWidth', 1, 'Parent', handles.image_axes);
        end
    end
end
off_channel = setdiff([1 2 3], channel);
if ~isempty(off_channel)
    for i = 1:numel(off_channel)
        eval(['set(handles.radiobutton' num2str(off_channel(i)) ', ''Visible'', ''off'');']);
    end
end
clear im temp im_display;

function colorchannel_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxValue_slider_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
maxvalue = get(handles.MaxValue_slider, 'Value');
minvalue = get(handles.MinValue_slider, 'Value');
if maxvalue <= minvalue
    maxvalue = minvalue+1/2^8;
    set(handles.MaxValue_slider, 'Value', maxvalue);
end
set(handles.maxvalue_text, 'String', num2str(maxvalue));
temp = get(handles.imhist_axes, 'UserData');
if isempty(temp)
    return;
end
hhigh = temp(2);
set(hhigh, 'XData', [maxvalue maxvalue]);
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
handles.display_ranges(channel, 2, imageid) = maxvalue;
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function MaxValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function MinValue_slider_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
maxvalue = get(handles.MaxValue_slider, 'Value');
minvalue = get(handles.MinValue_slider, 'Value');
if minvalue >= maxvalue
    minvalue = maxvalue-1/2^8;
    set(handles.MinValue_slider, 'Value', minvalue);
end
set(handles.minvalue_text, 'String', num2str(minvalue));
temp = get(handles.imhist_axes, 'UserData');
if isempty(temp)
    return;
end
hlow = temp(1);
set(hlow, 'XData', [minvalue minvalue]);
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
if ~get(handles.BBox_check, 'Value')
    handles.display_ranges(channel, 1, imageid) = minvalue;
    update_parameters(hObject, eventdata, handles);
end
colorchannel_list_Callback(hObject, eventdata, handles);

function MinValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function rotation_text_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
rotation_angle = str2double(get(handles.rotation_text, 'String'));
handles.rotation_angles(imageid) = rotation_angle;
set(handles.BBox_check, 'Value', 0);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
im_display = handles.image_axes.Children.CData;
[height, width, ~] = size(im_display);
im_display(round(height/2-height/500/2):round(height/2+height/500/2), :, :) = 1;
im_display(:, round(width/2-width/500/2):round(width/2+width/500/2), :) = 1;
handles.image_axes.Children.CData = im_display;
intensity_histogram(hObject, eventdata, handles);

function rotation_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotateImage_Callback(hObject, eventdata, handles)

function ImageDir_Callback(hObject, eventdata, handles)

function ImageDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radiobutton1_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 1;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function radiobutton2_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 2;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function radiobutton3_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 3;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function maxvalue_text_Callback(hObject, eventdata, handles)
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
if maxvalue > get(handles.MaxValue_slider, 'Max')
    maxvalue = get(handles.MaxValue_slider, 'Max');
    set(handles.maxvalue_text, 'String', num2str(maxvalue));
end
if maxvalue < get(handles.MaxValue_slider, 'Min')
    maxvalue = get(handles.MaxValue_slider, 'Min');
    set(handles.maxvalue_text, 'String', num2str(maxvalue));
end
set(handles.MaxValue_slider, 'Value', maxvalue);
MaxValue_slider_Callback(hObject, eventdata, handles);

function maxvalue_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minvalue_text_Callback(hObject, eventdata, handles)
minvalue = str2double(get(handles.minvalue_text, 'String'));
if minvalue < get(handles.MinValue_slider, 'Min')
    minvalue = get(handles.MinValue_slider, 'Min');
    set(handles.minvalue_text, 'String', num2str(minvalue));
end
if minvalue > get(handles.MinValue_slider, 'Max')
    minvalue = get(handles.MinValue_slider, 'Max');
    set(handles.minvalue_text, 'String', num2str(minvalue));
end
set(handles.MinValue_slider, 'Value', minvalue);
MinValue_slider_Callback(hObject, eventdata, handles);

function minvalue_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function highres_check_Callback(hObject, eventdata, handles)
image_list_Callback(hObject, eventdata, handles);

function SelectROI_Callback(hObject, eventdata, handles)
if numel(get(handles.image_list, 'Value')) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
himrect = imrect(handles.image_axes);
fcn = makeConstrainToRectFcn('imrect', get(handles.image_axes,'XLim'),get(handles.image_axes, 'YLim'));
setPositionConstraintFcn(himrect, fcn);
position = wait(himrect);
position = round(position);
delete(himrect);
hold(handles.image_axes, 'on');
rectangle('Position', position, 'Parent', handles.image_axes, 'EdgeColor', [1 1 1], 'LineWidth', 1);
handles.position = position;
ShowCropImage(handles);

function figure1_CreateFcn(hObject, eventdata, handles)

function SaveImage_Callback(hObject, eventdata, handles)
fcolor = get(handles.figure1, 'Color');
set(handles.figure1, 'Color', [0 0 0]);
im = getframe(handles.image_axes);
set(handles.figure1, 'Color', fcolor);
curpwd = pwd;
try
    cd(handles.masterhandles.DataDir);
end
[file, path] = uiputfile('*.tiff','Save');
if file ~= 0
    imwrite(im.cdata, [path file], 'TIFF', 'Compression', 'none');
end
if handles.STP
    if numel(get(handles.image_list, 'Value')) ~= 1
        return;
    end
    image_name = get(handles.image_list, 'String');
    image_name = image_name{get(handles.image_list, 'Value')};
    if ~isfile([handles.DataDir '\' image_name '.tiff'])
        imwrite(handles.im, [handles.DataDir '\' image_name '.tiff'], 'TIFF', 'Compression', 'none');
        msgbox('Done !');
    end
end
cd(curpwd);

function Magic_Callback(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
overflow_pixels = 10;
fthreshold = 10;
se = strel('disk', 3);
image_name = get(handles.image_list, 'String');
imageid = get(handles.image_list, 'Value');
channel = get(handles.colorchannel_list, 'Value');
if isempty(channel)
    errordlg('Please pick up the channels!')
    cd(curpwd);
    return;
end
try
    load('Adjustment_Parameters.mat');
catch
    errordlg('Adjustment_Parameters.mat is missing!')
    cd(curpwd);
    return;
end
[file, path] = uiputfile('*.avi','Save');
if file == 0
    cd(curpwd);
    return;
end
boundingbox_all = zeros(numel(imageid), 4);
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(imageid)
    disp(image_name{imageid(i)});
    try
        iminfo = imfinfo([image_name{imageid(i)} '.tiff'], 'tiff');
        im = zeros(iminfo(1).Height, iminfo(1).Width, numel(iminfo), 'uint16');
        for j = 1:numel(iminfo)
            im(:,:,j) = imread([image_name{imageid(i)} '.tiff'], 'tiff', 'Index', j);
        end
    catch
        load([image_name{imageid(i)} '.mat']);
    end
    if size(im, 3) < max(channel)
        continue;
    end
    thr = thrs(imageid(i), channel);
    thr = thr*(2^16-1);
    if numel(channel) ~= 1
        if isempty(find(~isnan(thr), 1))
            temp = im(:, :, handles.radiobutton_selected);
            thr = NaN;
        else
            temp = im(:, :, channel(find(~isnan(thr), 1)));
            thr = thr(find(~isnan(thr), 1));
        end
    else
        temp = im(:, :, channel);
    end
    if rotation_angles(imageid(i)) ~= 0
        temp = imrotate(temp, rotation_angles(imageid(i)), 'bicubic');
    end
    if mirror_image(imageid(i)) == 1
        temp = fliplr(temp);
    end
    if isnan(thr)
        [count, x] = imhist(temp, 2^16);
        count(1:2) = [];
        x(1:2) = [];
        for j = 1:numel(count)
            if count(end-j+1) >= fthreshold
                break;
            end
        end
        xfit = x(1:end-j+1);
        countfit = count(1:end-j+1);
        xlow = xfit(min(find(countfit == max(countfit))));
        countfit = countfit.*xfit;
        %     options = fitoptions('gauss2', 'Lower', [0 min(xfit) 0 0 xlow 0], 'Upper', [max(countfit) max(xfit) inf max(countfit) max(xfit) inf]);
        %     [curve, gof, o] = fit(xfit, countfit, 'gauss2', options);
        [curve, gof, o] = fit(xfit, countfit, 'gauss2');
        coeff = coeffvalues(curve);
        b1= coeff(2);
        if b1 >= max(xfit)
            b1 = 0;
        end
        b2= coeff(5);
        if b2 >= max(xfit)
            b2 = 0;
        end
        countfit = smooth(countfit, 5);
        countrange = countfit(xfit >= xlow & xfit <= max(b1, b2));
        xrange = xfit(xfit >= xlow & xfit <= max(b1, b2));
        thr = xrange(max(find(countrange == min(countrange))));
    end
%     figure;
%     plot(xfit, countfit, '-k');
%     hold on;
%     plot(curve);
%     figure;
%     plot(xrange, countrange);
    mask = temp > thr;
%     figure;
%     imshow(mask);
    mask = imerode(mask, se);
    mask = imdilate(mask, se);
%     figure;
%     imshow(mask);
    stats = regionprops(mask, 'Area', 'BoundingBox');
    area = zeros(1, numel(stats));
    boundingbox = zeros(numel(stats), 4);
    for j = 1:numel(stats)
        area(j) = stats(j).Area;
        boundingbox(j, :) = stats(j).BoundingBox;
        boundingbox(j, 3:4) = boundingbox(j, 1:2)+boundingbox(j, 3:4);
    end
    maxarea = max(area);
    boundingbox(area < maxarea*area_ratios(imageid(i)), :) = [];
    if size(boundingbox, 1) ~= 1
        boundingbox_all(i, 1:2) = min(boundingbox(:, 1:2));
        boundingbox_all(i, 3:4) = max(boundingbox(:, 3:4));
    else
        boundingbox_all(i, :) = boundingbox;
    end
    boundingbox_all(i, 3:4) = boundingbox_all(i, 3:4)-boundingbox_all(i, 1:2);
%     hold on;
%     rectangle('Position', boundingbox_all(i, :), 'EdgeColor', [1 1 1], 'LineWidth', 1);
    workbar(i/2/numel(imageid), 'Computing Ongoing...', 'Progress'); 
end

FrameRate = 5;
vidObj = VideoWriter([path file], 'Motion JPEG AVI');
set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
open(vidObj);
width = max(boundingbox_all(:, 3))+2*overflow_pixels;
height = max(boundingbox_all(:, 4))+2*overflow_pixels;
for i = 1:numel(imageid)
    try
        iminfo = imfinfo([image_name{imageid(i)} '.tiff'], 'tiff');
        im = zeros(iminfo(1).Height, iminfo(1).Width, numel(iminfo), 'uint16');
        for j = 1:numel(iminfo)
            im(:,:,j) = imread([image_name{imageid(i)} '.tiff'], 'tiff', 'Index', j);
        end
    catch
        load([image_name{imageid(i)} '.mat']);
    end
    if size(im, 3) < max(channel)
        continue;
    end
    for j = 1:numel(channel)
        temp = im(:, :, channel(j));
        temp = double(temp)/(2^16-1);
        if rotation_angles(imageid(i)) ~= 0
            temp = imrotate(temp, rotation_angles(imageid(i)), 'bicubic');
        end
        if mirror_image(imageid(i)) ~= 0
            temp = fliplr(temp);
        end
        temp = imadjust(temp, [max(display_ranges(channel(j), 1, imageid(i)), min(temp(:))) min(display_ranges(channel(j), 2, imageid(i)), max(temp(:)))], [0 1]);
        if j == 1
            im_display = zeros([size(temp), 3]);
        end
        im_display(:, :, handles.RGBorder(channel(j))) = temp;
    end
    boundingbox = boundingbox_all(i, :);
    top = round(boundingbox(2));
    left = round(boundingbox(1));
    bottom = round(boundingbox(2)+boundingbox(4))-1;
    right = round(boundingbox(1)+boundingbox(3))-1;
    top = top-ceil((height-boundingbox(4))/2);
    bottom = bottom+floor((height-boundingbox(4))/2);
    left = left-ceil((width-boundingbox(3))/2);
    right = right+floor((width-boundingbox(3))/2);
    w = size(im_display, 2);
    h = size(im_display, 1);
    imsave = im_display(max(1, top):min(bottom, h), max(1, left):min(right, w), :);
    if size(imsave, 1) < height && size(imsave, 2) < width
        imsave(height, width, :) = 0;
        imsave = circshift(imsave, [max(0, -top+1) max(0, -left+1)]);
    else
        if size(imsave, 1) < height
            imsave(height, :, :) = 0;
            imsave = circshift(imsave, [max(0, -top+1) 0]);
        end
        if size(imsave, 2) < width
            imsave(:, width, :) = 0;
            imsave = circshift(imsave, [0 max(0, -left+1)]);
        end
    end 
%     figure;
%     imshow(imsave);
    writeVideo(vidObj, imsave);
    workbar((i+numel(imageid))/2/numel(imageid), 'Computing Ongoing...', 'Progress'); 
end
cd(curpwd);
close(vidObj);
clear imsave;
disp('Done!');

function MirrorImage_check_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
handles.mirror_image(imageid) = get(handles.MirrorImage_check, 'Value');
set(handles.BBox_check, 'Value', 0);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function arearatio_text_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
area_ratio = str2double(get(handles.arearatio_text, 'String'));
handles.area_ratios(imageid) = area_ratio;
set(handles.BBox_check, 'Value', 1);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function arearatio_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BBox_check_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    set(handles.BBox_check, 'Value', 0);
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    set(handles.BBox_check, 'Value', 0);
    errordlg('Please choose only 1 Channel!');
    return;
end
if get(handles.BBox_check, 'Value')
    maxvalue = get(handles.MaxValue_slider, 'Value');
    set(handles.BBox_check, 'UserData', maxvalue);
    set(handles.image_list, 'Enable', 'off');
    set(handles.colorchannel_list, 'Enable', 'off');
    set(handles.MinValue_slider, 'Enable', 'off');
    set(handles.minvalue_text, 'Enable', 'off');
    colorchannel_list_Callback(hObject, eventdata, handles);
else
    set(handles.image_list, 'Enable', 'on');
    set(handles.colorchannel_list, 'Enable', 'on');
    set(handles.MinValue_slider, 'Enable', 'on');
    set(handles.minvalue_text, 'Enable', 'on');
    maxvalue = get(handles.BBox_check, 'UserData');
    set(handles.MaxValue_slider, 'Value', maxvalue);
    MaxValue_slider_Callback(hObject, eventdata, handles);
end

function WriteThr_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 channel!');
    return;
end
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
handles.thrs(imageid, channel) = maxvalue;
update_parameters(hObject, eventdata, handles);
msgbox('Done!');

function update_parameters(hObject, eventdata, handles)
guidata(hObject, handles);
display_ranges = handles.display_ranges;
rotation_angles = handles.rotation_angles;
mirror_image = handles.mirror_image;
area_ratios = handles.area_ratios;
thrs = handles.thrs;
cell_detection = handles.cell_detection;
curpwd = pwd;
cd(handles.DataDir);
save('Adjustment_Parameters.mat', 'rotation_angles', 'display_ranges', 'mirror_image', 'area_ratios', 'thrs', 'cell_detection');
cd(curpwd);

function Range_Apply_All_button_Callback(hObject, eventdata, handles)
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
minvalue = str2double(get(handles.minvalue_text, 'String'));
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
handles.display_ranges(channel, 1, :) = minvalue;
handles.display_ranges(channel, 2, :) = maxvalue;
update_parameters(hObject, eventdata, handles);

function Area_Apply_All_button_Callback(hObject, eventdata, handles)
area_ratios = str2double(get(handles.arearatio_text, 'String'));
handles.area_ratios(:) = area_ratios;
update_parameters(hObject, eventdata, handles);

function popupmenu1_Callback(hObject, eventdata, handles)
handles.RGBorder(1) = get(handles.popupmenu1, 'Value');
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu2_Callback(hObject, eventdata, handles)
handles.RGBorder(2) = get(handles.popupmenu2, 'Value');
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu3_Callback(hObject, eventdata, handles)
handles.RGBorder(3) = get(handles.popupmenu3, 'Value');
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mark_Cells_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
handles.figure1.UserData = 1;
while handles.figure1.UserData
    h = drawpoint(handles.image_axes, 'Color', [1 1 1]);
    handles.hpoint = [handles.hpoint h];
end
guidata(hObject, handles);

function quit_marking(hObject, eventdata, handles)
if isequal(eventdata.Key, 'escape')
    hObject.UserData = 0;
end

function Save_Positions_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
temp = [];
for i = 1:numel(handles.hpoint)
    if isvalid(handles.hpoint(i))
        temp = [temp; handles.hpoint(i).Position];
    end
end
handles.cell_detection{channel, imageid} = temp;
update_parameters(hObject, eventdata, handles);
delete(handles.hpoint);
handles.hpoint = [];
guidata(hObject, handles);
msgbox([num2str(size(temp, 1)) ' Cell Position Saved']);

function Recover_Markers_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
if ~isempty(handles.cell_detection{channel, imageid})
    handles.hpoint = images.roi.Point.empty;
    temp = handles.cell_detection{channel, imageid};
    for i = 1:size(temp, 1)
        handles.hpoint(i) = drawpoint(handles.image_axes, 'Color', [1 1 1], 'Position', temp(i, :));
    end
end
guidata(hObject, handles);
