function varargout = ShowCropImage(varargin)
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
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShowCropImage_OpeningFcn, ...
                   'gui_OutputFcn',  @ShowCropImage_OutputFcn, ...
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

function ShowCropImage_OpeningFcn(hObject, eventdata, handles, varargin)
handles.masterhandles = varargin{1};
imageid = get(handles.masterhandles.image_list, 'Value');
channel = get(handles.masterhandles.colorchannel_list, 'Value');
im = handles.masterhandles.im;
handles.RGBorder = handles.masterhandles.RGBorder;
color_channel = cell(1, size(im, 3));
for i = 1:numel(channel)
    display_range = handles.masterhandles.display_ranges(channel(i), :, imageid);
    eval(['set(handles.radiobutton' num2str(channel(i)) ', ''UserData'', display_range);']);
    color_channel{channel(i)} = ['Channel ' num2str(channel(i))];
end
set(handles.colorchannel_list, 'String', color_channel);
set(handles.colorchannel_list, 'Value', channel);
handles = Show_Image(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

handles.output = hObject;
guidata(hObject, handles);

function handles = Show_Image(hObject, eventdata, handles)
if numel(handles.masterhandles.image_axes.Children) == 1
    cla(handles.image_axes);
    return;
end
imageid = get(handles.masterhandles.image_list, 'Value');
channel = get(handles.colorchannel_list, 'Value');
channel_names = get(handles.colorchannel_list, 'String');
for i = 1:numel(channel)
    if isempty(channel_names{channel(i)})
        return;
    end
end
im = handles.masterhandles.im;
for i = 1:numel(channel)
    temp = im(:, :, channel(i));
    temp = double(temp)/(2^16-1);
    if handles.masterhandles.rotation_angles(imageid) ~= 0
        temp = imrotate(temp, handles.masterhandles.rotation_angles(imageid), 'bicubic');
    end
    if handles.masterhandles.mirror_image(imageid) ~= 0
        temp = fliplr(temp);
    end
    eval(['display_range = get(handles.radiobutton' num2str(channel(i)) ', ''UserData'');']);
    temp = imadjust(temp, display_range, [0 1]);
    if i == 1
        im_display = zeros([handles.masterhandles.position(4:-1:3)+1, 3]);
    end
    im_display(:, :, handles.RGBorder(channel(i))) = imcrop(temp, handles.masterhandles.position);
    eval(['set(handles.radiobutton' num2str(channel(i)) ', ''Visible'', ''on'');']);
end
handles.im_display = im_display;
guidata(hObject, handles);
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
end
off_channel = setdiff([1 2 3], channel);
if ~isempty(off_channel)
    for i = 1:numel(off_channel)
        eval(['set(handles.radiobutton' num2str(off_channel(i)) ', ''Visible'', ''off'');']);
    end
end
clear im temp im_display;

function intensity_histogram(hObject, eventdata, handles)
cla(handles.imhist_axes);
im = handles.masterhandles.im;
imageid = get(handles.masterhandles.image_list, 'Value');
channel = get(handles.masterhandles.colorchannel_list, 'Value');
hcolors(1:3) = {'r' 'g' 'b'};
if get(handles.radiobutton1, 'Value')
    if ~isempty(intersect(1, channel))
        im = im(:, :, 1);
        display_range = get(handles.radiobutton1, 'UserData');
        hcolor = hcolors{handles.RGBorder(1)};
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton2, 'Value')
    if ~isempty(intersect(2, channel))
        im = im(:, :, 2);
        display_range = get(handles.radiobutton2, 'UserData');
        hcolor = hcolors{handles.RGBorder(2)};
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton3, 'Value')
    if ~isempty(intersect(3, channel))
        im = im(:, :, 3);
        display_range = get(handles.radiobutton3, 'UserData');
        hcolor = hcolors{handles.RGBorder(3)};
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
im = double(im)/(2^16-1);
if handles.masterhandles.rotation_angles(imageid) ~= 0
    im = imrotate(im, handles.masterhandles.rotation_angles(imageid), 'bicubic');
end
if handles.masterhandles.mirror_image(imageid) ~= 0
    im = fliplr(im);
end
im = imcrop(im, handles.masterhandles.position);
[count, x] = imhist(im);
x = x(min(find(count ~= 0)):max(find(count ~= 0)));
count = count(min(find(count ~= 0)):max(find(count ~= 0)));
count = count.*x;
stem(handles.imhist_axes, x, count, hcolor, 'Marker', 'none');
ylim(handles.imhist_axes, [0 max(count)]);
xlim(handles.imhist_axes, [x(1) x(end)]);
hold(handles.imhist_axes, 'on');
hlow = plot(handles.imhist_axes, [max(display_range(1), x(1)) max(display_range(1), x(1))], [0 max(count)], [':' hcolor]);
hhigh = plot(handles.imhist_axes, [min(display_range(2), x(end)) min(display_range(2), x(end))], [0 max(count)], [':' hcolor]);
set(handles.imhist_axes, 'UserData', [hlow hhigh]);
set(handles.MinValue_slider, 'Min', x(1));
set(handles.MinValue_slider, 'Max', x(end));
set(handles.MinValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MinValue_slider, 'Value', max(display_range(1), x(1)));
set(handles.minvalue_text, 'String', num2str(max(display_range(1), x(1))));
set(handles.MaxValue_slider, 'Min', x(1));
set(handles.MaxValue_slider, 'Max', x(end));
set(handles.MaxValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MaxValue_slider, 'Value', min(display_range(2), x(end)));
set(handles.maxvalue_text, 'String', num2str(min(display_range(2), x(end))));
clear im;

function varargout = ShowCropImage_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function MaxValue_slider_Callback(hObject, eventdata, handles)
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
    display_range = get(handles.radiobutton1, 'UserData');
    display_range(2) = maxvalue;
    set(handles.radiobutton1, 'UserData', display_range);
end
if get(handles.radiobutton2, 'Value')
    display_range = get(handles.radiobutton2, 'UserData');
    display_range(2) = maxvalue;
    set(handles.radiobutton2, 'UserData', display_range);
end
if get(handles.radiobutton3, 'Value')
    display_range = get(handles.radiobutton3, 'UserData');
    display_range(2) = maxvalue;
    set(handles.radiobutton3, 'UserData', display_range);
end
Show_Image(hObject, eventdata, handles);

function MaxValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function MinValue_slider_Callback(hObject, eventdata, handles)
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
    display_range = get(handles.radiobutton1, 'UserData');
    display_range(1) = minvalue;
    set(handles.radiobutton1, 'UserData', display_range);
end
if get(handles.radiobutton2, 'Value')
    display_range = get(handles.radiobutton2, 'UserData');
    display_range(1) = minvalue;
    set(handles.radiobutton2, 'UserData', display_range);
end
if get(handles.radiobutton3, 'Value')
    display_range = get(handles.radiobutton3, 'UserData');
    display_range(1) = minvalue;
    set(handles.radiobutton3, 'UserData', display_range);
end
Show_Image(hObject, eventdata, handles);

function MinValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

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

function radiobutton1_Callback(hObject, eventdata, handles)
intensity_histogram(hObject, eventdata, handles);

function radiobutton2_Callback(hObject, eventdata, handles)
intensity_histogram(hObject, eventdata, handles);

function radiobutton3_Callback(hObject, eventdata, handles)
intensity_histogram(hObject, eventdata, handles);

function SaveImage_Callback(hObject, eventdata, handles)
curpwd = pwd;
try
    cd(handles.masterhandles.DataDir);
end
[file, path] = uiputfile('*.tiff','Save');
if file ~= 0
	imwrite(handles.im_display, [path file], 'TIFF', 'Compression', 'none');
end
cd(curpwd);

function colorchannel_list_Callback(hObject, eventdata, handles)
Show_Image(hObject, eventdata, handles);

function colorchannel_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
