function varargout = Object_Tracking_Copy(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Object_Tracking_Copy_OpeningFcn, ...
                   'gui_OutputFcn',  @Object_Tracking_Copy_OutputFcn, ...
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

function Object_Tracking_Copy_OpeningFcn(hObject, eventdata, handles, varargin)
handles.mouseline = 'Fezf2, 0.2x speed'; % for publication
handles.thr = 1; % brightness 0.6

handles.v = varargin{1};
handles.rect = [];
lightframeid = find(handles.v.data(:, 1) == 1);
Frames = handles.v.FrameRange(1)+numel(lightframeid)+handles.v.FrameRange(2);
skipframes = lightframeid(1)-handles.v.FrameRange(1)-1;
videoinfo = info(handles.v.FileReader);
FrameRate = videoinfo.VideoFrameRate;
Height = videoinfo.VideoSize(2);
Width = videoinfo.VideoSize(1);
handles.Frames = Frames;
handles.skipframes = skipframes;
handles.FrameRate = FrameRate;
im_all = single(zeros(Height, Width, 3, Frames)/0);
sim_all = single(zeros(Height, Width, 3, Frames)/0);
trajectory = zeros(Frames, 2)/0;
for i = 1:skipframes
    step(handles.v.FileReader);
end
im_all(:, :, :, 1) = step(handles.v.FileReader);
him = imshow(boost_IM_brightness(im_all(:, :, :, 1), handles.thr), [], 'Parent', handles.Movie_axes);
temp = handles.v.FileName;
textid = find(temp == '\');
temp(textid) = ':';
title(handles.Movie_axes, temp(textid(end-2)+1:end), 'interpreter', 'none');
handles.rotate_times = 0;
handles.currentindex = 1;
handles.currentim = im_all(:, :, :, 1);
handles.him = him;
handles.method = 1;
handles.im_all = im_all;
handles.sim_all = sim_all;
handles.trajectory = trajectory;
set(handles.Stop_button, 'UserData', 0);
cla(handles.plot_axes);
hold(handles.plot_axes, 'on');
xlabel(handles.plot_axes, 'dx (pixels)');
ylabel(handles.plot_axes, 'dy (pixels)');
title(handles.plot_axes, 'Trajectory');
handles.output = hObject;
guidata(hObject, handles);
disp(['Current Frame: ' num2str(1)]);
disp('Done!');

function varargout = Object_Tracking_Copy_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles;
zoominh = findall(handles.figure1,'tag','Exploration.ZoomIn');
set(zoominh, 'State', 'on');
xu_putdowntext('zoomin', zoominh);

function Play_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
skipframes = handles.skipframes;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = handles.FrameRate;
end
interval = 1/FrameRate;
Frames = handles.Frames;
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex;
rotate_times = handles.rotate_times;
im = im_all(:, :, :, index);
if rotate_times ~= 0
    im = imrotate(im, -90*rotate_times);
end
hblob = vision.BlobAnalysis('AreaOutputPort', false, 'CentroidOutputPort', true);
set(handles.Stop_button, 'UserData', 0);
stop_status = 0;
if get(handles.FT_radiobutton, 'Value')
    switch handles.method
        case 1
            points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 2
            points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 3
            points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 4
            points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 5
            points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 6
            points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
    end
    pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
    points = points.Location;
    initialize(pointTracker, points, im);
    oldPoints = points;
    
    x = handles.rect(1, 1);
    y = handles.rect(1, 2);
    w = handles.rect(1, 3);
    h = handles.rect(1, 4);
    bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
    
    while stop_status == 0 && index <= Frames-1
        time0 = tic;
        stop_status = get(handles.Stop_button, 'UserData');
        
        % get the next frame
        index = index+1;
        if isnan(im_all(1, 1, 1, index))
            im = step(handles.v.FileReader);
            im_all(:, :, :, index) = im;
        else
            im = im_all(:, :, :, index);
        end
        if rotate_times ~= 0
            im = imrotate(im, -90*rotate_times);
        end
        
        % Track the points. Note that some points may be lost.
        [points, isFound] = step(pointTracker, im);
        visiblePoints = points(isFound, :);
        oldInliers = oldPoints(isFound, :);
        
        if size(visiblePoints, 1) >= 2 % need at least 2 points
            
            % Estimate the geometric transformation between the old points
            % and the new points and eliminate outliers
            [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
            
            % Apply the transformation to the bounding box
            [bboxPolygon(1:2:end), bboxPolygon(2:2:end)] ...
                = transformPointsForward(xform, bboxPolygon(1:2:end), bboxPolygon(2:2:end));
            
            BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(im, 1), size(im, 2));
            centroid = step(hblob, BW);
            trajectory(index, :) = centroid;
            simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
            if handles.v.data(skipframes+index, 1)
                simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
                simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
            elseif skipframes+index < find(handles.v.data(:, 1), 1)
                simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
            else
                simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
            end
            sim_all(:, :, :, index) = simage;
            
            % Reset the points
            oldPoints = visiblePoints;
            setPoints(pointTracker, oldPoints);
        else
            release(pointTracker);
            currentindex = index-1;
            handles.im_all = im_all;
            handles.sim_all = sim_all;
            handles.trajectory = trajectory;
            refresh_handles(hObject, handles, trajectory, currentindex, im_all);
            disp('Points are not enough!');
            return;
        end
        
        % Display the annotated video frame using the video player object
        set(handles.him, 'CData', simage);
        drawnow;
        
        currentindex = index;
        handles.currentindex = currentindex;
        if get(handles.RT_check, 'Value')
            handles.trajectory = trajectory;
            plot_trajectory(handles);
        end
        
        while toc(time0) < interval
        end
        disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
    end
    release(pointTracker);
elseif get(handles.CC_radiobutton, 'Value')
    frame = rgb2hsv(im);
    frame = frame(:, :, 3);
    ROI = imcrop(frame, handles.rect);
    [rROI, cROI] = size(ROI);
    lastframe = frame;
    while stop_status == 0 && index <= Frames-1
        time0 = tic;
        stop_status = get(handles.Stop_button, 'UserData');
        
        % get the next frame
        index = index+1;
        if isnan(im_all(1, 1, 1, index))
            im = step(handles.v.FileReader);
            im_all(:, :, :, index) = im;
        else
            im = im_all(:, :, :, index);
        end
        if rotate_times ~= 0
            im = imrotate(im, -90*rotate_times);
        end
        frame = rgb2hsv(im);
        frame = frame(:, :, 3);
        crr = xcorr2(frame-mean2(frame), ROI-mean2(lastframe));
        if isfield(handles, 'Mask')
            if size(crr, 1) > size(handles.Mask, 1)
                diff = size(crr, 1)-size(handles.Mask, 1);
                handles.Mask(end+1:end+diff, :) = false;
            elseif size(crr, 1) < size(handles.Mask, 1)
                diff = size(handles.Mask, 1)-size(crr, 1);
                handles.Mask(1:diff,:) = [];
            end
            if size(crr, 2) > size(handles.Mask, 2)
                diff = size(crr, 2)-size(handles.Mask, 2);
                handles.Mask(:, end+1:end+diff) = false;
            elseif size(crr, 2) < size(handles.Mask, 2)
                diff = size(handles.Mask, 2)-size(crr, 2);
                handles.Mask(:, 1:diff) = [];
            end
            temp = crr.*handles.Mask;
        else
            temp = crr;
        end
        [~, maxindex] = max(temp(:));
        [row, column] = ind2sub(size(temp), maxindex);
        ROI = frame(row-rROI+1:row, column-cROI+1:column);
        lastframe = frame;
        BW = false(size(frame));
        BW(row-rROI+1:row, column-cROI+1:column) = 1;
        centroid = step(hblob, BW);
        trajectory(index, :) = centroid;
        simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
        if handles.v.data(skipframes+index, 1)
            simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        elseif skipframes+index < find(handles.v.data(:, 1), 1)
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        else
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        end
        sim_all(:, :, :, index) = simage;
        
        % Display the annotated video frame using the video player object
        set(handles.him, 'CData', simage);
        drawnow;
        
        currentindex = index;
        handles.currentindex = currentindex;
        if get(handles.RT_check, 'Value')
            handles.trajectory = trajectory;
            plot_trajectory(handles);
        end
        
        while toc(time0) < interval
        end
        disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
    end
end
if isDone(handles.v.FileReader)
    release(handles.v.FileReader);
end
handles.im_all = im_all;
handles.sim_all = sim_all;
handles.trajectory = trajectory;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
disp('Done!');
if ~get(handles.RT_check, 'Value')
    plot_trajectory(handles);
end

function LastFrame_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-1;
simage = sim_all(:, :, :, index);
set(handles.him, 'CData', simage);
currentindex = index;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
plot_trajectory(handles);

function refresh_handles(hObject, handles, trajectory, currentindex, im_all)
w = handles.rect(1, 3);
h = handles.rect(1, 4);
handles.rect = [trajectory(currentindex, 1)-w/2 trajectory(currentindex, 2)-h/2 w h];
handles.currentindex = currentindex;
handles.currentim = im_all(:, :, :, currentindex);
guidata(hObject, handles);
disp(['Current Frame: ' num2str(currentindex)]);

function NextFrame_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
skipframes = handles.skipframes;
Frames = handles.Frames;
currentindex = handles.currentindex;
if currentindex == Frames
    return;
end
index = currentindex;
rotate_times = handles.rotate_times;
im = im_all(:, :, :, index);
if rotate_times ~= 0
    im = imrotate(im, -90*rotate_times);
end
hblob = vision.BlobAnalysis('AreaOutputPort', false, 'CentroidOutputPort', true);
if get(handles.FT_radiobutton, 'Value')
    switch handles.method
        case 1
            points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 2
            points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 3
            points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 4
            points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 5
            points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
        case 6
            points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(handles.rect));
    end
    pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
    points = points.Location;
    initialize(pointTracker, points, im);
    oldPoints = points;
    
    x = handles.rect(1, 1);
    y = handles.rect(1, 2);
    w = handles.rect(1, 3);
    h = handles.rect(1, 4);
    bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
    
    % get the next frame
    index = index+1;
    if isnan(im_all(1, 1, 1, index))
        im = step(handles.v.FileReader);
        im_all(:, :, :, index) = im;
    else
        im = im_all(:, :, :, index);
    end
    if rotate_times ~= 0
        im = imrotate(im, -90*rotate_times);
    end
    
    % Track the points. Note that some points may be lost.
    [points, isFound] = step(pointTracker, im);
    visiblePoints = points(isFound, :);
    oldInliers = oldPoints(isFound, :);
    
    if size(visiblePoints, 1) >= 2 % need at least 2 points
        
        % Estimate the geometric transformation between the old points
        % and the new points and eliminate outliers
        [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
            oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
        
        % Apply the transformation to the bounding box
        [bboxPolygon(1:2:end), bboxPolygon(2:2:end)] ...
            = transformPointsForward(xform, bboxPolygon(1:2:end), bboxPolygon(2:2:end));
        
        BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(im, 1), size(im, 2));
        centroid = step(hblob, BW);
        trajectory(index, :) = centroid;
        simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
        if handles.v.data(skipframes+index, 1)
            simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        elseif skipframes+index < find(handles.v.data(:, 1), 1)
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        else
            simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
        end
        sim_all(:, :, :, index) = simage;
        
        % Display the annotated video frame using the video player object
        set(handles.him, 'CData', simage);
        drawnow;
    else
        disp('Points are not enough!');
        release(pointTracker);
        handles.im_all = im_all;
        guidata(hObject, handles);
        return;
    end
    release(pointTracker);
elseif get(handles.CC_radiobutton, 'Value')
    frame = rgb2hsv(im);
    frame = frame(:, :, 3);
    ROI = imcrop(frame, handles.rect);
    [rROI, cROI] = size(ROI);
    lastframe = frame;
    
    % get the next frame
    index = index+1;
    if isnan(im_all(1, 1, 1, index))
        im = step(handles.v.FileReader);
        im_all(:, :, :, index) = im;
    else
        im = im_all(:, :, :, index);
    end
    if rotate_times ~= 0
        im = imrotate(im, -90*rotate_times);
    end
    frame = rgb2hsv(im);
    frame = frame(:, :, 3);
    crr = xcorr2(frame-mean2(frame), ROI-mean2(lastframe));
    if isfield(handles, 'Mask')
        temp = crr.*handles.Mask;
    else
        temp = crr;
    end
    [~, maxindex] = max(temp(:));
    [row, column] = ind2sub(size(temp), maxindex);
    BW = false(size(frame));
    BW(row-rROI+1:row, column-cROI+1:column) = 1;
    centroid = step(hblob, BW);
    trajectory(index, :) = centroid;
    simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
    if handles.v.data(skipframes+index, 1)
        simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
        simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
    elseif skipframes+index < find(handles.v.data(:, 1), 1)
        simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
    else
        simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
    end
    sim_all(:, :, :, index) = simage;
    
    % Display the annotated video frame using the video player object
    set(handles.him, 'CData', simage);
    drawnow;
end
if isDone(handles.v.FileReader)
    release(handles.v.FileReader);
end
currentindex = index;
handles.im_all = im_all;
handles.sim_all = sim_all;
handles.trajectory = trajectory;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
disp('Done');
plot_trajectory(handles);

function FrameRate_text_Callback(hObject, eventdata, handles)

function FrameRate_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rotate_button_Callback(hObject, eventdata, handles)
if handles.currentindex == 1
    handles.rotate_times = handles.rotate_times+1;
    if handles.rotate_times == 4
        handles.rotate_times = 0;
    end
    simage = imrotate(boost_IM_brightness(handles.currentim, handles.thr), -90*handles.rotate_times);
    set(handles.him, 'CData', simage);
    drawnow;
    guidata(hObject, handles);
end

function SetEnd_button_Callback(hObject, eventdata, handles)
currentindex = handles.currentindex;
set(handles.SetEnd_button, 'UserData', currentindex);

function SetStart_button_Callback(hObject, eventdata, handles)
currentindex = handles.currentindex;
set(handles.SetStart_button, 'UserData', currentindex);

function Stop_button_Callback(hObject, eventdata, handles)
set(handles.Stop_button, 'UserData', 1);

function VideoName_text_Callback(hObject, eventdata, handles)

function VideoName_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LastFrame_button_ButtonDownFcn(hObject, eventdata, handles)

function Playback_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = handles.FrameRate;
end
interval = 1/FrameRate;
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-1;
set(handles.Stop_button, 'UserData', 0);
stop_status = 0;
while stop_status == 0 && index >= 1
    time0 = tic;
    stop_status = get(handles.Stop_button, 'UserData');
    simage = sim_all(:, :, :, index);
    set(handles.him, 'CData', simage);
    drawnow;
    index = index-1;
    currentindex = index+1;
    handles.currentindex = currentindex;
    if get(handles.RT_check, 'Value')
        plot_trajectory(handles);
    end
    while toc(time0) < interval
    end
    disp(['Real FrameRate is ' num2str(round(1/toc(time0)), '%d')]);
end
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
disp('Done!');
if ~get(handles.RT_check, 'Value')
    plot_trajectory(handles);
end

function FirstFrame_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = 1;
simage = sim_all(:, :, :, index);
set(handles.him, 'CData', simage);
currentindex = index;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
plot_trajectory(handles);

function SaveVideo_button_Callback(hObject, eventdata, handles)
sim_all = handles.sim_all;
FileName = handles.v.FileName;
if strcmp(handles.v.bodypart, 'left forelimb')
    if exist([FileName '_LeftForelimb_annotated.avi'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
                return;
            case 'yes'
        end
    end
end
if strcmp(handles.v.bodypart, 'right forelimb')
    if exist([FileName '_RightForelimb_annotated.avi'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
                return;
            case 'yes'
        end
    end
end
if strcmp(handles.v.bodypart, 'jaw')
    if exist([FileName '_Jaw_annotated.avi'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
                return;
            case 'yes'
        end
    end
end
FrameRate = str2double(get(handles.FrameRate_text, 'String'));
if isnan(FrameRate)
    FrameRate = handles.FrameRate;
end
if strcmp(handles.v.bodypart, 'left forelimb')
    vidObj = VideoWriter([FileName '_LeftForelimb_annotated'], 'Motion JPEG AVI');
elseif strcmp(handles.v.bodypart, 'right forelimb')
    vidObj = VideoWriter([FileName '_RightForelimb_annotated'], 'Motion JPEG AVI');
elseif strcmp(handles.v.bodypart, 'jaw')
    vidObj = VideoWriter([FileName '_Jaw_annotated'], 'Motion JPEG AVI');
end
set(vidObj, 'FrameRate', FrameRate, 'Quality', 75);
open(vidObj);
if ~get(handles.bin_check, 'Value')
    imvideo = sim_all;
    writeVideo(vidObj, imvideo);
else
    bin = str2double(get(handles.bin_text, 'String'));
    Height = size(sim_all, 1);
    Width = size(sim_all, 2);
    if mod(Height, bin) || mod(Width, bin)
        errordlg('Check your Bin !', 'ERROR');
        close(vidObj);
        return;
    end
    imvideo = single(zeros(Height/bin, Width/bin, 3, size(sim_all, 4)));
    temp = single(zeros(size(imvideo)));
    for m = 1:bin
        for n = 1:bin
            temp = temp+sim_all(m:bin:end, n:bin:end, :, :);
        end
    end
    temp = temp/bin^2;
    imvideo = temp;
    writeVideo(vidObj, imvideo);
end
close(vidObj);
clear imvideo;
disp('Done!');

function StepBack1_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
FrameRate = handles.FrameRate;
currentindex = handles.currentindex;
if currentindex == 1
    return;
end
index = currentindex-FrameRate;
if index <= 0
    index = 1;
end
simage = sim_all(:, :, :, index);
set(handles.him, 'CData', simage);
currentindex = index;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
plot_trajectory(handles);

function bin_text_Callback(hObject, eventdata, handles)

function bin_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bin_check_Callback(hObject, eventdata, handles)

function ROI_button_Callback(hObject, eventdata, handles)
sim_all = handles.sim_all;
trajectory = handles.trajectory;
skipframes = handles.skipframes;
im = imrotate(handles.currentim, -90*handles.rotate_times);
set(handles.him, 'CData', boost_IM_brightness(im, handles.thr));
[~, rect] = imcrop(handles.him);
x = rect(1, 1);
y = rect(1, 2);
w = rect(1, 3);
h = rect(1, 4);
bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(im, 1), size(im, 2));
hblob = vision.BlobAnalysis('AreaOutputPort', false, 'CentroidOutputPort', true);
centroid = step(hblob, BW);
trajectory(handles.currentindex, :) = centroid;
simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [centroid, 5], 'Color', [0 1 0], 'Opacity', 1);
if handles.v.data(handles.currentindex+skipframes, 1)
    simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
elseif skipframes+handles.currentindex < find(handles.v.data(:, 1), 1)
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
else
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
end
sim_all(:, :, :, handles.currentindex) = simage;
set(handles.him, 'CData', simage);
handles.rect = rect;
handles.sim_all = sim_all;
handles.trajectory = trajectory;
set(handles.him.Parent, 'XLim', [0.5 size(simage, 2)+0.5]);
set(handles.him.Parent, 'YLim', [0.5 size(simage, 1)+0.5]);
zoominh = findall(handles.figure1,'tag','Exploration.ZoomIn');
set(zoominh, 'State', 'on');
xu_putdowntext('zoomin', zoominh);
guidata(hObject, handles);
plot_trajectory(handles);

function plot_trajectory(handles)
lightindex = handles.v.data(:, 1);
skipframes = handles.skipframes;
trajectory = handles.trajectory;
currentindex = handles.currentindex;
lightindex = lightindex(1:skipframes+currentindex);
tj = [zeros(skipframes, 2)/0; trajectory(1:currentindex, :)];
tj(:, 1) = tj(:, 1)-tj(skipframes+1, 1);
tj(:, 2) = tj(skipframes+1, 2)-tj(:, 2);
tj_nolight = tj(lightindex == 0, :);
tj_light = tj(lightindex == 1, :);
tj1 = circshift(tj, [-1, 0]);
u = tj1(:, 1) - tj(:, 1);
v = tj1(:, 2) - tj(:, 2);
u_nolight = u(lightindex == 0);
u_light = u(lightindex == 1);
v_nolight = v(lightindex == 0);
v_light = v(lightindex == 1);
cla(handles.plot_axes);
if max(lightindex) == 0
    quiver(handles.plot_axes, tj_nolight(1:end-1, 1), tj_nolight(1:end-1, 2), u_nolight(1:end-1), v_nolight(1:end-1), 0,  'Color', 'black', 'LineWidth', 1);
else
    if currentindex+skipframes > find(lightindex == 1, 1, 'last')
        quiver(handles.plot_axes, tj_nolight(1:end-1, 1), tj_nolight(1:end-1, 2), u_nolight(1:end-1), v_nolight(1:end-1), 0,  'Color', 'black', 'LineWidth', 1);
        p1 = quiver(handles.plot_axes, tj_light(1:end-1, 1), tj_light(1:end-1, 2), u_light(1:end-1), v_light(1:end-1), 0,  'Color', 'blue', 'LineWidth', 1);
        quiver(handles.plot_axes, tj_light(end, 1), tj_light(end, 2), u_light(end), v_light(end), 0,  'Color', 'black', 'LineWidth', 1);
    else
        quiver(handles.plot_axes, tj_nolight(:, 1), tj_nolight(:, 2), u_nolight, v_nolight, 0,  'Color', 'black', 'LineWidth', 1);
        p1 = quiver(handles.plot_axes, tj_light(1:end-1, 1), tj_light(1:end-1, 2), u_light(1:end-1), v_light(1:end-1), 0,  'Color', 'blue', 'LineWidth', 1);
    end
end
p2 = plot(handles.plot_axes, tj(1, 1), tj(1, 2), 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
if currentindex+skipframes == size(handles.v.data, 1)
%     plot(handles.plot_axes, tj_nolight(2:end-1, 1), tj_nolight(2:end-1, 2), 'ok', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    p3 = plot(handles.plot_axes, tj_nolight(end, 1), tj_nolight(end, 2), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
else
%     plot(handles.plot_axes, tj_nolight(2:end, 1), tj_nolight(2:end, 2), 'ok', 'MarkerSize', 1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
end
if max(lightindex) == 1
    p4 = plot(handles.plot_axes, tj_light(1, 1), tj_light(1, 2), 'sg', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
    if currentindex+skipframes >= find(lightindex == 1, 1, 'last')
%         plot(handles.plot_axes, tj_light(2:end-1, 1), tj_light(2:end-1, 2), 'ob', 'MarkerSize', 1, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        p5 = plot(handles.plot_axes, tj_light(end, 1), tj_light(end, 2), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    else
%         plot(handles.plot_axes, tj_light(2:end, 1), tj_light(2:end, 2), 'ob', 'MarkerSize', 1, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    end
end
% if currentindex == size(trajectory, 1)
%     xrange = max(trajectory(:, 1)) - min(trajectory(:, 1));
%     yrange = max(trajectory(:, 2)) - min(trajectory(:, 2));
%     box(handles.plot_axes, 'off');
%     legend(handles.plot_axes, [p1 p2 p3 p4 p5], 'Light On', 'Video Start', 'Video End', 'Light Start', 'Light End', 'Location', 'BestOutside');
%     set(handles.plot_axes, 'XLim', [min(trajectory(:, 1))-0.05*xrange max(trajectory(:, 1))+0.05*xrange], 'YLim', [min(trajectory(:, 2))-0.05*yrange max(trajectory(:, 2))+0.05*yrange]);
% end

function Mask_button_Callback(hObject, eventdata, handles)
im = imrotate(handles.currentim, -90*handles.rotate_times);
frame = rgb2hsv(im);
frame = frame(:, :, 3);
ROI = imcrop(frame, handles.rect);
crr = xcorr2(frame-mean2(frame), ROI-mean2(frame));
figure;imshow(crr, [])
h = imfreehand('Closed', 'True');
wait(h);
BW = createMask(h);
close;
handles.Mask = BW;
guidata(hObject, handles);

function FT_radiobutton_Callback(hObject, eventdata, handles)
set(handles.Mask_button, 'Visible', 'off');

function CC_radiobutton_Callback(hObject, eventdata, handles)
set(handles.Mask_button, 'Visible', 'on');

function MTrack_button_Callback(hObject, eventdata, handles)
im_all = handles.im_all;
sim_all = handles.sim_all;
trajectory = handles.trajectory;
skipframes = handles.skipframes;
currentindex = handles.currentindex;
if currentindex == handles.Frames
    return;
end
index = currentindex+1;
if isnan(im_all(1, 1, 1, index))
    im = step(handles.v.FileReader);
    im_all(:, :, :, index) = im;
else
    im = im_all(:, :, :, index);
end
set(handles.him, 'CData', boost_IM_brightness(im, handles.thr));
[x, y] = ginput(1);
trajectory(index, :) = [x y];
simage = insertShape(boost_IM_brightness(im, handles.thr), 'FilledCircle', [x, y, 5], 'Color', [0 1 0], 'Opacity', 1);
if handles.v.data(skipframes+index, 1)
    simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Stimulation On'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
elseif skipframes+index < find(handles.v.data(:, 1), 1)
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Pre Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
else
    simage = insertText(simage, [0 size(im, 1)], [handles.mouseline ', Post Stimulation'], 'FontSize', 18, 'BoxColor', [1 1 1], 'BoxOpacity', 0.4, 'TextColor', 'white', 'AnchorPoint', 'LeftBottom');
end
sim_all(:, :, :, index) = simage;
set(handles.him, 'CData', simage);
if isDone(handles.v.FileReader)
    release(handles.v.FileReader);
end
currentindex = index;
handles.im_all = im_all;
handles.sim_all = sim_all;
handles.trajectory = trajectory;
refresh_handles(hObject, handles, trajectory, currentindex, im_all);
plot_trajectory(handles);

function RT_check_Callback(hObject, eventdata, handles)

function plot_axes_ButtonDownFcn(hObject, eventdata, handles)
hf = figure('Color', [1 1 1]);
copyobj(handles.plot_axes, hf);
set(gca, 'Position', [0.13 0.11 0.775 0.815]);

function SaveTrajec_button_Callback(hObject, eventdata, handles)
trajectory_all = zeros(size(handles.v.data, 1), 2)/0;
trajectory = handles.trajectory;
trajectory_all(1:size(trajectory, 1), :) = trajectory;
trajectory_all = circshift(trajectory_all, handles.skipframes);
trajectory = trajectory_all;
FileName = handles.v.FileName;
if strcmp(handles.v.bodypart, 'left forelimb')
    if exist([FileName '_trajectory.mat'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
            case 'yes'
                save([FileName '_trajectory.mat'], 'trajectory');
        end
    else
        save([FileName '_trajectory.mat'], 'trajectory');
    end
end
if strcmp(handles.v.bodypart, 'right forelimb')
    if exist([FileName '_rightpaw_trajectory.mat'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
            case 'yes'
                save([FileName '_rightpaw_trajectory.mat'], 'trajectory');
        end
    else
        save([FileName '_rightpaw_trajectory.mat'], 'trajectory');
    end
end
if strcmp(handles.v.bodypart, 'jaw')
    if exist([FileName '_jaw_trajectory.mat'], 'file')
        user_response = modaldlg;
        switch lower(user_response)
            case 'no'
                % take no action
            case 'yes'
                save([FileName '_jaw_trajectory.mat'], 'trajectory');
        end
    else
        save([FileName '_jaw_trajectory.mat'], 'trajectory');
    end
end
disp('Done');
