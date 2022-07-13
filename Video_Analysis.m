%%
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June. 2015
% xan@cshl.edu
%*---------------------------------------------------------------------*
[videofile, videopath] = uigetfile('*.avi', 'Select the video to be analyzed');
if videofile == 0
    return;
end
[parameterfile, parameterpath] = uigetfile('*.mat', 'Select the parameter file');
if parameterfile == 0
    return;
end

load([parameterpath parameterfile]);
pre = APre;
post = APost;
rduration = ADuration;

readerobj = VideoReader([videopath videofile]);
FrameRate = readerobj.FrameRate;
Frames = readerobj.NumberOfFrames;
centroid_all = [];

videoFileReader = vision.VideoFileReader([videopath videofile]);
videoFrame       = step(videoFileReader);
figure; imshow(videoFrame, []);
[ROI, rect] = imcrop;
close gcf;
x = rect(1, 1);
y = rect(1, 2);
w = rect(1, 3);
h = rect(1, 4);
bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'Color', [1 1 0]);

BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(videoFrame, 1), size(videoFrame, 2));
hblob = vision.BlobAnalysis('AreaOutputPort', false, 'CentroidOutputPort', true);
centroid = step(hblob, BW);
hf = figure('Position', [20+size(videoFrame, 2)+30, 100, size(videoFrame, 2), size(videoFrame, 1)]);
plot(centroid(1), centroid(2), 'sg', 'MarkerSize', 4, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
haxes = gca;
axis(haxes, 'equal');
box(haxes, 'off');
set(gca, 'YDir', 'reverse');
hold(gca, 'on');
centroid_all = centroid;

Method = 1;
switch Method
    case 1
        points = detectMinEigenFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
    case 2
        points = detectHarrisFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
    case 3
        points = detectFASTFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
    case 4
        points = detectBRISKFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
    case 5
        points = detectSURFFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
    case 6
        points = detectMSERFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
end
figure, imshow(videoFrame, []), hold on, title('Detected features in the ROI');
plot(points);
figure(hf);

pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
points = points.Location;
initialize(pointTracker, points, videoFrame);

videoPlayer  = vision.VideoPlayer('Position',...
    [20 100 [size(videoFrame, 2), size(videoFrame, 1)]+30]);

oldPoints = points;

while ~isDone(videoFileReader)
    % get the next frame
    videoFrame = step(videoFileReader);

    % Track the points. Note that some points may be lost.
    [points, isFound] = step(pointTracker, videoFrame);
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
        
        % Insert a bounding box around the object being tracked
        videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon);
                
        % Display tracked points
        videoFrame = insertMarker(videoFrame, visiblePoints, '+', ...
            'Color', 'green');
        
        BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(videoFrame, 1), size(videoFrame, 2));
        centroid = step(hblob, BW);
        videoFrame = insertShape(videoFrame, 'FilledCircle', [centroid 3], 'Color', [1 0 0]);
        plot(haxes, centroid(1), centroid(2), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        centroid_all = [centroid_all; centroid];
        
        % Reset the points
        oldPoints = visiblePoints;
        setPoints(pointTracker, oldPoints);        
    end
    
    % Display the annotated video frame using the video player object
    step(videoPlayer, videoFrame);
end
cla(haxes);
% plot(haxes, centroid_all(:, 1), centroid_all(:, 2), '-k');
% plot(haxes, centroid_all(2:end-1, 1), centroid_all(2:end-1, 2), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(haxes, centroid_all(1, 1), centroid_all(1, 2), 'sg', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
% plot(haxes, centroid_all(end, 1), centroid_all(end, 2), 'sr', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% 
% plot(haxes, centroid_all(round(pre*FrameRate), 1), centroid_all(round(pre*FrameRate), 2), 'og', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
% plot(haxes, centroid_all(round((rduration-post)*FrameRate+1), 1), centroid_all(round((rduration-post)*FrameRate+1), 2), 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

lsindex = round(pre*FrameRate);
leindex = round((rduration-post)*FrameRate)+1;
xvs2ls = centroid_all(1:lsindex, 1);
yvs2ls = centroid_all(1:lsindex, 2);
xle2ve = centroid_all(leindex:end, 1);
yle2ve = centroid_all(leindex:end, 2);
xls2le = centroid_all(lsindex:leindex, 1);
yls2le = centroid_all(lsindex:leindex, 2);
plot(haxes, xvs2ls, yvs2ls, '-k');
plot(haxes, xle2ve, yle2ve, '-k');
plot(haxes, xvs2ls(2:end-1), yvs2ls(2:end-1), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
plot(haxes, xle2ve(2:end-1), yle2ve(2:end-1), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
plot(haxes, xls2le, yls2le, '-b');
p1 = plot(haxes, xls2le(2:end-1), yls2le(2:end-1), 'ob', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b');
p2 = plot(haxes, centroid_all(1, 1), centroid_all(1, 2), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
p3 = plot(haxes, centroid_all(end, 1), centroid_all(end, 2), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
p4 = plot(haxes, xls2le(1), yls2le(1), 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
p5 = plot(haxes, xls2le(end), yls2le(end), 'sg', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
legend([p1 p2 p3 p4 p5], 'Light On', 'Video Start', 'Video End', 'Light Start', 'Light End', 'Location', 'BestOutside');
% legend(haxes, 'boxoff');
xrange = max(centroid_all(:, 1)) - min(centroid_all(:, 1));
yrange = max(centroid_all(:, 2)) - min(centroid_all(:, 2));
set(haxes, 'XLim', [min(centroid_all(:, 1))-0.05*xrange max(centroid_all(:, 1))+0.05*xrange], 'YLim', [min(centroid_all(:, 2))-0.05*yrange max(centroid_all(:, 2))+0.05*yrange]);

coord1 = centroid_all;
coord2 = circshift(centroid_all, -1);
x = coord2(:, 1) - coord1(:, 1);
y = coord2(:, 2) - coord1(:, 2);
figure;
quiver(coord1(1:lsindex-1, 1), coord1(1:lsindex-1, 2), x(1:lsindex-1), y(1:lsindex-1), 0,  'Color', 'black', 'LineWidth', 1);
hold on;
quiver(coord1(leindex:end-1, 1), coord1(leindex:end-1, 2), x(leindex:end-1), y(leindex:end-1), 0,  'Color', 'black', 'LineWidth', 1);
p1 = quiver(coord1(lsindex:leindex-1, 1), coord1((lsindex:leindex-1), 2), x(lsindex:leindex-1), y(lsindex:leindex-1), 0,  'Color', 'blue', 'LineWidth', 1);
p2 = plot(coord1(1, 1), coord1(1, 2), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
p3 = plot(coord1(end, 1), coord1(end, 2), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
p4 = plot(coord1(lsindex, 1), coord1(lsindex, 2), 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
p5 = plot(coord1(leindex, 1), coord1(leindex, 2), 'sg', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
set(gca, 'YDir', 'reverse');
axis(gca, 'equal');
box(gca, 'off');
legend([p1 p2 p3 p4 p5], 'Light On', 'Video Start', 'Video End', 'Light Start', 'Light End', 'Location', 'BestOutside');
% legend(haxes, 'boxoff');
set(gca, 'XLim', [min(centroid_all(:, 1))-0.05*xrange max(centroid_all(:, 1))+0.05*xrange], 'YLim', [min(centroid_all(:, 2))-0.05*yrange max(centroid_all(:, 2))+0.05*yrange]);

% Clean up
release(videoFileReader);
release(videoPlayer);
release(pointTracker);

implay([videopath videofile]);