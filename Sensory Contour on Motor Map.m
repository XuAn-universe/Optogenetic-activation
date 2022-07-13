%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
filtersize = 5;
h = fspecial('gaussian', filtersize, filtersize/5);
filtered_image = imfilter(FrameImage,h);
figure;
imshow(filtered_image, [mean2(filtered_image)-3*std2(filtered_image) mean2(filtered_image)+3*std2(filtered_image)]);
BW = filtered_image <= mean2(filtered_image)-2*std2(filtered_image);
% sorted_image = sort(filtered_image(:));
% forelimb_area = filtered_image <= sorted_image(round(numel(sorted_image)*0.05));
% figure;
% imshow(forelimb_area, [])
figure;
imshow(BW, [])
se = strel('disk', 3); %5
BW = imdilate(BW, se);
figure;
imshow(BW, []);
BW = imfill(BW, 'holes');
figure;
imshow(BW, []);
BW = imclearborder(BW, 8);
figure;
imshow(BW, []);
se = strel('disk', 1);
BW = imerode(BW, se);
temp = bwconncomp(BW);
while temp.NumObjects > 1
    BW = imerode(BW, se);
    temp = bwconncomp(BW);
end
figure;
imshow(BW, []);
B = bwboundaries(BW);
prop = regionprops(BW, 'Centroid');
centroid = prop.Centroid;

%% Contour on 1D maps
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;
center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
disp(center);

scalefactor = 50;
mappixel_resolution = 375/scalefactor;
ML = (ML-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
AP = (8*scalefactor+scalefactor/2+0.5)-AP*1000/mappixel_resolution;
center(1) = (center(1)-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
center(2) = (8*scalefactor+scalefactor/2+0.5)-center(2)*1000/mappixel_resolution;

hold on;
plot(ML, AP, 'k', 'LineWidth', 2);
plot(center(1), center(2), 'wo', 'MarkerFaceColor', 'w');

%% Contour on 2D vector maps
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;
center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
disp(center);

refx0 = -0.5;
refy0 = 7.5;
ML = ML*1000/375+refx0;
AP = AP*1000/375+refy0;
center(1) = center(1)*1000/375+refx0;
center(2) = center(2)*1000/375+refy0;

hold on;
plot(ML, AP, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
% plot(ML, AP, '-g', 'LineWidth', 1.5);
plot(center(1), center(2), 'ro', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
% plot(center(1), center(2), 'go', 'MarkerFaceColor', 'g');

%% Forelimb Sensory on 2D map
n = 6;
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
BW = BW1;
hold on;
for i = 1:n
    eval(['temp = BW' num2str(i) ';']);
    prop = regionprops(temp, 'Centroid');
    centroid = prop.Centroid;
    center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
    center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
    disp(center);
    
    refx0 = -0.5;
    refy0 = 7.5;
    center(1) = center(1)*1000/375+refx0;
    center(2) = center(2)*1000/375+refy0;
    
    plot(center(1), center(2), 'ro', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    BW = BW|temp;
end
B = bwboundaries(BW);
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;

refx0 = -0.5;
refy0 = 7.5;
ML = ML*1000/375+refx0;
AP = AP*1000/375+refy0;

plot(ML, AP, '-', 'Color', [1 0 0], 'LineWidth', 1);

%% Hindlimb Sensory on 2D map
n = 6;
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
HBW = HBW1;
hold on;
for i = 1:n
    eval(['temp = HBW' num2str(i) ';']);
    prop = regionprops(temp, 'Centroid');
    centroid = prop.Centroid;
    center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
    center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
    disp(center);
    
    refx0 = -0.5;
    refy0 = 7.5;
    center(1) = center(1)*1000/375+refx0;
    center(2) = center(2)*1000/375+refy0;
    
    plot(center(1), center(2), 'go', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
    HBW = HBW|temp;
end
B = bwboundaries(HBW);
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;

refx0 = -0.5;
refy0 = 7.5;
ML = ML*1000/375+refx0;
AP = AP*1000/375+refy0;

plot(ML, AP, '-', 'Color', [0 1 0], 'LineWidth', 1);

%% Contour on 1D maps
n = 6;
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
scalefactor = 50;
mappixel_resolution = 375/scalefactor;
BW = BW1;
hold on;
for i = 1:n
    eval(['temp = BW' num2str(i) ';']);
    prop = regionprops(temp, 'Centroid');
    centroid = prop.Centroid;
    center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
    center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
    disp(center);
    
    center(1) = (center(1)-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
    center(2) = (8*scalefactor+scalefactor/2+0.5)-center(2)*1000/mappixel_resolution;
    
    plot(center(1), center(2), 'wo', 'MarkerFaceColor', 'w');
    
    BW = BW|temp;
end
B = bwboundaries(BW);
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;

ML = (ML-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
AP = (8*scalefactor+scalefactor/2+0.5)-AP*1000/mappixel_resolution;

plot(ML, AP, 'w', 'LineWidth', 2);

%% Contour on 1D maps
n = 6;
CView_Bottom2bregma = 0.5;
CView_Middle2bregma = 0; %(anterior- posterior+)
camera_resolution = [1312 1082];
pixel_resolution = 2.63;
top2bregma = (camera_resolution(2)-NumData.ROI(2))*pixel_resolution/1000+CView_Bottom2bregma;
left2bregma = (camera_resolution(1)/2-NumData.ROI(1))*pixel_resolution/1000-CView_Middle2bregma;
scalefactor = 50;
mappixel_resolution = 375/scalefactor;
HBW = HBW1;
hold on;
for i = 1:n
    eval(['temp = HBW' num2str(i) ';']);
    prop = regionprops(temp, 'Centroid');
    centroid = prop.Centroid;
    center(1) = top2bregma-centroid(2)*pixel_resolution/1000;
    center(2) = left2bregma-centroid(1)*pixel_resolution/1000;
    disp(center);

    center(1) = (center(1)-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
    center(2) = (8*scalefactor+scalefactor/2+0.5)-center(2)*1000/mappixel_resolution;

    plot(center(1), center(2), 'ko', 'MarkerFaceColor', 'k');
    
    HBW = HBW|temp;
end
B = bwboundaries(HBW);
ML = B{1}(:, 1);
AP = B{1}(:, 2);
ML = top2bregma-ML*pixel_resolution/1000;
AP = left2bregma-AP*pixel_resolution/1000;

ML = (ML-0.375)*1000/mappixel_resolution+(scalefactor/2+0.5);
AP = (8*scalefactor+scalefactor/2+0.5)-AP*1000/mappixel_resolution;

plot(ML, AP, 'k', 'LineWidth', 2);
