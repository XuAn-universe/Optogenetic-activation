%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% Grip Force
VperG = 10/600;
fid = fopen('G:\Free-moving Feeding Data\032520_EmxttA-cannulaR6_GripForce\Session003.bin','r');
[data,count] = fread(fid,[2,inf],'double');
fclose(fid);
timestamps = data(1, :);
voltage = data(2, :);
SampleRate = 1000;
freq = 10;
[b, a] = butter(6, freq/(SampleRate/2), 'low');
voltage_filtered = filtfilt(b, a, voltage);

figure;
plot(timestamps, voltage/VperG, '-k');
hold on;
plot(timestamps, voltage_filtered/VperG, '-r');
xlabel('Time (s)');
ylabel('gram');
box off;

figure;
plot(timestamps, voltage_filtered/VperG, '-r');
xlabel('Time (s)');
ylabel('gram');
box off;

%%
gca_findpeak(timestamps, voltage_filtered/VperG);

%% Bite Force
VperG = 10/2000;
fid = fopen('G:\Free-moving Feeding Data\032520_EmxttA-cannulaR6_BiteForce\Session004.bin','r');
[data,count] = fread(fid,[2,inf],'double');
fclose(fid);
timestamps = data(1, :);
voltage = data(2, :);

figure;
plot(timestamps, voltage/VperG, '-k');
xlabel('Time (s)');
ylabel('gram');
box off;

%%
gca_findpeak(timestamps, voltage/VperG);