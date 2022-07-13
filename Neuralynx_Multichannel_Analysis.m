function varargout = Neuralynx_Multichannel_Analysis(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Neuralynx_Multichannel_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @Neuralynx_Multichannel_Analysis_OutputFcn, ...
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

function Neuralynx_Multichannel_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = Neuralynx_Multichannel_Analysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function dir_button_Callback(hObject, eventdata, handles)
if ~isfield(handles, 'DataDir')
    DataDir = uigetdir('C:\', 'Select a data folder to analyze');
else
    DataDir = uigetdir(handles.DataDir, 'Select a data folder to analyze');
end
if DataDir ~= 0
    handles.DataDir = DataDir;
    set(handles.dir_text, 'String', DataDir);
    RefreshLists(hObject, eventdata, handles)
end
guidata(hObject, handles);

function RefreshLists(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
DirList = dir;
nCSC = 0;
for n = 3:size(DirList)
    if strncmp(DirList(n).name, 'CSC', 3)
        nCSC = nCSC+1;
        channels{nCSC} = DirList(n).name;
    end
end
cd(curpwd);
chaennelID = zeros(nCSC, 1);
for n = 1:nCSC
    channelID(n) = str2num(channels{n}(4:min(find(channels{n} == '.' | channels{n} == '_'))-1));
end
[~, order] = sort(channelID);
set(handles.Channel_List, 'Value', 1);
set(handles.Channel_List, 'String', channels(order));

function dir_text_Callback(hObject, eventdata, handles)

function dir_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Cell_List_Callback(hObject, eventdata, handles)

function Cell_List_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Repeat_List_Callback(hObject, eventdata, handles)

function Repeat_List_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Event_List_Callback(hObject, eventdata, handles)
global data
Events = get(handles.Event_List, 'String');
Event = Events{get(handles.Event_List, 'Value')};
EventTimestamps = [];
n = 0;
trialnumber_last = 0;
RepeatString = [];
for i = 1:length(data.EventStrings)
    if strcmp(data.EventStrings{i}, Event)
        n = n+1;
        EventTimestamps(n) = data.Timestamps_EV(i);
        trialnumber = find(data.Timestamps_startstopEV > EventTimestamps(n), 1)/2;
        if trialnumber ~= trialnumber_last
            repeat = 0;
            trialnumber_last = trialnumber;
        end
        repeat = repeat+1;
        RepeatString{n} = ['Trial' num2str(trialnumber, '%03d') ' Repeat' num2str(repeat, '%03d')];
    end
end
set(handles.Repeat_List, 'Value', []);
set(handles.Repeat_List, 'String', RepeatString);
data.Event = Event;
data.EventTimestamps = EventTimestamps;

function Event_List_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Channel_List_Callback(hObject, eventdata, handles)
global data
DataDir = handles.DataDir;
[Timestamps_EV, EventIDs, EventStrings, Header_EV] = Nlx2MatEV([DataDir '\Events.nev'], [1 1 0 0 1], 1, 1, [] );
CSCnames = get(handles.Channel_List, 'String');
CSCname = CSCnames{get(handles.Channel_List, 'Value')};
set(handles.Channel_List, 'Enable', 'off');
pause(0.01);
[Timestamps, Samples, Header] = Nlx2MatCSC([DataDir '\' CSCname], [1 0 0 0 1], 1, 1, []);
ADBitVolts = Header{17};
ADBitVolts(1:find(ADBitVolts == ' ')) = [];
ADBitVolts = str2double(ADBitVolts);
SampleFrequency = Header{15};
SampleFrequency(1:find(SampleFrequency == ' ')) = [];
SampleFrequency = str2double(SampleFrequency);
if Timestamps == 0
    set(handles.Channel_List, 'Enable', 'on');
    return;
end
packagesize = size(Samples, 1);
Events = EventStrings(Timestamps_EV >= Timestamps(1) & Timestamps_EV <= Timestamps(end)+(packagesize-1)*1/SampleFrequency*10^6 & EventIDs ~= 19);
uniqueEvents{1} = Events{1};
for i = 1:length(Events)
    n = 0;
    for j = 1:length(uniqueEvents)
        if ~strcmp(Events{i}, uniqueEvents{j})
            n = n+1;
        end
    end
    if n == j
        uniqueEvents{j+1} = Events{i};
    end
end
set(handles.Event_List, 'Value', 1);
set(handles.Event_List, 'String', uniqueEvents);
set(handles.Repeat_List, 'Value', []);
set(handles.Repeat_List, 'String', []);
set(handles.Cell_List, 'Value', []);
set(handles.Cell_List, 'String', []);

data.SampleFrequency = SampleFrequency;
data.Timestamps_EV = Timestamps_EV(Timestamps_EV >= Timestamps(1) & Timestamps_EV <= Timestamps(end)+(packagesize-1)*1/SampleFrequency*10^6 & EventIDs ~= 19);
data.Timestamps_startstopEV = Timestamps_EV(EventIDs == 19);
data.EventStrings = Events;
Timestamps = repmat(Timestamps, packagesize, 1)+(0:1:packagesize-1)'*ones(1, size(Samples, 2))*1/SampleFrequency*10^6;
data.Timestamps = Timestamps(:);
data.Samples = Samples(:);
data.ADBitVolts = ADBitVolts;
set(handles.Channel_List, 'Enable', 'on');

function Channel_List_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Post_text_Callback(hObject, eventdata, handles)

function Post_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pre_text_Callback(hObject, eventdata, handles)

function Pre_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CSCPlot_Callback(hObject, eventdata, handles)
global data
SampleFrequency = data.SampleFrequency;
ADBitVolts = data.ADBitVolts;
EventTimestamps = data.EventTimestamps;
Timestamps = data.Timestamps;
Samples = data.Samples;
pre = str2double(get(handles.Pre_text, 'String'));
if isnan(pre)
    errordlg('Please input time for Pre Event!', 'Error');
    return;
end
pre = pre*10^6;
post = str2double(get(handles.Post_text, 'String'));
if isnan(post)
    errordlg('Please input time for Post Event!', 'Error');
    return;
end
post = post*10^6;
Repeat = get(handles.Repeat_List, 'Value');
if numel(Repeat) ~= 1
    linecolors = jet(256);
    figure;
    ha = axes;
    hold(ha, 'on');
    xlabel(ha, 'Time (s)');
    ylabel(ha, 'Voltage (mV)');
end

if get(handles.notchfilter60_check, 'Value')
    notchfilter60 = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',SampleFrequency);
end

highpass = str2double(get(handles.highpass_text, 'String'));
lowpass = str2double(get(handles.lowpass_text, 'String'));
if ~isnan(highpass) && ~isnan(lowpass)
    [b, a] = ellip(4, 0.1, 40, [highpass lowpass]*2/SampleFrequency);
elseif ~isnan(lowpass)
    [b, a] = ellip(4, 0.1, 40, lowpass*2/SampleFrequency);
elseif ~isnan(highpass)
    [b, a] = ellip(4, 0.1, 40, highpass*2/SampleFrequency);
end
           
for i = 1:length(Repeat)
    Timestamps_Repeat = Timestamps(Timestamps >= EventTimestamps(Repeat(i))-pre & Timestamps <= EventTimestamps(Repeat(i))+post);
    Timestamps_Repeat_ALL(:, i) = Timestamps_Repeat-EventTimestamps(Repeat(i));
    Samples_Repeat = Samples(Timestamps >= EventTimestamps(Repeat(i))-pre & Timestamps <= EventTimestamps(Repeat(i))+post);
    
    if get(handles.notchfilter60_check, 'Value')
        Samples_Repeat = filtfilt(notchfilter60, Samples_Repeat);
    end
    if ~isnan(highpass) && ~isnan(lowpass)
        Samples_Repeat = filtfilt(b, a, Samples_Repeat);
    elseif ~isnan(lowpass)
        Samples_Repeat = filtfilt(b, a, Samples_Repeat);
    elseif ~isnan(highpass)
        Samples_Repeat = Samples_Repeat - filtfilt(b, a, Samples_Repeat);
    end
    Samples_Repeat_ALL(:, i) = Samples_Repeat;
    
    if ~get(handles.CSCAverage_check, 'Value')
        figure;
        plot((Timestamps_Repeat-EventTimestamps(Repeat(i)))/10^6, Samples_Repeat*ADBitVolts*10^3, '-k');
        box off;
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        hold on;
        yrange = ylim;
        plot([0 0], yrange, ':k');
        ylim(yrange);
        title(['Repeat ' num2str(Repeat(i))]);
    end
    if numel(Repeat) ~= 1
        plot(ha, (Timestamps_Repeat-EventTimestamps(Repeat(i)))/10^6, Samples_Repeat*ADBitVolts*10^3, '-', 'Color', linecolors(round(i/length(Repeat)*256), :));
    end
end
if numel(Repeat) ~= 1
    yrange = ylim(ha);
    plot(ha, [0 0], yrange, ':k');
    ylim(ha, yrange);
    title(ha, 'Waveform Overlay');
    
    figure;
    Timestamps_Repeat_Mean = mean(Timestamps_Repeat_ALL, 2)/10^6;
    Samples_Repeat_Mean = mean(Samples_Repeat_ALL*ADBitVolts*10^3, 2);
    Samples_Repeat_STD = std(Samples_Repeat_ALL*ADBitVolts*10^3, 0, 2);
    Samples_Repeat_SEM = Samples_Repeat_STD/sqrt(size(Samples_Repeat_ALL, 2));
    patch([Timestamps_Repeat_Mean; Timestamps_Repeat_Mean(end:-1:1)], [Samples_Repeat_Mean+Samples_Repeat_SEM; Samples_Repeat_Mean(end:-1:1)-Samples_Repeat_SEM(end:-1:1)],...
        [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    hold on;
    plot(Timestamps_Repeat_Mean, Samples_Repeat_Mean, '-k');
    box off;
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    yrange = ylim;
    plot([0 0], yrange, ':k');
    ylim(yrange);
    title('Waveform Average with SEM');
    
    peaktstart = str2double(get(handles.peaktstart_text, 'String'))/10^3;
    peaktend = str2double(get(handles.peaktend_text, 'String'))/10^3;
    MeanAmpPre = mean(Samples_Repeat_Mean(Timestamps_Repeat_Mean <= 0))*10^3;
    PeakAmp = max(Samples_Repeat_Mean(Timestamps_Repeat_Mean >= peaktstart & Timestamps_Repeat_Mean <= peaktend))*10^3;
    disp(['Mean Amplitude for Pre = ' num2str(MeanAmpPre, '%.2f') ' micro V']);
    if ~isnan(peaktstart) && ~isnan(peaktend)
        disp(['Peak Amplitude = ' num2str(PeakAmp, '%.2f') ' micro V']);
        disp(['Peak Amplitude-Mean Amplitude = ' num2str(PeakAmp-MeanAmpPre, '%.2f') ' micro V']);
    end
end

function CSCPlot_full_Callback(hObject, eventdata, handles)
global data
SampleFrequency = data.SampleFrequency;
ADBitVolts = data.ADBitVolts;
Timestamps = data.Timestamps;
Samples = data.Samples;
highpass = str2double(get(handles.highpass_text, 'String'));
lowpass = str2double(get(handles.lowpass_text, 'String'));
if ~isnan(highpass) && ~isnan(lowpass)
    [b, a] = ellip(4, 0.1, 40, [highpass lowpass]*2/SampleFrequency);
    Samples = filtfilt(b, a, Samples);
elseif ~isnan(lowpass)
    [b, a] = ellip(4, 0.1, 40, lowpass*2/SampleFrequency);
    Samples = filtfilt(b, a, Samples);
elseif ~isnan(highpass)
    [b, a] = ellip(4, 0.1, 40, highpass*2/SampleFrequency);
    Samples = Samples - filtfilt(b, a, Samples);
end
if get(handles.notchfilter60_check, 'Value')
    notchfilter60 = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',SampleFrequency);
    Samples = filtfilt(notchfilter60, Samples);
end
figure;
plot(Timestamps/10^6, Samples*ADBitVolts*10^3, '-k');
box off;
xlabel('Time (s)');
ylabel('Voltage (mV)');
xlim([min(Timestamps)/10^6-1/SampleFrequency max(Timestamps)/10^6+1/SampleFrequency]);

function CloseALL_Callback(hObject, eventdata, handles)
set(handles.figure1, 'HandleVisibility', 'off');
close all;
set(handles.figure1, 'HandleVisibility', 'on');

function CSCAverage_check_Callback(hObject, eventdata, handles)

function SpikeSorting_Callback(hObject, eventdata, handles)
DataDir = handles.DataDir;
CSCnames = get(handles.Channel_List, 'String');
CSCname = CSCnames{get(handles.Channel_List, 'Value')};
CSC.pathname = [DataDir '\'];
CSC.filename = CSCname;
set(handles.figure1, 'Visible', 'off');
wave_clus(CSC);
set(handles.figure1, 'Visible', 'on');
set(handles.Cell_List, 'Value', []);
set(handles.Cell_List, 'String', []);

function notchfilter60_check_Callback(hObject, eventdata, handles)

function SpikeWaveform_Callback(hObject, eventdata, handles)
global data
if isfield(data, 'SpikeSorting')
    cellid = get(handles.Cell_List, 'Value');
    SampleFrequency = data.SampleFrequency;
    cluster_class = data.SpikeSorting.cluster_class;
    spikes = data.SpikeSorting.spikes;
    waves = spikes(cluster_class == cellid, :);
    figure;
    plot((1:size(waves, 2))/SampleFrequency*10^3, waves.', '-k');
    box off;
    xlim([1 size(waves, 2)]/SampleFrequency*10^3);
    xlabel('Time (ms)');
    ylabel('Voltage (micro V)');
    axis square;
    title(['N = ' num2str(size(waves, 1)) ' spikes']);
end

function ScatterPlot_Callback(hObject, eventdata, handles)
global data
if ~isfield(data, 'SpikeSorting')
    return;
end
EventTimestamps = data.EventTimestamps;
EventTimestamps = EventTimestamps/10^3;
cluster_class = data.SpikeSorting.cluster_class;
cellid = get(handles.Cell_List, 'Value');
SpikeTimestamps = cluster_class(cluster_class(:, 1) == cellid, 2);
pre = str2double(get(handles.Pre_text, 'String'));
if isnan(pre)
    errordlg('Please input time for Pre Event!', 'Error');
    return;
end
pre = pre*10^3;
post = str2double(get(handles.Post_text, 'String'));
if isnan(post)
    errordlg('Please input time for Post Event!', 'Error');
    return;
end
post = post*10^3;
bin = str2double(get(handles.Bin_text, 'String'));
if isnan(bin)
    errordlg('Please input time for Bin!', 'Error');
    return;
end
Repeat = get(handles.Repeat_List, 'Value');
figure;
axis on;
hold on;
for i = 1:length(Repeat)
    Timestamps_Repeat = SpikeTimestamps(SpikeTimestamps >= EventTimestamps(Repeat(i))-pre & SpikeTimestamps <= EventTimestamps(Repeat(i))+post);
    Timestamps_Repeat = Timestamps_Repeat-EventTimestamps(Repeat(i));
    line([Timestamps_Repeat.'/10^3; Timestamps_Repeat.'/10^3], [ones(1, numel(Timestamps_Repeat))*(i-1); ones(1, numel(Timestamps_Repeat))*i], 'Color', 'k');
    
    edges = EventTimestamps(Repeat(i))-pre:bin:EventTimestamps(Repeat(i))+post;
    spikerate(i, :) = histcounts(SpikeTimestamps, edges)/bin*10^3;
end
xlim([-pre/10^3 post/10^3]);
xlabel('Time (s)');
ylabel('Trial');
line([0 0], [0 length(Repeat)], 'LineStyle', ':', 'Color', 'k');
if numel(Repeat) ~= 1
    figure;
    spikerate_Mean = mean(spikerate);
    spikerate_STD = std(spikerate);
    spikerate_SEM = spikerate_STD/sqrt(size(spikerate, 1));
    patch([-pre+bin/2:bin:post-bin/2 fliplr(-pre+bin/2:bin:post-bin/2)], [spikerate_Mean+spikerate_SEM spikerate_Mean(end:-1:1)-spikerate_SEM(end:-1:1)],...
        [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    hold on;
    plot(-pre+bin/2:bin:post-bin/2, spikerate_Mean, '-k');
    box off;
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    yrange = ylim;
    plot([0 0], yrange, ':k');
    ylim(yrange);
    title('Mean Firing Rate with SEM');
end

function UpdateCell_Callback(hObject, eventdata, handles)
global data
DataDir = handles.DataDir;
CSCnames = get(handles.Channel_List, 'String');
CSCname = CSCnames{get(handles.Channel_List, 'Value')};
try
    data.SpikeSorting = load([DataDir '\times_' CSCname(1:end-4) '.mat']);
    cellnumber = max(data.SpikeSorting.cluster_class(:, 1));
    set(handles.Cell_List, 'Value', 1);
    set(handles.Cell_List, 'String', num2str((1:cellnumber)'));
end

function Bin_text_Callback(hObject, eventdata, handles)

function Bin_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function highpass_text_Callback(hObject, eventdata, handles)

function highpass_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lowpass_text_Callback(hObject, eventdata, handles)

function lowpass_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function peaktstart_text_Callback(hObject, eventdata, handles)

function peaktstart_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function peaktend_text_Callback(hObject, eventdata, handles)

function peaktend_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
