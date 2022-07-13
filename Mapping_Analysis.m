function varargout = Mapping_Analysis(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June. 2015
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mapping_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @Mapping_Analysis_OutputFcn, ...
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

function Mapping_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = Mapping_Analysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function dir_button_Callback(hObject, eventdata, handles)
try
    DataDir = uigetdir(get(handles.dir_text, 'String'), 'Select a animal folder to analysis');
catch
    DataDir = uigetdir('C:\', 'Select a animal folder to analysis');
end
if DataDir ~= 0
    handles.DataDir = DataDir;
    set(handles.dir_text, 'String', DataDir);
    RefreshLists(hObject, eventdata, handles);
    if exist([DataDir '\ReferencePoint.mat'], 'file')
        load([DataDir '\ReferencePoint.mat']);
        set(handles.PG1Y0_text, 'String', num2str(Rf.PG1Y0));
        set(handles.PG1Z0_text, 'String', num2str(Rf.PG1Z0));
        set(handles.PG2X0_text, 'String', num2str(Rf.PG2X0));
        set(handles.PG2Z0_text, 'String', num2str(Rf.PG2Z0));
    else
        set(handles.PG1Y0_text, 'String', '');
        set(handles.PG1Z0_text, 'String', '');
        set(handles.PG2X0_text, 'String', '');
        set(handles.PG2Z0_text, 'String', '');
    end
    set(handles.ProcessTJ_button, 'UserData', []);
end

function RefreshLists(hObject, eventdata, handles)
handles.DataDir = get(handles.dir_text, 'String');
LastExpFolder = FindLastExpFolder(handles.DataDir);
ExpNumber = str2double(LastExpFolder(4:6));
set(handles.Exp_list, 'Value', []);
set(handles.Exp_list, 'String', num2str([1:ExpNumber]'));

function dir_text_Callback(hObject, eventdata, handles)

function dir_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Close_All_button_Callback(hObject, eventdata, handles)
set(handles.figure1, 'HandleVisibility', 'off');
try
    set(handles.TrackingGUI_Copy.figure1, 'HandleVisibility', 'off');
    set(handles.TrackingGUI.figure1, 'HandleVisibility', 'off');
end
close all;
set(handles.figure1, 'HandleVisibility', 'on');
try
    set(handles.TrackingGUI_Copy.figure1, 'HandleVisibility', 'on');
    set(handles.TrackingGUI.figure1, 'HandleVisibility', 'on');
end

function ALeft_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
handles.parameter_location = [DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat'];
load(handles.parameter_location);
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    handles.RepeatNumber = RepeatNumber(i);
    try
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mat']);
        yleft = y(:, 1:3);
        lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat'];
        Compute_Accelerometer(x, yleft, handles, lightpath, SP);
%         if get(handles.Save_check, 'Value')
%             saveas(hfamplitude, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Accelerometer_Amplitude(L)'], 'fig');
%             if ~isempty(hfangle)
%                 saveas(hfangle, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Accelerometer_Direction(L)'], 'fig');
%             end
%         end
    catch
        errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function  [sdelay, threshold] = Light_Delay(x, y, APre, SampleRate, SP)
if  isfield(SP, 'PulseWidth')
    threshold = 2.5;
else
    bmean = mean(y(x < APre));
    bsd = std(y(x < APre));
    threshold = bmean+2*bsd;
end
temp = find(y > threshold);
temp = temp(temp > APre*SampleRate);
if ~isempty(temp)
    startpointn = temp(1);
    sdelay = (x(startpointn)-APre)*1000;
else
    sdelay = 0;
end
    
function Compute_Accelerometer(x, y, handles, lightpath, SP)
load(handles.parameter_location);
sdelay = str2double(get(handles.SDelay_text, 'String'));
if isnan(sdelay)
    try
        ldata = load(lightpath);
        if isfield(SP, 'PulseWidth')
            [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
        else
            [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
        end
    catch
        sdelay = 0;
    end
    if sdelay == 0
        warndlg('Light Delay cannot be calculated, so is set to 0.', ['Repeat' num2str(handles.RepeatNumber)]);
    end
end
voltage = str2double(get(handles.Voltage_text, 'String'));
sensitivity = str2double(get(handles.Sensitivity_text, 'String'));
g = (y(:, 1:3)-voltage/2)/sensitivity;
g = -g;

xphi = atan2d(g(:, 2), g(:, 3));
ytheta = atand(-g(:, 1)./sqrt(g(:, 2).^2+g(:, 3).^2));
figure('Color', [1 1 1]);
subplot(3, 3, 1);
hp = plot(x, xphi, '-k');
hold on;
yrange = max(xphi)-min(xphi);
ylim = [min(xphi)-0.05*yrange max(xphi)+0.05*yrange];
plot([APre+sdelay/1000 APre+sdelay/1000], ylim, ':k');
plot([ADuration-APost+sdelay/1000 ADuration-APost+sdelay/1000], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'phi');
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

subplot(3, 3, 2);
hp = plot(x, ytheta, '-k');
hold on;
yrange = max(ytheta)-min(ytheta);
ylim = [min(ytheta)-0.05*yrange max(ytheta)+0.05*yrange];
plot([APre+sdelay/1000 APre+sdelay/1000], ylim, ':k');
plot([ADuration-APost+sdelay/1000 ADuration-APost+sdelay/1000], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'theta');
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

npre = round((sdelay/1000+APre)*SampleRate);
xphi0 = mean(xphi(1:npre));
ytheta0 = mean(ytheta(1:npre));
gstandard = zeros(size(g));
for i = 1:length(x)
    gstandard(i, :) = ([cosd(ytheta0) 0 -sind(ytheta0); 0 1 0; sind(ytheta0) 0 cosd(ytheta0)]^(-1)*[1 0 0; 0 cosd(xphi0) sind(xphi0); 0 -sind(xphi0) cosd(xphi0)]^(-1)*[g(i, 1); g(i, 2); g(i, 3)])';
end
g0 = mean(gstandard(1:npre, 3));
gstandard(:, 3) = gstandard(:, 3)-g0;
subplot(3, 3, 5);
hp = plot(x, gstandard(:, 1), '-r', x, gstandard(:, 2), '-g', x, gstandard(:, 3), '-b');
hold on;
yrange = max(max(gstandard)) - min(min(gstandard));
ylim = [min(min(gstandard))-0.05*yrange max(max(gstandard))+0.05*yrange];
plot([APre+sdelay/1000 APre+sdelay/1000], ylim, ':k');
plot([ADuration-APost+sdelay/1000 ADuration-APost+sdelay/1000], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'x', 'y', 'z');
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

amplitude = sqrt(sum(g.^2, 2));
threshold = 5*std(amplitude(1:npre));
amean = mean(amplitude(1:npre));
temp = find(amplitude > amean+threshold);
temp = temp(temp > npre);
if ~isempty(temp)
    startpointp = temp(1);
else
    startpointp = [];
end
temp = find(amplitude < amean-threshold);
temp = temp(temp > npre);
if ~isempty(temp)
    startpointn = temp(1);
else
    startpointn = [];
end
if ~isempty(startpointp) && ~isempty(startpointn)
    startpoint = min(startpointp, startpointn);
elseif ~isempty(startpointp) && isempty(startpointn)
    startpoint = startpointp;
elseif isempty(startpointp) && ~isempty(startpointn)
    startpoint = startpointn;
else
    startpoint = [];
end

subplot(3, 3, 4);
hp = plot(x, amplitude, '-k');
hold on;
plot([x(1) x(end)], [amean+threshold amean+threshold], ':k');
plot([x(1) x(end)], [amean-threshold amean-threshold], ':k');
yrange = max(amplitude) - min(amplitude);
ylim = [min(amplitude)-0.05*yrange max(amplitude)+0.05*yrange];
plot([APre+sdelay/1000 APre+sdelay/1000], ylim, ':k');
plot([ADuration-APost+sdelay/1000 ADuration-APost+sdelay/1000], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'Amplitude');
if ~isempty(startpoint)
    title(['Movement starts at ' num2str((x(startpoint)-APre-sdelay/1000)*1000) ' ms']);
end
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

Intperiod = str2double(get(handles.Integration_text, 'String'));
if ~isempty(startpoint)
    gx = gstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 1);
    gy = gstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 2);
    gz = gstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 3);
    subplot(3, 3, 6);
    hp1 = plot(gx, '-r');
    hold on;
    hp2 = plot(gy, '-g');
    hp3 = plot(gz, '-b');
    yrange = max([gx; gy; gz]) - min([gx; gy; gz]);
    ylim = [min([gx; gy; gz])-0.05*yrange max([gx; gy; gz])+0.05*yrange];
    set(gca, 'YLim', ylim);
    box off;
    ylabel('Acceleration (g)');
    legend([hp1 hp2 hp3], 'x', 'y', 'z');
    gx = sum(gx);
    gy = sum(gy);
    gz = sum(gz);
    ahorizontal = atan2d(gy, gx);
    aelevation = atand(gz/sqrt(gx^2+gy^2));
    title(['The angle of horizontal is ' num2str(ahorizontal) ' degree and elevation is ' num2str(aelevation) ' degree']);
    set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
end

signal = amplitude(1:npre);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.999
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
subplot(3, 3, 7);
plot(f(1:i), MX(1:i), '-k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

nlightend = round((ADuration-APost+sdelay/1000)*SampleRate);
signal = amplitude(npre+1:nlightend);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.999
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
subplot(3, 3, 8);
plot(f(1:i), MX(1:i), '-k');
hold on;
if SP.Frequency == 60
    SP.Frequency = 120;
end
plot([SP.Frequency SP.Frequency], [min(MX(1:i)) max(MX(1:i))], ':k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);
set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));

function myplot(hObject, eventdata, ha)
hf = figure('Color', [1 1 1]);
copyobj(ha, hf);
set(gca, 'Position', [0.13 0.11 0.775 0.815]);

function VideoPlay_button_Callback(hObject, eventdata, handles)
handles.VideoSource = 'Point Grey';
PlayVideo(hObject, eventdata, handles);

function PlayVideo(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
try
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
catch
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\exp' num2str(ExpNumber, '%03d') '.mat']);
    if isfield(StimulusParameter, 'Order')
        SP.Order = StimulusParameter.Order;
    end
end
for i = 1:length(RepeatNumber)
    if isfield(StimulusParameter, 'Order')
        TrialNumber = find(SP.Order == SpotNumber);
        TrialNumber = TrialNumber(RepeatNumber(i));
    else
        TrialNumber = RepeatNumber(i);
    end
    flag = 0;
    switch handles.VideoSource
        case 'Point Grey'
            try
                mh = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mp4']);
                text = ['exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d')];
                mh.Parent.Name = text;
                flag = 1;
            end
            try
                mh = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.avi']);
                text = ['exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d')];
                mh.Parent.Name = text;
                flag = 1;
            end
            try
                parametertext = ['Laser Knob: ' num2str(SP.LaserKnob, '%.2f')...
                    '; Pulse Width: ' num2str(SP.PulseWidth) ' ms; Frequency: ' num2str(SP.Frequency)...
                    ' Hz; Duration: ' num2str(SP.Duration) ' s'];
                mh1 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\C1_trial' num2str(TrialNumber, '%03d') '.avi']);
                text1 = ['exp' num2str(ExpNumber, '%03d') '\C1_trial' num2str(TrialNumber, '%03d') ' | ' parametertext];
                mh1.Parent.Name = text1;
                mh2 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\C2_trial' num2str(TrialNumber, '%03d') '.avi']);
                text2 = ['exp' num2str(ExpNumber, '%03d') '\C2_trial' num2str(TrialNumber, '%03d') ' | ' parametertext];
                mh2.Parent.Name = text2;
                flag = 1;
            end
            try
                parametertext = ['Laser Knob: ' num2str(StimulusParameter.LaserKnob, '%.2f')...
                    '; Pulse Width: ' num2str(StimulusParameter.PulseWidth) ' ms; Frequency: ' num2str(StimulusParameter.Frequency)...
                    ' Hz; Duration: ' num2str(StimulusParameter.Duration) ' s'];
                mh1 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
                text1 = ['exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber) ' | ' parametertext];
                mh1.Parent.Name = text1;
                mh1.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
                mh2 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
                text2 = ['exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) ' | ' parametertext];
                mh2.Parent.Name = text2;
                mh2.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
                try
                    mh3 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\LGwebcam' num2str(TrialNumber) '.avi']);
                    text3 = ['exp' num2str(ExpNumber, '%03d') '\LGwebcam' num2str(TrialNumber) ' | ' parametertext];
                    mh3.Parent.Name = text3;
                    mh3.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
                end
                flag = 1;
            end
        case 'LG Webcam'
            try
                mh = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.mp4']);
                text = ['exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d')];
                mh.Parent.Name = text;
                flag = 1;
            end
            try
                mh = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.avi']);
                text = ['exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d')];
                mh.Parent.Name = text;
                flag = 1;
            end
    end
    if flag == 0
        errordlg(['Cannot play exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function VideoAnalysis_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
elseif length(RepeatNumber) ~= 1
    errordlg('You can only choose one repeat !', 'ERROR');
    return;
end
FramesBL = str2double(get(handles.FramesBL_text, 'String'));
FramesAL = str2double(get(handles.FramesAL_text, 'String'));
try
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
catch
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\exp' num2str(ExpNumber, '%03d') '.mat']);
    if isfield(StimulusParameter, 'Order')
        SP.Order = StimulusParameter.Order;
    end
end
for i = 1:length(RepeatNumber)
    if isfield(StimulusParameter, 'Order')
        TrialNumber = find(SP.Order == SpotNumber);
        TrialNumber = TrialNumber(RepeatNumber(i));
    else
        TrialNumber = RepeatNumber(i);
    end
    v.data = load([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
    v.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
    v.FileName = [DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber)];
    v.FrameRange = [FramesBL FramesAL];
    if get(handles.LForelimb_radiobutton, 'Value')
        v.bodypart = 'left forelimb';
    elseif get(handles.RForelimb_radiobutton, 'Value')
        v.bodypart = 'right forelimb';
    elseif get(handles.Jaw_radiobutton, 'Value')
        v.bodypart = 'jaw';
    end
    handles.TrackingGUI_Copy = Object_Tracking_Copy(v);
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        v.data = load([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
        v.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
        v.FileName = [DataDir '\exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber)];
        v.FrameRange = [FramesBL FramesAL];
        if get(handles.LForelimb_radiobutton, 'Value')
            v.bodypart = 'left forelimb';
        elseif get(handles.RForelimb_radiobutton, 'Value')
            v.bodypart = 'right forelimb';
        end
        handles.TrackingGUI = Object_Tracking(v);
    end
end   
guidata(hObject, handles);

% load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
% pre = APre;
% post = APost;
% rduration = ADuration;
% for i = 1:length(RepeatNumber)
%     TrialNumber = find(SP.Order == SpotNumber);
%     TrialNumber = TrialNumber(RepeatNumber(i));
%     try
%         lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat'];
%         sdelay = str2double(get(handles.SDelay_text, 'String'));
%         if isnan(sdelay)
%             try
%                 ldata = load(lightpath);
%                 if isfield(SP, 'PulseWidth')
%                     [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
%                 else
%                     [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
%                 end
%             catch
%                 sdelay = 0;
%             end
%             if sdelay == 0
%                 warndlg('Light Delay cannot be calculated, so is set to 0.', ['Repeat' num2str(RepeatNumber(i))]);
%             end
%         end
%         load([DataDir '\exp' num2str(ExpNumber, '%03d') '\timestamp' num2str(TrialNumber, '%03d') '.mat']);
%         
%         try
%             readerobj = VideoReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mp4']);
%         catch
%             readerobj = VideoReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.avi']);
%         end
%         FrameRate = readerobj.FrameRate;
%         Frames = readerobj.NumberOfFrames;
%         centroid_all = [];
%         
%         try
%             videoFileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mp4']);
%         catch
%             videoFileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.avi']);
%         end
%         videoFrame = step(videoFileReader);
%         figure; imshow(videoFrame, []);
%         [ROI, rect] = imcrop;
%         close gcf;
%         x = rect(1, 1);
%         y = rect(1, 2);
%         w = rect(1, 3);
%         h = rect(1, 4);
%         bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
%         videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'Color', [1 1 0]);
%         
%         BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(videoFrame, 1), size(videoFrame, 2));
%         hblob = vision.BlobAnalysis('AreaOutputPort', false, 'CentroidOutputPort', true);
%         centroid = step(hblob, BW);
%         hf = figure('Position', [20+size(videoFrame, 2)+30, 100, size(videoFrame, 2), size(videoFrame, 1)]);
%         plot(centroid(1), centroid(2), 'sg', 'MarkerSize', 4, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         haxes = gca;
%         axis(haxes, 'equal');
%         box(haxes, 'off');
%         set(gca, 'YDir', 'reverse');
%         hold(gca, 'on');
%         centroid_all = centroid;
%         
%         Method = 6;
%         switch Method
%             case 1
%                 points = detectMinEigenFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%             case 2
%                 points = detectHarrisFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%             case 3
%                 points = detectFASTFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%             case 4
%                 points = detectBRISKFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%             case 5
%                 points = detectSURFFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%             case 6
%                 points = detectMSERFeatures(rgb2gray(videoFrame), 'ROI', ceil(rect));
%         end
%         figure, imshow(videoFrame, []), hold on, title('Detected features in the ROI');
%         plot(points);
%         figure(hf);
%         
%         pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
%         points = points.Location;
%         initialize(pointTracker, points, videoFrame);
%         
%         videoPlayer  = vision.VideoPlayer('Position',...
%             [20 100 [size(videoFrame, 2), size(videoFrame, 1)]+30]);
%         
%         oldPoints = points;
%         
%         while ~isDone(videoFileReader)
%             % get the next frame
%             videoFrame = step(videoFileReader);
%             
%             % Track the points. Note that some points may be lost.
%             [points, isFound] = step(pointTracker, videoFrame);
%             visiblePoints = points(isFound, :);
%             oldInliers = oldPoints(isFound, :);
%             
%             if size(visiblePoints, 1) >= 2 % need at least 2 points
%                 
%                 % Estimate the geometric transformation between the old points
%                 % and the new points and eliminate outliers
%                 [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
%                     oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
%                 
%                 % Apply the transformation to the bounding box
%                 [bboxPolygon(1:2:end), bboxPolygon(2:2:end)] ...
%                     = transformPointsForward(xform, bboxPolygon(1:2:end), bboxPolygon(2:2:end));
%                 
%                 % Insert a bounding box around the object being tracked
%                 videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon);
%                 
%                 % Display tracked points
%                 videoFrame = insertMarker(videoFrame, visiblePoints, '+', ...
%                     'Color', 'green');
%                 
%                 BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(videoFrame, 1), size(videoFrame, 2));
%                 centroid = step(hblob, BW);
%                 videoFrame = insertShape(videoFrame, 'FilledCircle', [centroid 3], 'Color', [1 0 0]);
%                 plot(haxes, centroid(1), centroid(2), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%                 centroid_all = [centroid_all; centroid];
%                 
%                 % Reset the points
%                 oldPoints = visiblePoints;
%                 setPoints(pointTracker, oldPoints);
%             end
%             
%             % Display the annotated video frame using the video player object
%             step(videoPlayer, videoFrame);
%         end
%         cla(haxes);
%         % plot(haxes, centroid_all(:, 1), centroid_all(:, 2), '-k');
%         % plot(haxes, centroid_all(2:end-1, 1), centroid_all(2:end-1, 2), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%         % plot(haxes, centroid_all(1, 1), centroid_all(1, 2), 'sg', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         % plot(haxes, centroid_all(end, 1), centroid_all(end, 2), 'sr', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         %
%         % plot(haxes, centroid_all(round(pre*FrameRate), 1), centroid_all(round(pre*FrameRate), 2), 'og', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         % plot(haxes, centroid_all(round((rduration-post)*FrameRate+1), 1), centroid_all(round((rduration-post)*FrameRate+1), 2), 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         
%         lsindex = find(timeStamp > pre+sdelay/1000, 1 );
%         leindex = find(timeStamp <= rduration-post+sdelay/1000, 1, 'last' );
%         videodata.trajectory = centroid_all;
%         videodata.lightframeid = [lsindex leindex];
%         xvs2ls = centroid_all(1:lsindex, 1);
%         yvs2ls = centroid_all(1:lsindex, 2);
%         xle2ve = centroid_all(leindex:end, 1);
%         yle2ve = centroid_all(leindex:end, 2);
%         xls2le = centroid_all(lsindex:leindex, 1);
%         yls2le = centroid_all(lsindex:leindex, 2);
%         plot(haxes, xvs2ls, yvs2ls, '-k');
%         plot(haxes, xle2ve, yle2ve, '-k');
%         plot(haxes, xvs2ls(2:end-1), yvs2ls(2:end-1), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
%         plot(haxes, xle2ve(2:end-1), yle2ve(2:end-1), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
%         plot(haxes, xls2le, yls2le, '-b');
%         p1 = plot(haxes, xls2le(2:end-1), yls2le(2:end-1), 'ob', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b');
%         p2 = plot(haxes, centroid_all(1, 1), centroid_all(1, 2), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         p3 = plot(haxes, centroid_all(end, 1), centroid_all(end, 2), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         p4 = plot(haxes, xls2le(1), yls2le(1), 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         p5 = plot(haxes, xls2le(end), yls2le(end), 'sg', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         legend([p1 p2 p3 p4 p5], 'Light On', 'Video Start', 'Video End', 'Light Start', 'Light End', 'Location', 'BestOutside');
%         % legend(haxes, 'boxoff');
%         xrange = max(centroid_all(:, 1)) - min(centroid_all(:, 1));
%         yrange = max(centroid_all(:, 2)) - min(centroid_all(:, 2));
%         set(haxes, 'XLim', [min(centroid_all(:, 1))-0.05*xrange max(centroid_all(:, 1))+0.05*xrange], 'YLim', [min(centroid_all(:, 2))-0.05*yrange max(centroid_all(:, 2))+0.05*yrange]);
%         distance = sum(sqrt(diff(xls2le).^2+diff(yls2le).^2));
%         title(haxes, ['The total distance travelled during light on is ' num2str(distance) ' pixels']);
%         
%         coord1 = centroid_all;
%         coord2 = circshift(centroid_all, -1);
%         x = coord2(:, 1) - coord1(:, 1);
%         y = coord2(:, 2) - coord1(:, 2);
%         hftrajectory = figure;
%         quiver(coord1(1:lsindex-1, 1), coord1(1:lsindex-1, 2), x(1:lsindex-1), y(1:lsindex-1), 0,  'Color', 'black', 'LineWidth', 1);
%         hold on;
%         quiver(coord1(leindex:end-1, 1), coord1(leindex:end-1, 2), x(leindex:end-1), y(leindex:end-1), 0,  'Color', 'black', 'LineWidth', 1);
%         p1 = quiver(coord1(lsindex:leindex-1, 1), coord1((lsindex:leindex-1), 2), x(lsindex:leindex-1), y(lsindex:leindex-1), 0,  'Color', 'blue', 'LineWidth', 1);
%         p2 = plot(coord1(1, 1), coord1(1, 2), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         p3 = plot(coord1(end, 1), coord1(end, 2), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%         p4 = plot(coord1(lsindex, 1), coord1(lsindex, 2), 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         p5 = plot(coord1(leindex, 1), coord1(leindex, 2), 'sg', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         set(gca, 'YDir', 'reverse');
%         axis(gca, 'equal');
%         box(gca, 'off');
%         legend([p1 p2 p3 p4 p5], 'Light On', 'Video Start', 'Video End', 'Light Start', 'Light End', 'Location', 'BestOutside');
%         % legend(haxes, 'boxoff');
%         set(gca, 'XLim', [min(centroid_all(:, 1))-0.05*xrange max(centroid_all(:, 1))+0.05*xrange], 'YLim', [min(centroid_all(:, 2))-0.05*yrange max(centroid_all(:, 2))+0.05*yrange]);
%         
%         if get(handles.Save_check, 'Value')
%             saveas(hf, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Moving_Points'], 'fig');
%             saveas(hftrajectory, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Moving_Trajectory'], 'fig');
%         end
%         save([DataDir '\exp' num2str(ExpNumber, '%03d') '\VideoData' num2str(TrialNumber, '%03d') '.mat'], 'videodata');
%         
%         % Clean up
%         release(videoFileReader);
%         release(videoPlayer);
%         release(pointTracker);
%     catch
%         errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
%     end
% end

function ARight_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
handles.parameter_location = [DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat'];
load(handles.parameter_location);
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    handles.RepeatNumber = RepeatNumber(i);
    try
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mat']);
        yright = y(:, 4:6);
        lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat'];
        Compute_Accelerometer(x, yright, handles, lightpath, SP);
%         if get(handles.Save_check, 'Value')
%             saveas(hfamplitude, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Accelerometer_Amplitude(R)'], 'fig');
%             if ~isempty(hfangle)
%                 saveas(hfangle, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Accelerometer_Direction(R)'], 'fig');
%             end
%         end
    catch
        errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function Repeat_list_Callback(hObject, eventdata, handles)

function Repeat_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Exp_list_Callback(hObject, eventdata, handles)
set(handles.Spot_list, 'Value', []);
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = 0;
for i = 1:length(ExpNumber)
    ExpFolder = ['exp' num2str(ExpNumber(i), '%03d')];
    curpwd = pwd;
    cd([DataDir '\' ExpFolder]);
    try
        load VRDLPParameter.mat;
        cd(curpwd);
        SpotNumber = max([max(SP.Order) SpotNumber]);
    catch
        load(['exp' num2str(ExpNumber(i), '%03d') '.mat']);
        cd(curpwd);
        if isfield(StimulusParameter, 'Order')
            SpotNumber = max([max(StimulusParameter.Order) SpotNumber]);
        else
            SpotNumber = 1;
        end
    end
end
set(handles.Spot_list, 'String', num2str([1:SpotNumber]'));

function Exp_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Spot_list_Callback(hObject, eventdata, handles)
set(handles.Repeat_list, 'Value', []);
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
RepeatNumber = 0;
for i = 1:length(ExpNumber)
    ExpFolder = ['exp' num2str(ExpNumber(i), '%03d')];
    curpwd = pwd;
    cd([DataDir '\' ExpFolder]);
    try
        load VRDLPParameter.mat;
        cd(curpwd);
        RepeatNumber = max([SP.Repeat RepeatNumber]);
    catch
        load(['exp' num2str(ExpNumber(i), '%03d') '.mat']);
        cd(curpwd);
        RepeatNumber = max([StimulusParameter.Repeat RepeatNumber]);
    end
end
set(handles.Repeat_list, 'String', num2str([1:RepeatNumber]'));

function Spot_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Voltage_text_Callback(hObject, eventdata, handles)

function Voltage_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Sensitivity_text_Callback(hObject, eventdata, handles)

function Sensitivity_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Integration_text_Callback(hObject, eventdata, handles)

function Integration_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SDelay_text_Callback(hObject, eventdata, handles)

function SDelay_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Signal_Length_text_Callback(hObject, eventdata, handles)

function Signal_Length_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Repeat_Exclude_text_Callback(hObject, eventdata, handles)

function Repeat_Exclude_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Go_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
voltage = str2double(get(handles.Voltage_text, 'String'));
sensitivity = str2double(get(handles.Sensitivity_text, 'String'));
Intperiod = str2double(get(handles.Integration_text, 'String'));
signal_length = str2double(get(handles.Signal_Length_text, 'String'));
if isnan(signal_length)
    signal_length = SP.Duration*1000;
end
repeatexc = get(handles.Repeat_Exclude_text, 'String');
if ~isempty(repeatexc)
    repeatexc = strs2numbs(repeatexc);
end
LastTrial = FindLastTrial([DataDir '\exp' num2str(ExpNumber, '%03d')]);
TrialNumber = str2double(LastTrial(6:8));
lsdratio_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
rsdratio_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
lamplitude_mean_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
ramplitude_mean_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
lamplitude_max_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
ramplitude_max_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
langlee_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
langleh_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
ranglee_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
rangleh_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
ldelay_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
rdelay_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width, SP.Repeat-numel(repeatexc));
repeats_all = zeros(SP.ROI_Height/SP.Unit_Height, SP.ROI_Width/SP.Unit_Width);
no_sdelay = [];
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:TrialNumber
    if any(repeatexc == ceil(i/max(SP.Order)))
        WaitSecs(0.01);
        workbar(i/TrialNumber, ['Trial' num2str(i) '/' num2str(TrialNumber)], 'Progress');
        continue;
    end
    lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(i, '%03d') '.mat'];
    sdelay = str2double(get(handles.SDelay_text, 'String'));
    if isnan(sdelay)
        try
            ldata = load(lightpath);
            if isfield(SP, 'PulseWidth')
                [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
            else
                [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
            end
        catch
            sdelay = 0;
        end
        if sdelay == 0
            no_sdelay = [no_sdelay 'Spot' num2str(SP.Order(i)) 'Repeat' num2str(ceil(i/max(SP.Order))) ';'];
        end
    end
    
    index = SP.Order(i);
    row = mod(index, SP.ROI_Height/SP.Unit_Height);
    if row == 0
        row = SP.ROI_Height/SP.Unit_Height;
    end
    column = ceil(index/SP.ROI_Height*SP.Unit_Height);
    if mod(index, SP.ROI_Height/SP.Unit_Height) == 0
        if column ~= index/SP.ROI_Height*SP.Unit_Height
            column = round(index/SP.ROI_Height*SP.Unit_Height);
        end
    end
    
    try
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(i, '%03d') '.mat']);
    catch
        lsdratio_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        ldelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        lamplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        lamplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        langleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        langlee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        rsdratio_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        rdelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        ramplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        ramplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        rangleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        ranglee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column) = repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1;
        workbar(i/TrialNumber, ['Trial' num2str(i) '/' num2str(TrialNumber)], 'Progress');
        continue;
    end
        
    npre = round((APre+sdelay/1000)*SampleRate);
    nlightend = round((APre+sdelay/1000+signal_length/1000)*SampleRate);
    lg = (y(:, 1:3)-voltage/2)/sensitivity;
    lg = -lg;
    lamplitude = sqrt(sum(lg.^2, 2));
    lsd0 = std(lamplitude(1:npre));
    lsd = std(lamplitude(npre+1:nlightend));
    lsdratio = (lsd-lsd0)/lsd0;
    if lsdratio < 0
        lsdratio = 0;
    end
    lsdratio_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = lsdratio;

    lthreshold = 5*lsd0;
    lamean = mean(lamplitude(1:npre));
    temp = find(lamplitude > lamean+lthreshold);
    temp = temp(temp > npre);
    if ~isempty(temp)
        startpointp = temp(1);
    else
        startpointp = [];
    end
    temp = find(lamplitude < lamean-lthreshold);
    temp = temp(temp > npre);
    if ~isempty(temp)
        startpointn = temp(1);
    else
        startpointn = [];
    end
    if ~isempty(startpointp) && ~isempty(startpointn)
        startpoint = min(startpointp, startpointn);
    elseif ~isempty(startpointp) && isempty(startpointn)
        startpoint = startpointp;
    elseif isempty(startpointp) && ~isempty(startpointn)
        startpoint = startpointn;
    else
        startpoint = [];
    end
    if isempty(startpoint)
        ldelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
    else
        ldelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = (x(startpoint)-APre-sdelay/1000)*1000;
    end
    
    if ~isempty(startpoint) && x(startpoint) < ADuration-APost+sdelay/1000
        lxphi = atan2d(lg(:, 2), lg(:, 3));
        lytheta = atand(-lg(:, 1)./sqrt(lg(:, 2).^2+lg(:, 3).^2));
        lxphi0 = mean(lxphi(1:npre));
        lytheta0 = mean(lytheta(1:npre));
        lgstandard = zeros(size(lg));
        for j = 1:length(x)
            lgstandard(j, :) = ([cosd(lytheta0) 0 -sind(lytheta0); 0 1 0; sind(lytheta0) 0 cosd(lytheta0)]^(-1)*[1 0 0; 0 cosd(lxphi0) sind(lxphi0); 0 -sind(lxphi0) cosd(lxphi0)]^(-1)*[lg(j, 1); lg(j, 2); lg(j, 3)])';
        end
        lg0 = mean(lgstandard(1:npre, 3));
        lgstandard(:, 3) = lgstandard(:, 3)-lg0;
        lgx = lgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 1);
        lgy = lgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 2);
        lgz = lgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 3);
        lamplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = mean(sqrt(lgx.^2+lgy.^2+lgz.^2));
        lamplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = max(sqrt(lgx.^2+lgy.^2+lgz.^2));
        
        lahorizontal = atan2d(sum(lgy), sum(lgx));
        laelevation = atand(sum(lgz)/sqrt(sum(lgx)^2+sum(lgy)^2));
        langleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = lahorizontal;
        langlee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = laelevation;
    else
        lamplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = 0;
        lamplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = 0;
        langleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        langlee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
    end
    
    rg = (y(:, 4:6)-voltage/2)/sensitivity;
    rg = -rg;
    ramplitude = sqrt(sum(rg.^2, 2));
    rsd0 = std(ramplitude(1:npre));
    rsd = std(ramplitude(npre+1:nlightend));
    rsdratio = (rsd-rsd0)/rsd0;
    if rsdratio < 0
        rsdratio = 0;
    end
    rsdratio_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = rsdratio;

    rthreshold = 5*rsd0;
    ramean = mean(ramplitude(1:npre));
    temp = find(ramplitude > ramean+rthreshold);
    temp = temp(temp > npre);
    if ~isempty(temp)
        startpointp = temp(1);
    else
        startpointp = [];
    end
    temp = find(ramplitude < ramean-rthreshold);
    temp = temp(temp > npre);
    if ~isempty(temp)
        startpointn = temp(1);
    else
        startpointn = [];
    end
    if ~isempty(startpointp) && ~isempty(startpointn)
        startpoint = min(startpointp, startpointn);
    elseif ~isempty(startpointp) && isempty(startpointn)
        startpoint = startpointp;
    elseif isempty(startpointp) && ~isempty(startpointn)
        startpoint = startpointn;
    else
        startpoint = [];
    end
    if isempty(startpoint)
        rdelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
    else
        rdelay_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = (x(startpoint)-APre-sdelay/1000)*1000;
    end
    
    if ~isempty(startpoint) && x(startpoint) < ADuration-APost+sdelay/1000
        rxphi = atan2d(rg(:, 2), rg(:, 3));
        rytheta = atand(-rg(:, 1)./sqrt(rg(:, 2).^2+rg(:, 3).^2));
        rxphi0 = mean(rxphi(1:npre));
        rytheta0 = mean(rytheta(1:npre));
        rgstandard = zeros(size(rg));
        for j = 1:length(x)
            rgstandard(j, :) = ([cosd(rytheta0) 0 -sind(rytheta0); 0 1 0; sind(rytheta0) 0 cosd(rytheta0)]^(-1)*[1 0 0; 0 cosd(rxphi0) sind(rxphi0); 0 -sind(rxphi0) cosd(rxphi0)]^(-1)*[rg(j, 1); rg(j, 2); rg(j, 3)])';
        end
        rg0 = mean(rgstandard(1:npre, 3));
        rgstandard(:, 3) = rgstandard(:, 3)-rg0;
        rgx = rgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 1);
        rgy = rgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 2);
        rgz = rgstandard(x>=x(startpoint)&x<=x(startpoint)+Intperiod/1000, 3);
        ramplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = mean(sqrt(rgx.^2+rgy.^2+rgz.^2));
        ramplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = max(sqrt(rgx.^2+rgy.^2+rgz.^2));
        
        rahorizontal = atan2d(sum(rgy), sum(rgx));
        raelevation = atand(sum(rgz)/sqrt(sum(rgx)^2+sum(rgy)^2));
        rangleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = rahorizontal;
        ranglee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = raelevation;
    else
        ramplitude_mean_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = 0;
        ramplitude_max_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = 0;
        rangleh_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
        ranglee_all(SP.ROI_Height/SP.Unit_Height-row+1, column, repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1) = NaN;
    end
    
    repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column) = repeats_all(SP.ROI_Height/SP.Unit_Height-row+1, column)+1;
    workbar(i/TrialNumber, ['Trial' num2str(i) '/' num2str(TrialNumber)], 'Progress');
end
if ~isempty(no_sdelay)
    warndlg(['For ' no_sdelay ' the light Delay cannot be calculated, so is set to 0.'], 'Warning');
end
UserData.lsdratio_all = lsdratio_all;
UserData.rsdratio_all = rsdratio_all;
UserData.lamplitude_mean_all = lamplitude_mean_all;
UserData.ramplitude_mean_all = ramplitude_mean_all;
UserData.lamplitude_max_all = lamplitude_max_all;
UserData.ramplitude_max_all = ramplitude_max_all;
UserData.langleh_all = langleh_all;
UserData.langlee_all = langlee_all;
UserData.rangleh_all = rangleh_all;
UserData.ranglee_all = ranglee_all;
UserData.ldelay_all = ldelay_all;
UserData.rdelay_all = rdelay_all;
UserData.ExpNumber = ExpNumber;
UserData.SP = SP;
set(handles.Go_button, 'UserData', UserData);

function myimshow(hObject, eventdata, him)
figure('Color', [1 1 1]);
axes;
set(gca, 'YDir', 'reverse');
copyobj(him, gca);
axis image;
axis off;
UserData = get(him, 'UserData');
switch UserData{1}
    case 1
        colormap jet;
        colorbar;
    case 2
        colormap hsv;
        try
            colorbar('Ticks', [0 0.5 1], 'TickLabels', {'0', '180', '360'});
        catch
            colorbar('YTick', [0.5 32.5 64.5], 'YTickLabel', {'0', '180', '360'});
        end
    case 3
        colormap hsv;
        try
            colorbar('Limits', [0 0.5], 'Ticks', [0 0.25 0.5], 'TickLabels', {'-90', '0', '90'});
        catch
            colorbar('YLimit', [0 32], 'YTick', [0.5 16.5 32.5], 'YTickLabel', {'-90', '0', '90'});
        end
    case 4
        colormap gray;
        hcb = colorbar;
        try
            set(hcb, 'Ticks', [0 0.5 1], 'TickLabels', {num2str(max(UserData{3}(:))) num2str(max(UserData{3}(:))/2+min(UserData{3}(:))/2) num2str(min(UserData{3}(:)))});
        catch
            set(hcb, 'YTick', [0.5 32.5 64.5], 'YTickLabel', {num2str(max(UserData{3}(:))) num2str(max(UserData{3}(:))/2+min(UserData{3}(:))/2) num2str(min(UserData{3}(:)))});
        end
    case 5
        clf;
        copyobj(UserData{3}, gcf);
        set(gca, 'Position', [0.13 0.11 0.775 0.815]);
        colormap jet;
        hcb = colorbar;
        title(hcb, UserData{4});
end
title(UserData{2});

function WDuration_text_Callback(hObject, eventdata, handles)

function WDuration_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DInterval_text_Callback(hObject, eventdata, handles)

function DInterval_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FAnalysis_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
handles.parameter_location = [DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat'];
load(handles.parameter_location);
wduration = str2double(get(handles.WDuration_text, 'String'));
wduration = round(wduration*SampleRate);
dinterval = str2double(get(handles.DInterval_text, 'String'));
dinterval = round(dinterval/1000*SampleRate);
winfun = hanning(wduration);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(wduration);          % Next highest power of 2 greater than length(x).
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    try
        lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat'];
        sdelay = str2double(get(handles.SDelay_text, 'String'));
        if isnan(sdelay)
            try
                ldata = load(lightpath);
                if isfield(SP, 'PulseWidth')
                    [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
                else
                    [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
                end
            catch
                sdelay = 0;
            end
            if sdelay == 0
                warndlg('Light Delay cannot be calculated, so is set to 0.', ['Repeat' num2str(RepeatNumber(i))]);
            end
        end
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mat']);
        if Accelerometer_check
            seeg = y(:, 7)/Gain;
        else
            seeg = y/Gain;
        end
        if get(handles.EEG_Reverse_check, 'Value')
            seeg = -seeg;
        end
        figure('Color', [1 1 1]);
        subplot(2, 3, 1);
        hp = plot(x, seeg, '-k');
        hold on;
        yrange = max(seeg)-min(seeg);
        ylim = [min(seeg)-0.05*yrange max(seeg)+0.05*yrange];
        plot([APre+sdelay/1000 APre+sdelay/1000], ylim, ':k');
        plot([ADuration-APost+sdelay/1000 ADuration-APost+sdelay/1000], ylim, ':k');
        set(gca, 'YLim', ylim, 'XLim', [min(x) max(x)]);
        box off;
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend(hp, 'EEG');
        set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
        
        signal = seeg.*hanning(length(seeg));
        NFFT0 = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
        FFTX = fft(signal, NFFT0);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
        NumUniquePts = ceil((NFFT0+1)/2);
        FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
        MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
        MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
        MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
        MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
        MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
        MX = MX.^2;
        f = (0:NumUniquePts-1)*2*Fn/NFFT0;
        percentp = MX/sum(MX)*100;
        subplot(2, 3, 2);
        hp = plot(f(f <= 130), percentp(f <= 130), '-k');
        yrange = max(percentp(f <= 130))-min(percentp(f <= 130));
        ylim = [min(percentp(f <= 130))-0.05*yrange max(percentp(f <= 130))+0.05*yrange];
        set(gca, 'YLim', ylim, 'XLim', [0 130]);
        box off;
        xlabel('Frequency (Hz)');
        ylabel('Percentage of Power (%)');
        legend(hp, 'Percentage of Power');
        set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
        
        delta =(f < 4);
        theta = (f >= 4 & f < 8);
        alpha = (f >= 8 & f < 14);
        beta = (f >= 14 & f < 30);
        gama = (f >= 30 & f < 130);
        percentb(1) = sum(MX(delta))/sum(MX);
        percentb(2) = sum(MX(theta))/sum(MX);
        percentb(3) = sum(MX(alpha))/sum(MX);
        percentb(4) = sum(MX(beta))/sum(MX);
        percentb(5) = sum(MX(gama))/sum(MX);
        subplot(2, 3, 3);
        bar(percentb, 'FaceColor', 'black', 'EdgeColor', 'white', 'BarWidth', 0.8);
        box off;
        set(gca, 'XTick', [1:5], 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gama'}, 'TickDir', 'out');
        ylabel('Percentage (%)');
        set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
        
        fmap = [];
        n = 0;
        for j = 1:dinterval:length(seeg)-wduration+1
            signal = seeg(j:j+wduration-1);
            signal = signal.*winfun;
            FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
            NumUniquePts = ceil((NFFT+1)/2);
            FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
            MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
            MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
            MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
            MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
            MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
            MX = MX.^2;
            f = (0:NumUniquePts-1)*2*Fn/NFFT;
            fmap(:, n+1) = MX(f <= 130);
            n = n+1;
        end
        xtime = x(1+floor(wduration/2):dinterval:length(seeg)-wduration+1+floor(wduration/2))';
        yfrequency = f(f <= 130)';
        subplot(2, 3, 4);
        him = imshow(flipud(10*log10(fmap)), []);
        colormap jet;
        hcb = colorbar;
        title(hcb, 'Power (dB)');
        axis on;
        set(gca, 'XTick', [1:floor(numel(xtime)/6):numel(xtime)], 'XTickLabel', num2str(xtime([1:floor(numel(xtime)/6):numel(xtime)])));
        set(gca, 'YTick', fliplr(numel(yfrequency)+1-[1:floor(numel(yfrequency)/6):numel(yfrequency)]), 'YTickLabel', num2str(flipud(yfrequency([1:floor(numel(yfrequency)/6):numel(yfrequency)])), '%.1f'));
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Power Spectrum');
        set(him, 'UserData', {5, 'Power Spectrum', gca, 'Power (dB)'});
        set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
        
        subplot(2, 3, 5);
        fmean = mean(fmap, 2);
        fmean = fmean/sum(fmean)*100;
        hp = plot(f(f <= 130), fmean, '-k');
        yrange = max(fmean)-min(fmean);
        ylim = [min(fmean)-0.05*yrange max(fmean)+0.05*yrange];
        set(gca, 'YLim', ylim, 'XLim', [0 130]);
        box off;
        xlabel('Frequency (Hz)');
        ylabel('Percentage of Power (%)');
        legend(hp, 'Percentage of Power');
        set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
        
        percentb(1) = sum(fmean(f < 4))/sum(fmean);
        percentb(2) = sum(fmean(f >= 4 & f < 8))/sum(fmean);
        percentb(3) = sum(fmean(f >= 8 & f < 14))/sum(fmean);
        percentb(4) = sum(fmean(f >= 14 & f < 30))/sum(fmean);
        percentb(5) = sum(fmean(f >= 30 & f < 130))/sum(fmean);
        subplot(2, 3, 6);
        bar(percentb, 'FaceColor', 'black', 'EdgeColor', 'white', 'BarWidth', 0.8);
        box off;
        set(gca, 'XTick', [1:5], 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gama'}, 'TickDir', 'out');
        ylabel('Percentage (%)');
        set(gca, 'ButtonDownFcn', @(hObject, eventdata)myplot(hObject, eventdata, gca));
    catch
        errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function CAverage_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
handles.parameter_location = [DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat'];
load(handles.parameter_location);
seeg_all = [];
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    try
        lightpath = [DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat'];
        sdelay = str2double(get(handles.SDelay_text, 'String'));
        if isnan(sdelay)
            try
                ldata = load(lightpath);
                if isfield(SP, 'PulseWidth')
                    [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
                else
                    [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
                end
            catch
                sdelay = 0;
            end
            if sdelay == 0
                warndlg('Light Delay cannot be calculated, so is set to 0.', ['Repeat' num2str(RepeatNumber(i))]);
            end
        end
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\trial' num2str(TrialNumber, '%03d') '.mat']);
        if Accelerometer_check
            seeg = y(:, 7)/Gain;
        else
            seeg = y/Gain;
        end
        if get(handles.EEG_Reverse_check, 'Value')
            seeg = -seeg;
        end
        seeg(x < sdelay/1000) = [];
        seeg = [seeg; zeros(numel(x)-numel(seeg), 1)/0];
        seeg_all = [seeg_all seeg];
    catch
        errordlg(['Cannot add exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end
if numel(RepeatNumber) == 1
    figure;
    hp = plot(x, seeg_all, '-k');
    hold on;
    yrange = max(seeg_all)-min(seeg_all);
    ylim = [min(seeg_all)-0.05*yrange max(seeg_all)+0.05*yrange];
    plot([APre APre], ylim, ':k');
    plot([ADuration-APost ADuration-APost], ylim, ':k');
    set(gca, 'YLim', ylim);
    box off;
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend(hp, 'Average EEG');
else
    seeg_mean = mean(seeg_all, 2);
    seeg_sem = std(seeg_all, 0, 2)/sqrt(size(seeg_all, 2));
    temp = isnan(seeg_mean);
    seeg_mean(temp) = [];
    seeg_sem(temp) = [];
    x(temp) = [];
    figure;
    patch([x; x(end:-1:1)], [seeg_mean+seeg_sem; seeg_mean(end:-1:1)-seeg_sem(end:-1:1)], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    hold on;
    hp = plot(x, seeg_mean, '-k');
    yrange = max(seeg_mean+seeg_sem)-min(seeg_mean-seeg_sem);
    ylim = [min(seeg_mean-seeg_sem)-0.05*yrange max(seeg_mean+seeg_sem)+0.05*yrange];
    plot([APre APre], ylim, ':k');
    plot([ADuration-APost ADuration-APost], ylim, ':k');
    set(gca, 'YLim', ylim);
    box off;
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend(hp, 'Average EEG');
end

function RefreshLists_Callback(hObject, eventdata, handles)
RefreshLists(hObject, eventdata, handles);

function Save_check_Callback(hObject, eventdata, handles)

function MapPlot_Button_Callback(hObject, eventdata, handles)
UserData = get(handles.Go_button, 'UserData');
lsdratio_all = UserData.lsdratio_all;
rsdratio_all = UserData.rsdratio_all;
lamplitude_mean_all = UserData.lamplitude_mean_all;
ramplitude_mean_all = UserData.ramplitude_mean_all;
lamplitude_max_all = UserData.lamplitude_max_all;
ramplitude_max_all = UserData.ramplitude_max_all;
langleh_all = UserData.langleh_all;
langlee_all = UserData.langlee_all;
rangleh_all = UserData.rangleh_all;
ranglee_all = UserData.ranglee_all;
ldelay_all = UserData.ldelay_all;
rdelay_all = UserData.rdelay_all;
ExpNumber = UserData.ExpNumber;
SP = UserData.SP;

lsdratio_mean = compute_mean(lsdratio_all);
rsdratio_mean = compute_mean(rsdratio_all);
hf1 = figure('NumberTitle', 'off', 'Name', ['Exp' num2str(ExpNumber, '%03d')], 'Color', [1 1 1]);
subplot(2, 3, 1, 'Parent', hf1);
him = imshow(imresize(lsdratio_mean, [size(lsdratio_mean, 1)*SP.Unit_Height size(lsdratio_mean, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Left Mean SD Ratio Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Left Mean SD Ratio Map');
subplot(2, 3, 4, 'Parent', hf1);
him = imshow(imresize(rsdratio_mean, [size(rsdratio_mean, 1)*SP.Unit_Height size(rsdratio_mean, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Right Mean SD Ratio Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Right Mean SD Ratio Map');

lamplitude_mean = compute_mean(lamplitude_mean_all);
ramplitude_mean = compute_mean(ramplitude_mean_all);
subplot(2, 3, 2, 'Parent', hf1);
him = imshow(imresize(lamplitude_mean, [size(lamplitude_mean, 1)*SP.Unit_Height size(lamplitude_mean, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Left Mean Acceleration Amplitude Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Left Mean Acceleration Amplitude Map');
subplot(2, 3, 5, 'Parent', hf1);
him = imshow(imresize(ramplitude_mean, [size(ramplitude_mean, 1)*SP.Unit_Height size(ramplitude_mean, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Right Mean Acceleration Amplitude Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Right Mean Acceleration Amplitude Map');

lamplitude_max = compute_mean(lamplitude_max_all);
ramplitude_max = compute_mean(ramplitude_max_all);
subplot(2, 3, 3, 'Parent', hf1);
him = imshow(imresize(lamplitude_max, [size(lamplitude_max, 1)*SP.Unit_Height size(lamplitude_max, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Left Max Acceleration Amplitude Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Left Max Acceleration Amplitude Map');
subplot(2, 3, 6, 'Parent', hf1);
him = imshow(imresize(ramplitude_max, [size(ramplitude_max, 1)*SP.Unit_Height size(ramplitude_max, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {1, 'Right Max Acceleration Amplitude Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap jet;
colorbar;
title('Right Max Acceleration Amplitude Map');

langleh_mean = zeros(size(langleh_all, 1), size(langleh_all, 2));
langlee_mean = zeros(size(langlee_all, 1), size(langlee_all, 2));
rangleh_mean = zeros(size(rangleh_all, 1), size(rangleh_all, 2));
ranglee_mean = zeros(size(ranglee_all, 1), size(ranglee_all, 2));
for i = 1:size(langleh_mean, 1)
    for j = 1:size(langleh_mean, 2)
        temp = squeeze(langleh_all(i, j, :));
        temp(temp < 0) = temp(temp < 0)+360;
        temp(isnan(temp)) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        langleh_mean(i, j) = temp;
        
        temp = squeeze(langlee_all(i, j, :));
        temp(isnan(temp)) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        langlee_mean(i, j) = temp;
        
        temp = squeeze(rangleh_all(i, j, :));
        temp(temp < 0) = temp(temp < 0)+360;
        temp(isnan(temp)) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        rangleh_mean(i, j) = temp;
        
        temp = squeeze(ranglee_all(i, j, :));
        temp(isnan(temp)) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        ranglee_mean(i, j) = temp;
    end
end
lhsvmaph(:, :, 1) = langleh_mean/360;
lhsvmaph(:, :, 2) = 1;
lhsvmaph(:, :, 3) = mat2gray(lamplitude_mean);
hf2 = figure('NumberTitle', 'off', 'Name', ['Exp' num2str(ExpNumber, '%03d')], 'Color', [1 1 1]);
subplot(2, 2, 1, 'Parent', hf2);
him = imshow(imresize(hsv2rgb(lhsvmaph), [size(lhsvmaph, 1)*SP.Unit_Height size(lhsvmaph, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {2, 'Left Horizontal Angle Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap hsv;
try
    colorbar('Ticks', [0 0.5 1], 'TickLabels', {'0', '180', '360'});
catch
    colorbar('YTick', [0.5 32.5 64.5], 'YTickLabel', {'0', '180', '360'});
end
title('Left Horizontal Angle Map');

lhsvmape(:, :, 1) = (langlee_mean+90)/360;
lhsvmape(:, :, 2) = 1;
lhsvmape(:, :, 3) = mat2gray(lamplitude_mean);
subplot(2, 2, 2, 'Parent', hf2);
him = imshow(imresize(hsv2rgb(lhsvmape), [size(lhsvmape, 1)*SP.Unit_Height size(lhsvmape, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {3, 'Left Elevation Angle Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap hsv;
try
    colorbar('Limits', [0 0.5], 'Ticks', [0 0.25 0.5], 'TickLabels', {'-90', '0', '90'});
catch
    colorbar('YLimit', [0 32], 'YTick', [0.5 16.5 32.5], 'YTickLabel', {'-90', '0', '90'});
end
title('Left Elevation Angle Map');

rhsvmaph(:, :, 1) = rangleh_mean/360;
rhsvmaph(:, :, 2) = 1;
rhsvmaph(:, :, 3) = mat2gray(ramplitude_mean);
subplot(2, 2, 3, 'Parent', hf2);
him = imshow(imresize(hsv2rgb(rhsvmaph), [size(rhsvmaph, 1)*SP.Unit_Height size(rhsvmaph, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {2, 'Right Horizontal Angle Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap hsv;
try
    colorbar('Ticks', [0 0.5 1], 'TickLabels', {'0', '180', '360'});
catch
    colorbar('YTick', [0.5 32.5 64.5], 'YTickLabel', {'0', '180', '360'});
end
title('Right Horizontal Angle Map');

rhsvmape(:, :, 1) = (ranglee_mean+90)/360;
rhsvmape(:, :, 2) = 1;
rhsvmape(:, :, 3) = mat2gray(ramplitude_mean);
subplot(2, 2, 4, 'Parent', hf2);
him = imshow(imresize(hsv2rgb(rhsvmape), [size(rhsvmape, 1)*SP.Unit_Height size(rhsvmape, 2)*SP.Unit_Width], 'nearest'), []);
set(him, 'UserData', {3, 'Right Elevation Angle Map'});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap hsv;
try
    colorbar('Limits', [0 0.5], 'Ticks', [0 0.25 0.5], 'TickLabels', {'-90', '0', '90'});
catch
    colorbar('YLimit', [0 32], 'YTick', [0.5 16.5 32.5], 'YTickLabel', {'-90', '0', '90'});
end
title('Right Elevation Angle Map');

ldelay_mean = zeros(size(ldelay_all, 1), size(ldelay_all, 2));
rdelay_mean = zeros(size(rdelay_all, 1), size(rdelay_all, 2));
for i = 1:size(ldelay_mean, 1)
    for j = 1:size(ldelay_mean, 2)
        temp = squeeze(ldelay_all(i, j, :));
        temp(isnan(temp)) = [];
        temp(temp/1000 >= SP.Duration) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        ldelay_mean(i, j) = temp;
        
        temp = squeeze(rdelay_all(i, j, :));
        temp(isnan(temp)) = [];
        temp(temp/1000 >= SP.Duration) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        rdelay_mean(i, j) = temp;
    end
end
temp = 1-mat2gray(ldelay_mean);
r = temp;
r(isnan(ldelay_mean)) = 1;
r = imresize(r, [size(r, 1)*SP.Unit_Height size(r, 2)*SP.Unit_Width], 'nearest');
g = temp;
g(isnan(ldelay_mean)) = 0;
g = imresize(g, [size(g, 1)*SP.Unit_Height size(g, 2)*SP.Unit_Width], 'nearest');
b = temp;
b(isnan(ldelay_mean)) = 0;
b = imresize(b, [size(b, 1)*SP.Unit_Height size(b, 2)*SP.Unit_Width], 'nearest');
ldelay_map(:, :, 1) = r;
ldelay_map(:, :, 2) = g;
ldelay_map(:, :, 3) = b;
hf3 = figure('NumberTitle', 'off', 'Name', ['Exp' num2str(ExpNumber, '%03d')], 'Color', [1 1 1]);
haxes = subplot(1, 2, 1, 'Parent', hf3);
axes(haxes);
him = imshow(ldelay_map);
set(him, 'UserData', {4, 'Left Movement Delay Map', ldelay_mean});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap gray;
hcb = colorbar;
try
    set(hcb, 'Ticks', [0 0.5 1], 'TickLabels', {num2str(max(ldelay_mean(:))) num2str(max(ldelay_mean(:))/2+min(ldelay_mean(:))/2) num2str(min(ldelay_mean(:)))});
catch
    set(hcb, 'YTick', [0.5 32.5 64.5], 'YTickLabel', {num2str(max(ldelay_mean(:))) num2str(max(ldelay_mean(:))/2+min(ldelay_mean(:))/2) num2str(min(ldelay_mean(:)))});
end
title('Left Movement Delay Map');

temp = 1-mat2gray(rdelay_mean);
r = temp;
r(isnan(rdelay_mean)) = 1;
r = imresize(r, [size(r, 1)*SP.Unit_Height size(r, 2)*SP.Unit_Width], 'nearest');
g = temp;
g(isnan(rdelay_mean)) = 0;
g = imresize(g, [size(g, 1)*SP.Unit_Height size(g, 2)*SP.Unit_Width], 'nearest');
b = temp;
b(isnan(rdelay_mean)) = 0;
b = imresize(b, [size(b, 1)*SP.Unit_Height size(b, 2)*SP.Unit_Width], 'nearest');
rdelay_map(:, :, 1) = r;
rdelay_map(:, :, 2) = g;
rdelay_map(:, :, 3) = b;
haxes = subplot(1, 2, 2, 'Parent', hf3);
axes(haxes);
him = imshow(rdelay_map);
set(him, 'UserData', {4, 'Right Movement Delay Map', rdelay_mean});
set(him, 'ButtonDownFcn', @(hObject, eventdata)myimshow(hObject, eventdata, him));
colormap gray;
hcb = colorbar;
try
    set(hcb, 'Ticks', [0 0.5 1], 'TickLabels', {num2str(max(rdelay_mean(:))) num2str(max(rdelay_mean(:))/2+min(rdelay_mean(:))/2) num2str(min(rdelay_mean(:)))});
catch
    set(hcb, 'YTick', [0.5 32.5 64.5], 'YTickLabel', {num2str(max(rdelay_mean(:))) num2str(max(rdelay_mean(:))/2+min(rdelay_mean(:))/2) num2str(min(rdelay_mean(:)))});
end
title('Right Movement Delay Map');

function data_mean = compute_mean(data)
data_mean = zeros(size(data, 1), size(data, 2));
for i = 1:size(data_mean, 1)
    for j = 1:size(data_mean, 2)
        temp = squeeze(data(i, j, :));
        temp(isnan(temp)) = [];
        if isempty(temp)
            temp = NaN;
        else
            temp = mean(temp);
        end
        data_mean(i, j) = temp;
    end
end

function LP_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
handles.parameter_location = [DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat'];
load(handles.parameter_location);
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    try
        load([DataDir '\exp' num2str(ExpNumber, '%03d') '\LData' num2str(TrialNumber, '%03d') '.mat']);
        if isfield(SP, 'PulseWidth')
            y = y(:, 2);
        else
            y = y(:, 1);
        end
        figure;
        hp = plot(x, y, '-k');
        hold on;
        yrange = max(y)-min(y);
        ylim = [min(y)-0.05*yrange max(y)+0.05*yrange];
        plot([APre APre], ylim, ':k');
        plot([ADuration-APost ADuration-APost], ylim, ':k');
        set(gca, 'YLim', ylim, 'XLim', [min(x) max(x)]);
        box off;
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend(hp, 'Light Profile');
        
        [sdelay, threshold] = Light_Delay(x, y, APre, SampleRate, SP);
        plot([x(1) x(end)], [threshold threshold], ':k');
        if sdelay ~= 0
            title(['Delay of the light is ' num2str(sdelay) ' ms']);
        end
    catch
        errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': trial' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function CutVideo_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if isempty(ExpNumber)
    errordlg('You forgot to choose experiment !', 'ERROR');
    return;
end
% threshold_const = 5;
for i = 1:length(ExpNumber)
    handles.parameter_location = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat'];
    load(handles.parameter_location);
    try
        videodir = [DataDir '\LG Webcam Videos\Video ' num2str(ExpNumber(i)) '.mp4'];
        readerobj = VideoReader(videodir);
    catch
        videodir = [DataDir '\LG Webcam Videos\Video ' num2str(ExpNumber(i)) '.wmv'];
        readerobj = VideoReader(videodir);
    end
    FrameRate = readerobj.FrameRate;
    FrameRate = round(FrameRate);
    if FrameRate == 15
        FrameRate = 30;
        msgbox(['FrameRate of Exp ' num2str(ExpNumber(i)) ' is changed from 15 to 30']);
    end
    Frames = readerobj.NumberOfFrames;
    readerobj = VideoReader(videodir);
    
    im = readFrame(readerobj);
    imhsv = rgb2hsv(im);
    n = 1;
    figure; imshow(im, []);
    [~, rect] = imcrop;
    ROI = imcrop(imhsv(:, :, 2), rect);
    close gcf;
    
    Frames_pre = APre*FrameRate;
    images(:, :, :, 1) = im;
    hsvvalue(1) = mean(ROI(:));
    frameid = 1;
    figure;
    plot(frameid, 0, 'k.');
    hold on;
    hp = gca;
    
    lightid = 0;
    threshold = [];
    lightonframe = [];
    
    while hasFrame(readerobj)
        frameid = frameid+1;
        n = n+1;
        im = readFrame(readerobj);
        imhsv = rgb2hsv(im);
        ROI = imcrop(imhsv(:, :, 2), rect);
        images(:, :, :, n) = im;
        hsvvalue(n) = mean(ROI(:));
        try
            plot(hp, frameid, hsvvalue(n)-hsvvalue(n-1), 'k.');
        end
%         drawnow;
        if lightid == 0
            if isempty(threshold) && n == 156
                threshold = 0.01;
                plot(hp, frameid, threshold, 'r.');
%                 drawnow;
            end
        else
            if isempty(threshold) && n == Frames_pre+1*FrameRate
                threshold = max(diff(hsvvalue));
                plot(hp, frameid, threshold, 'r.');
%                 drawnow;
            end
        end
        if ~isempty(threshold)
            if hsvvalue(n)-hsvvalue(n-1) >= threshold
                lightid = lightid+1;
                lightonframe(lightid) = n;
                threshold = [];
            end
        end
        if lightid ~= 0 && isempty(threshold)
            if n == lightonframe(lightid)+round((ADuration-APre)*FrameRate)-1
                tic
                imvideo = images;
                cutvideodir = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\webcam' num2str(lightid, '%03d') '.avi'];
                vidObj = VideoWriter(cutvideodir, 'Motion JPEG AVI');
                set(vidObj, 'FrameRate', FrameRate, 'Quality', 75);
                open(vidObj);
                writeVideo(vidObj, imvideo);
                close(vidObj);
                clear images satuation;
                n = 0;
                hsvvalue = [];
                toc
            end
        end
    end
    save([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\LGwebcam_lightonframe.mat'], 'lightonframe');
    disp('Done!');
end

function WCVideoPlay_button_Callback(hObject, eventdata, handles)
handles.VideoSource = 'LG Webcam';
PlayVideo(hObject, eventdata, handles);

function WCAnalysis_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
pre = APre;
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    try
        try
            readerobj = VideoReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.mp4']);
        catch
            readerobj = VideoReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.avi']);
        end
        FrameRate = readerobj.FrameRate;
        FrameRate = round(FrameRate);
        Frames = readerobj.NumberOfFrames;
        centroid_all = [];
        
        try
            videoFileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.mp4']);
        catch
            videoFileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber, '%03d') '\webcam' num2str(TrialNumber, '%03d') '.avi']);
        end
        videoFrame = step(videoFileReader);
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
        
        lsindex = pre*FrameRate+1;
        leindex = lightendframe(TrialNumber)-1;
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
        distance = sum(sqrt(diff(xls2le).^2+diff(yls2le).^2));
        title(haxes, ['The total distance travelled during light on is ' num2str(distance) ' pixels']);
        
        coord1 = centroid_all;
        coord2 = circshift(centroid_all, -1);
        x = coord2(:, 1) - coord1(:, 1);
        y = coord2(:, 2) - coord1(:, 2);
        hftrajectory = figure;
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
        
        if get(handles.Save_check, 'Value')
            saveas(hf, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Moving_Points'], 'fig');
            saveas(hftrajectory, ['Spot' num2str(SpotNumber, '%03d') '_Repeat' num2str(RepeatNumber(i), '%03d') '_Moving_Trajectory'], 'fig');
        end
        
        % Clean up
        release(videoFileReader);
        release(videoPlayer);
        release(pointTracker);
    catch
        errordlg(['Cannot analyze exp' num2str(ExpNumber, '%03d') ': webcam' num2str(TrialNumber, '%03d')], 'ERROR');
    end
end

function RenameFiles_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
curpwd = pwd;
cd([DataDir '\LG Webcam Videos']);
DirList = dir;
n= 0;
id = [];
for i=3:size(DirList)
    if strncmp(DirList(i).name, 'Video ', 6)
        n = n+1;
        videonames{n} = DirList(i).name;
        temp = DirList(i).name;
        temp(1:6) = [];
        temp(end-3:end) = [];
        id = [id str2double(temp)];
    end
end
cd(curpwd);
[~, index] = sort(id);
for i = 1:n
    movefile([DataDir '\LG Webcam Videos\' videonames{index(i)}], [DataDir '\LG Webcam Videos\Video ' num2str(i) '.mp4']);
end

function CAverageExp_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
seeg_all = [];
for i = 1:length(ExpNumber)
    handles.parameter_location = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat'];
    load(handles.parameter_location);
    if isempty(seeg_all)
        seeg_all = cell(1, max(SP.Order));
    end
    for j = 1:length(SP.Order)
        try
            lightpath = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\LData' num2str(j, '%03d') '.mat'];
            sdelay = str2double(get(handles.SDelay_text, 'String'));
            if isnan(sdelay)
                try
                    ldata = load(lightpath);
                    if isfield(SP, 'PulseWidth')
                        [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 2), APre, SampleRate, SP);
                    else
                        [sdelay, threshold] = Light_Delay(ldata.x, ldata.y(:, 1), APre, SampleRate, SP);
                    end
                catch
                    sdelay = 0;
                end
                if sdelay == 0
                    warndlg('Light Delay cannot be calculated, so is set to 0.', ['Exp' num2str(ExpNumber(i), '%03d') ': LData' num2str(j, '%03d')]);
                end
            end
            load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\trial' num2str(j, '%03d') '.mat']);
            if Accelerometer_check
                seeg = y(:, 7)/Gain;
            else
                seeg = y/Gain;
            end
            if get(handles.EEG_Reverse_check, 'Value')
                seeg = -seeg;
            end
            seeg(x < sdelay/1000) = [];
            seeg = [seeg; zeros(numel(x)-numel(seeg), 1)/0];
            seeg_all{SP.Order(j)} = [seeg_all{SP.Order(j)} seeg];
        catch
            errordlg(['Cannot add exp' num2str(ExpNumber(i), '%03d') ': trial' num2str(j, '%03d')], 'ERROR');
        end
    end
end

for i = 1:size(seeg_all, 2)
    seeg = seeg_all{i};
    if size(seeg, 2) == 1
        figure;
        hp = plot(x, seeg, '-k');
        hold on;
        yrange = max(seeg)-min(seeg);
        ylim = [min(seeg)-0.05*yrange max(seeg)+0.05*yrange];
        plot([APre APre], ylim, ':k');
        plot([ADuration-APost ADuration-APost], ylim, ':k');
        set(gca, 'YLim', ylim);
        box off;
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        title(['Spot ' num2str(i)]);
        legend(hp, 'Average EEG');
    else
        seeg_mean = mean(seeg, 2);
        seeg_sem = std(seeg, 0, 2)/sqrt(size(seeg, 2));
        temp = isnan(seeg_mean);
        seeg_mean(temp) = [];
        seeg_sem(temp) = [];
        x1 = x;
        x1(temp) = [];
        figure;
        patch([x1; x1(end:-1:1)], [seeg_mean+seeg_sem; seeg_mean(end:-1:1)-seeg_sem(end:-1:1)], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        hold on;
        hp = plot(x1, seeg_mean, '-k');
        yrange = max(seeg_mean+seeg_sem)-min(seeg_mean-seeg_sem);
        ylim = [min(seeg_mean-seeg_sem)-0.05*yrange max(seeg_mean+seeg_sem)+0.05*yrange];
        plot([APre APre], ylim, ':k');
        plot([ADuration-APost ADuration-APost], ylim, ':k');
        set(gca, 'YLim', ylim);
        box off;
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        title(['Spot ' num2str(i)]);
        legend(hp, 'Average EEG');
    end
end

function EEG_Reverse_check_Callback(hObject, eventdata, handles)

function PlotTrajectory_button_Callback(hObject, eventdata, handles)
clc;
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = get(handles.Spot_list, 'Value');
RepeatNumber = get(handles.Repeat_list, 'Value');
FramesBL = str2double(get(handles.FramesBL_text, 'String'));
PG1Y0 = str2double(get(handles.PG1Y0_text, 'String'));
PG1Z0 = str2double(get(handles.PG1Z0_text, 'String'));
PG2X0 = str2double(get(handles.PG2X0_text, 'String'));
PG2Z0 = str2double(get(handles.PG2Z0_text, 'String'));
if isnan(PG1Y0) || isnan(PG1Z0) || isnan(PG2X0) || isnan(PG2Z0)
    try
        load([DataDir '\ReferencePoint.mat']);
        PG1Y0 = Rf.PG1Y0;
        PG1Z0 = Rf.PG1Z0;
        PG2X0 = Rf.PG2X0;
        PG2Z0 = Rf.PG2Z0;
        set(handles.PG1Y0_text, 'String', num2str(PG1Y0));
        set(handles.PG1Z0_text, 'String', num2str(PG1Z0));
        set(handles.PG2X0_text, 'String', num2str(PG2X0));
        set(handles.PG2Z0_text, 'String', num2str(PG2Z0));
    catch
        errordlg('You need to set the Reference Point first!', 'Error');
        return;
    end
end
FrameRate = 100;
resolution = [512 640];
CP1 = load([DataDir '\Camera Calibration\Camera Calibration (PG1)\cameraParams.mat']);
CP1 = CP1.cameraParams;
CP2 = load([DataDir '\Camera Calibration\Camera Calibration (PG2)\cameraParams.mat']);
CP2 = CP2.cameraParams;
inch2mmcst = 25.4;
PG1REF = undistortPoints([PG1Y0 resolution(1)-PG1Z0], CP1);
PG1REF = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), PG1REF);
PG2REF = undistortPoints([PG2X0 resolution(1)-PG2Z0], CP1);
PG2REF = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), PG2REF);
linecolors = jet(256);
markersize = 8;
Fs = FrameRate;                             % sampling frequency
Fn = Fs/2;                                  % Nyquist frequency
for i = 1:length(SpotNumber)
    x_all = [];
    y_all = [];
    z1_all = [];
    z2_all = [];
    if FramesBL > 1
        x_baseline_all = [];
        y_baseline_all = [];
        z1_baseline_all = [];
        z2_baseline_all = [];
    end
    x_all_pixel = [];
    y_all_pixel = [];
    z1_all_pixel = [];
    z2_all_pixel = [];
    PG1repeat_all = [];
    PG2repeat_all = [];
    RepeatNumberAll = 0;
    for j = 1:length(ExpNumber)
        try
            load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\VRDLPParameter.mat']);
            Frames = ceil(FrameRate*SP.Duration);
        catch
            load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\exp' num2str(ExpNumber(j), '%03d') '.mat']);
            SP.Order = StimulusParameter.Order;
            Frames = ceil(FrameRate*StimulusParameter.Duration);
        end
        for k = 1:length(RepeatNumber)
            TrialNumber = find(SP.Order == SpotNumber(i));
            try
                TrialNumber = TrialNumber(RepeatNumber(k));
            catch
                continue;
            end
            try
                v.data = load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
                if get(handles.RForelimb_radiobutton, 'Value')
                    load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG1camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                else
                    load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG1camera' num2str(TrialNumber) '_trajectory.mat']);
                end
                temp = find(v.data(:, 1) == 1);
                y = trajectory(:, 1);
                y_baseline = y(temp(1)-FramesBL:temp(1)-1);
                y = y([temp(1)-1; temp]);
                z1 = trajectory(:, 2);
                z1_baseline = z1(temp(1)-FramesBL:temp(1)-1);
                z1 = z1([temp(1)-1; temp]);
                if numel(y) < Frames+1
                    if numel(y) == Frames
                        frameid = v.data(:, 2);
                        frameid = frameid([temp(1)-1; temp]);
                        frameid = frameid-frameid(1)+1;
                        lostframe = setdiff(1:Frames+1, frameid);
                        if numel(lostframe) ~= 1
                            continue;
                        end
                        if lostframe == Frames+1
                            y = [y; y(end)];
                            z1 = [z1; z1(end)];
                        else
                            yq = interp1(frameid, y, lostframe);
                            y = [y(1:lostframe-1); yq; y(lostframe:end)];
                            z1q = interp1(frameid, z1, lostframe);
                            z1 = [z1(1:lostframe-1); z1q; z1(lostframe:end)];
                        end
                    else
                        continue;
                    end
                elseif numel(y) > Frames+1
                    y = y(1:Frames+1);
                    z1 = z1(1:Frames+1);
                end
                % baseline
                if FramesBL > 1
                    temp = [y_baseline resolution(1)-z1_baseline];
                    temp = undistortPoints(temp, CP1);
                    temp = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), temp);
                    y_baseline = (temp(:, 1)-PG1REF(1))*inch2mmcst;
                    z1_baseline = (temp(:, 2)-PG1REF(2))*inch2mmcst;
                    y_baseline_all = [y_baseline_all y_baseline];
                    z1_baseline_all = [z1_baseline_all z1_baseline];
                end
                
                y_all_pixel = [y_all_pixel y];
                z1_all_pixel = [z1_all_pixel z1];
                temp = [y resolution(1)-z1];
                temp = undistortPoints(temp, CP1);
                temp = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), temp);
                y = (temp(:, 1)-PG1REF(1))*inch2mmcst;
                z1 = (temp(:, 2)-PG1REF(2))*inch2mmcst;
                y_all = [y_all y];
                z1_all = [z1_all z1];
                PG1repeat_all = [PG1repeat_all RepeatNumber(k)+RepeatNumberAll];
            catch
                disp(RepeatNumber(k)+RepeatNumberAll);
            end
            try
                v.data = load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
                if get(handles.LForelimb_radiobutton, 'Value')
                    load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG2camera' num2str(TrialNumber) '_trajectory.mat']);
                elseif get(handles.RForelimb_radiobutton, 'Value')
                    load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG2camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                elseif get(handles.Jaw_radiobutton, 'Value')
                    load([DataDir '\exp' num2str(ExpNumber(j), '%03d') '\PG2camera' num2str(TrialNumber) '_jaw_trajectory.mat']);
                end
                temp = find(v.data(:, 1) == 1);
                x = trajectory(:, 1);
                x_baseline = x(temp(1)-FramesBL:temp(1)-1);
                x = x([temp(1)-1; temp]);
                z2 = trajectory(:, 2);
                z2_baseline = z2(temp(1)-FramesBL:temp(1)-1);
                z2 = z2([temp(1)-1; temp]);
                if numel(x) < Frames+1
                    if numel(x) == Frames
                        frameid = v.data(:, 2);
                        frameid = frameid([temp(1)-1; temp]);
                        frameid = frameid-frameid(1)+1;
                        lostframe = setdiff(1:Frames+1, frameid);
                        if numel(lostframe) ~= 1
                            continue;
                        end
                        if lostframe == Frames+1
                            x = [x; x(end)];
                            z2 = [z2; z2(end)];
                        else
                            xq = interp1(frameid, x, lostframe);
                            x = [x(1:lostframe-1); xq; x(lostframe:end)];
                            z2q = interp1(frameid, z2, lostframe);
                            z2 = [z2(1:lostframe-1); z2q; z2(lostframe:end)];
                        end
                    else
                        continue;
                    end
                elseif numel(x) > Frames+1
                    x = x(1:Frames+1);
                    z2 = z2(1:Frames+1);
                end
                % baseline
                if FramesBL > 1
                    temp = [x_baseline resolution(1)-z2_baseline];
                    temp = undistortPoints(temp, CP2);
                    temp = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), temp);
                    x_baseline = (temp(:, 1)-PG2REF(1))*inch2mmcst;
                    z2_baseline = (temp(:, 2)-PG2REF(2))*inch2mmcst;
                    x_baseline_all = [x_baseline_all x_baseline];
                    z2_baseline_all = [z2_baseline_all z2_baseline];
                end
                
                x_all_pixel = [x_all_pixel x];
                z2_all_pixel = [z2_all_pixel z2];
                temp = [x resolution(1)-z2];
                temp = undistortPoints(temp, CP2);
                temp = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), temp);
                x = (temp(:, 1)-PG2REF(1))*inch2mmcst;
                z2 = (temp(:, 2)-PG2REF(2))*inch2mmcst;
                x_all = [x_all x];
                z2_all = [z2_all z2];
                PG2repeat_all = [PG2repeat_all RepeatNumber(k)+RepeatNumberAll];
            catch
                disp(RepeatNumber(k)+RepeatNumberAll);
            end
        end
        RepeatNumberAll = RepeatNumberAll+StimulusParameter.Repeat;
    end
    t = (-1:Frames-1)/FrameRate*1000;
    
    %1D trajectory
    if ~isempty(y_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1)
        py = zeros(1, size(y_all, 2));
        for j = 1:size(y_all, 2)
            py(j) = plot(t, y_all(:, j)-y_all(1, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
            set(py(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
            hold on;
        end
        if size(y_all, 2) > 1
            plot(t, mean(y_all, 2)-mean(y_all(1, :)), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('dY position (mm)');
        title(['N = ' num2str(size(y_all, 2)) ' | + Retraction; - Advance']);
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spy = zeros(1, size(y_all, 2));
        for j = 1:size(y_all, 2)
            spy(j) = plot(t, y_all(:, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
            set(spy(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
            hold on;
        end
        if size(y_all, 2) > 1
            plot(t, mean(y_all, 2), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('Y position (mm)');
        title(['N = ' num2str(size(y_all, 2)) ' | + Retraction; - Advance']);
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    if ~isempty(z1_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1);
        pz1 = zeros(1, size(z1_all, 2));
        for j = 1:size(z1_all, 2)
            pz1(j) = plot(t, z1_all(:, j)-z1_all(1, j), '-', 'Color', linecolors(round(j/size(z1_all, 2)*256), :), 'LineWidth', 2);
            set(pz1(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
            hold on;
        end
        if size(z1_all, 2) > 1
            plot(t, mean(z1_all, 2)-mean(z1_all(1, :)), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('dZ position (mm)');
        title(['N = ' num2str(size(z1_all, 2)) ' | + Elevation; - Depression']);
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spz1 = zeros(1, size(z1_all, 2));
        for j = 1:size(z1_all, 2)
            spz1(j) = plot(t, z1_all(:, j), '-', 'Color', linecolors(round(j/size(z1_all, 2)*256), :), 'LineWidth', 2);
            set(spz1(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
            hold on;
        end
        if size(z1_all, 2) > 1
            plot(t, mean(z1_all, 2), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('Z position (mm)');
        title(['N = ' num2str(size(z1_all, 2)) ' | + Elevation; - Depression']);
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    if ~isempty(x_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1);
        px = zeros(1, size(x_all, 2));
        for j = 1:size(x_all, 2)
            px(j) = plot(t, x_all(:, j)-x_all(1, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
            set(px(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(x_all, 2) > 1
            plot(t, mean(x_all, 2)-mean(x_all(1, :)), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('dX position (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | + Lateral; - Medial']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | Jaw: + Lateral; - Medial']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spx = zeros(1, size(x_all, 2));
        for j = 1:size(x_all, 2)
            spx(j) = plot(t, x_all(:, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
            set(spx(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(x_all, 2) > 1
            plot(t, mean(x_all, 2), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('X position (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | + Lateral; - Medial']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | Jaw: + Lateral; - Medial']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    if ~isempty(z2_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1);
        pz2 = zeros(1, size(z2_all, 2));
        for j = 1:size(z2_all, 2)
            pz2(j) = plot(t, z2_all(:, j)-z2_all(1, j), '-', 'Color', linecolors(round(j/size(z2_all, 2)*256), :), 'LineWidth', 2);
            set(pz2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(z2_all, 2) > 1
            plot(t, mean(z2_all, 2)-mean(z2_all(1, :)), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('dZ (backup) position (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | + Elevation; - Depression']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | Jaw: + Elevation; - Depression']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spz2 = zeros(1, size(z2_all, 2));
        for j = 1:size(z2_all, 2)
            spz2(j) = plot(t, z2_all(:, j), '-', 'Color', linecolors(round(j/size(z2_all, 2)*256), :), 'LineWidth', 2);
            set(spz2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(z2_all, 2) > 1
            plot(t, mean(z2_all, 2), '-k', 'LineWidth', 4);
        end
        box off;
        xlabel('Time after light on (ms)');
        ylabel('Z (backup) position (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | + Elevation; - Depression']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | Jaw: + Elevation; - Depression']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    
    % for publication start
    singlecolor = [0.85 0.85 0.85];
    if ~isempty(x_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1);
        px = zeros(1, size(x_all, 2));
        for j = 1:size(x_all, 2)
            px(j) = plot(t, x_all(:, j)-x_all(1, j), '-', 'Color', singlecolor, 'LineWidth', 1);
            set(px(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(x_all, 2) > 1
            plot(t, mean(x_all, 2)-mean(x_all(1, :)), '-k', 'LineWidth', 3);
        end
        box off;
        xlabel('Time (ms)');
        ylabel('dX (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | + Lateral; - Medial']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | Jaw: + Lateral; - Medial']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spx = zeros(1, size(x_all, 2));
        for j = 1:size(x_all, 2)
            spx(j) = plot(t, x_all(:, j), '-', 'Color', singlecolor, 'LineWidth', 1);
            set(spx(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(x_all, 2) > 1
            plot(t, mean(x_all, 2), '-k', 'LineWidth', 3);
        end
        box off;
        xlabel('Time (ms)');
        ylabel('X (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | + Lateral; - Medial']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(x_all, 2)) ' | Jaw: + Lateral; - Medial']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    if ~isempty(z2_all)
        figure('NumberTitle', 'off', 'Name', ['Spot ' num2str(SpotNumber(i))], 'Color', [1 1 1]);
        subplot(1, 2, 1);
        pz2 = zeros(1, size(z2_all, 2));
        for j = 1:size(z2_all, 2)
            pz2(j) = plot(t, z2_all(:, j)-z2_all(1, j), '-', 'Color', singlecolor, 'LineWidth', 1);
            set(pz2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(z2_all, 2) > 1
            plot(t, mean(z2_all, 2)-mean(z2_all(1, :)), '-k', 'LineWidth', 3);
        end
        box off;
        xlabel('Time (ms)');
        ylabel('dZ (backup) (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | + Elevation; - Depression']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | Jaw: + Elevation; - Depression']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
        subplot(1, 2, 2);
        spz2 = zeros(1, size(z2_all, 2));
        for j = 1:size(z2_all, 2)
            spz2(j) = plot(t, z2_all(:, j), '-', 'Color', singlecolor, 'LineWidth', 1);
            set(spz2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            hold on;
        end
        if size(z2_all, 2) > 1
            plot(t, mean(z2_all, 2), '-k', 'LineWidth', 3);
        end
        box off;
        xlabel('Time (ms)');
        ylabel('Z (backup) (mm)');
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | + Elevation; - Depression']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            title(['N = ' num2str(size(z2_all, 2)) ' | Jaw: + Elevation; - Depression']);
        end
        xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    end
    % for publication stop
    
%     p = 0.05;
%     delay = inf;
%     for j = 2:size(x_all, 1)
%         H = ttest(x_all(j, :), x_all(1, :), 'alpha', p);
%         if H == 1
%             delay = (j-2)/FrameRate*1000;
%             break;
%         end
%     end
%     for j = 2:size(z2_all, 1)
%         H = ttest(z2_all(j, :), z2_all(1, :), 'alpha', p);
%         if H == 1
%             delay = min(delay, (j-2)/FrameRate*1000);
%             break;
%         end
%     end
%     display(delay);
    
    %fft transform
    figure;
    subplot(2, 2, 1);
    fftpy = zeros(1, size(y_all, 2));
    for j = 1:size(y_all, 2)
        [f, MX] = FFT_trajectory(y_all(:, j), Fn);
        fftpy(j) = plot(f, MX, '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
        set(fftpy(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
    end
    box off;
    xlabel('Frequency (Hz)');
    ylabel('Power (mm^{2})');
    title(['N = ' num2str(size(y_all, 2)) ' | Y axis']);
    subplot(2, 2, 2);
    fftpz1 = zeros(1, size(z1_all, 2));
    for j = 1:size(z1_all, 2)
        [f, MX] = FFT_trajectory(z1_all(:, j), Fn);
        fftpz1(j) = plot(f, MX, '-', 'Color', linecolors(round(j/size(z1_all, 2)*256), :), 'LineWidth', 2);
        set(fftpz1(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
    end
    box off;
    xlabel('Frequency (Hz)');
    ylabel('Power (mm^{2})');
    title(['N = ' num2str(size(z1_all, 2)) ' | Z1 axis']);
    subplot(2, 2, 3);
    fftpx = zeros(1, size(x_all, 2));
    for j = 1:size(x_all, 2)
        [f, MX] = FFT_trajectory(x_all(:, j), Fn);
        fftpx(j) = plot(f, MX, '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
        set(fftpx(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
    end
    box off;
    xlabel('Frequency (Hz)');
    ylabel('Power (mm^{2})');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | X axis']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Jaw: X axis']);
    end
    subplot(2, 2, 4);
    fftpz2 = zeros(1, size(z2_all, 2));
    for j = 1:size(z2_all, 2)
        [f, MX] = FFT_trajectory(z2_all(:, j), Fn);
        fftpz2(j) = plot(f, MX, '-', 'Color', linecolors(round(j/size(z2_all, 2)*256), :), 'LineWidth', 2);
        set(fftpz2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
    end
    box off;
    xlabel('Frequency (Hz)');
    ylabel('Power (mm^{2})');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(z2_all, 2)) ' | Z2 axis']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(z2_all, 2)) ' | Jaw: Z2 axis']);
    end
    
    %2D trajectory
    figure;
    subplot(2, 2, 1);
    p2y = zeros(1, size(y_all, 2));
    for j = 1:size(y_all, 2)
        p2y(j) = plot(y_all(:, j)-y_all(1, j), z1_all(:, j)-z1_all(1, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
        set(p2y(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
        plot(y_all(end, j)-y_all(1, j), z1_all(end, j)-z1_all(1, j), 's', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
    end
    if size(y_all, 2) > 1
        plot(mean(y_all, 2)-mean(y_all(1, :)), mean(z1_all, 2)-mean(z1_all(1, :)), '-k', 'LineWidth', 2);
        plot(mean(y_all(end, :)-y_all(1, :)), mean(z1_all(end, :)-z1_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('dY position (mm)');
    ylabel('dZ position (mm)');
    title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
    subplot(2, 2, 2);
    sp2y = zeros(1, size(y_all, 2));
    for j = 1:size(y_all, 2)
        sp2y(j) = plot(y_all(:, j), z1_all(:, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
        set(sp2y(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
        plot(y_all(1, j), z1_all(1, j), 'o', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
        plot(y_all(end, j), z1_all(end, j), 's', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
    end
    if size(y_all, 2) > 1
        plot(mean(y_all, 2), mean(z1_all, 2), '-k', 'LineWidth', 2);
        plot(mean(y_all(1, :)), mean(z1_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
        plot(mean(y_all(end, :)), mean(z1_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('Y position (mm)');
    ylabel('Z position (mm)');
    title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
    subplot(2, 2, 3);
    p2x = zeros(1, size(x_all, 2));
    for j = 1:size(x_all, 2)
        p2x(j) = plot(x_all(:, j)-x_all(1, j), z2_all(:, j)-z2_all(1, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
        set(p2x(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
        plot(x_all(end, j)-x_all(1, j), z2_all(end, j)-z2_all(1, j), 's', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
    end
    if size(x_all, 2) > 1
        plot(mean(x_all, 2)-mean(x_all(1, :)), mean(z2_all, 2)-mean(z2_all(1, :)), '-k', 'LineWidth', 2);
        plot(mean(x_all(end, :)-x_all(1, :)), mean(z2_all(end, :)-z2_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('dX position (mm)');
    ylabel('dZ position (mm)');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Jaw: Front Camera']);
    end
    subplot(2, 2, 4);
    sp2x = zeros(1, size(x_all, 2));
    for j = 1:size(x_all, 2)
        sp2x(j) = plot(x_all(:, j), z2_all(:, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
        set(sp2x(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
        plot(x_all(1, j), z2_all(1, j), 'o', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
        plot(x_all(end, j), z2_all(end, j), 's', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
    end
    if size(x_all, 2) > 1
        plot(mean(x_all, 2), mean(z2_all, 2), '-k', 'LineWidth', 2);
        plot(mean(x_all(1, :)), mean(z2_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
        plot(mean(x_all(end, :)), mean(z2_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('X position (mm)');
    ylabel('Z position (mm)');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Jaw: Front Camera']);
    end
    
    % for publication start
    %2D trajectory
    singlecolor = [0.85 0.85 0.85];
    figure;
    subplot(2, 2, 1);
    p2y = zeros(1, size(y_all, 2));
    for j = 1:size(y_all, 2)
        p2y(j) = plot(y_all(:, j)-y_all(1, j), z1_all(:, j)-z1_all(1, j), '-', 'Color', singlecolor, 'LineWidth', 1);
        set(p2y(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
        plot(y_all(end, j)-y_all(1, j), z1_all(end, j)-z1_all(1, j), 's', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
    end
    if size(y_all, 2) > 1
        plot(mean(y_all, 2)-mean(y_all(1, :)), mean(z1_all, 2)-mean(z1_all(1, :)), '-k', 'LineWidth', 3);
        plot(mean(y_all(end, :)-y_all(1, :)), mean(z1_all(end, :)-z1_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('dY position (mm)');
    ylabel('dZ position (mm)');
    title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
    subplot(2, 2, 2);
    sp2y = zeros(1, size(y_all, 2));
    for j = 1:size(y_all, 2)
        sp2y(j) = plot(y_all(:, j), z1_all(:, j), '-', 'Color', singlecolor, 'LineWidth', 1);
        set(sp2y(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG1repeat_all(j)));
        hold on;
        plot(y_all(1, j), z1_all(1, j), 'o', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
        plot(y_all(end, j), z1_all(end, j), 's', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
    end
    if size(y_all, 2) > 1
        plot(mean(y_all, 2), mean(z1_all, 2), '-k', 'LineWidth', 3);
        plot(mean(y_all(1, :)), mean(z1_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
        plot(mean(y_all(end, :)), mean(z1_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('Y position (mm)');
    ylabel('Z position (mm)');
    title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
    subplot(2, 2, 3);
    p2x = zeros(1, size(x_all, 2));
    for j = 1:size(x_all, 2)
        p2x(j) = plot(x_all(:, j)-x_all(1, j), z2_all(:, j)-z2_all(1, j), '-', 'Color', singlecolor, 'LineWidth', 1);
        set(p2x(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
        plot(x_all(end, j)-x_all(1, j), z2_all(end, j)-z2_all(1, j), 's', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
    end
    if size(x_all, 2) > 1
        plot(mean(x_all, 2)-mean(x_all(1, :)), mean(z2_all, 2)-mean(z2_all(1, :)), '-k', 'LineWidth', 3);
        plot(mean(x_all(end, :)-x_all(1, :)), mean(z2_all(end, :)-z2_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('dX position (mm)');
    ylabel('dZ position (mm)');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Jaw: Front Camera']);
    end
    subplot(2, 2, 4);
    sp2x = zeros(1, size(x_all, 2));
    for j = 1:size(x_all, 2)
        sp2x(j) = plot(x_all(:, j), z2_all(:, j), '-', 'Color', singlecolor, 'LineWidth', 1);
        set(sp2x(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        hold on;
        plot(x_all(1, j), z2_all(1, j), 'o', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
        plot(x_all(end, j), z2_all(end, j), 's', 'Color', singlecolor, 'MarkerFaceColor', singlecolor, 'MarkerSize', markersize);
    end
    if size(x_all, 2) > 1
        plot(mean(x_all, 2), mean(z2_all, 2), '-k', 'LineWidth', 3);
        plot(mean(x_all(1, :)), mean(z2_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
        plot(mean(x_all(end, :)), mean(z2_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    end
    axis equal;
    box off;
    xlabel('X position (mm)');
    ylabel('Z position (mm)');
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);
    elseif get(handles.Jaw_radiobutton, 'Value')
        title(['N = ' num2str(size(x_all, 2)) ' | Jaw: Front Camera']);
    end
    % for publication stop
    
    commonrepeat = intersect(PG1repeat_all, PG2repeat_all);
    realrepeats = numel(commonrepeat);
    %3D trajectory
    figure;
    p3axish1 = subplot(1, 2, 1);
    common_xyz = zeros(size(x_all, 1), realrepeats, 3)/0;
    p3h1 = zeros(1, realrepeats);
    for j = 1:realrepeats
        p3h1(j) = plot3(x_all(:, PG2repeat_all == commonrepeat(j))-x_all(1, PG2repeat_all == commonrepeat(j)), ...
            y_all(:, PG1repeat_all == commonrepeat(j))-y_all(1, PG1repeat_all == commonrepeat(j)),...
            z1_all(:, PG1repeat_all == commonrepeat(j))-z1_all(1, PG1repeat_all == commonrepeat(j)),...
            '-', 'Color', linecolors(round(j/realrepeats*256), :), 'LineWidth', 2);
        set(p3h1(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, commonrepeat(j)));
        hold(gca, 'on');
        plot3(x_all(end, PG2repeat_all == commonrepeat(j))-x_all(1, PG2repeat_all == commonrepeat(j)),...
            y_all(end, PG1repeat_all == commonrepeat(j))-y_all(1, PG1repeat_all == commonrepeat(j)),...
            z1_all(end, PG1repeat_all == commonrepeat(j))-z1_all(1, PG1repeat_all == commonrepeat(j)),...
            's', 'Color', linecolors(round(j/realrepeats*256), :), 'MarkerSize', markersize, 'MarkerFaceColor', linecolors(round(j/realrepeats*256), :));
        common_xyz(:, j, 1) = x_all(:, PG2repeat_all == commonrepeat(j));
        common_xyz(:, j, 2) = y_all(:, PG1repeat_all == commonrepeat(j));
        common_xyz(:, j, 3) = z1_all(:, PG1repeat_all == commonrepeat(j));
    end
    axis equal;
    xlabel('dX position (mm)');
    ylabel('dY position (mm)');
    zlabel('dZ position (mm)');
    title(['N = ' num2str(realrepeats)]);
    
    p3axish2 = subplot(1, 2, 2);
    p3h2 = zeros(1, realrepeats);
    for j = 1:realrepeats
        p3h2(j) = plot3(x_all(:, PG2repeat_all == commonrepeat(j)), y_all(:, PG1repeat_all == commonrepeat(j)), z1_all(:, PG1repeat_all == commonrepeat(j)),...
            '-', 'Color', linecolors(round(j/realrepeats*256), :), 'LineWidth', 2);
        set(p3h2(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, commonrepeat(j)));
        hold(gca, 'on');
        plot3(x_all(1, PG2repeat_all == commonrepeat(j)), y_all(1, PG1repeat_all == commonrepeat(j)), z1_all(1, PG1repeat_all == commonrepeat(j)),...
            'o', 'Color', linecolors(round(j/realrepeats*256), :), 'MarkerSize', markersize, 'MarkerFaceColor', linecolors(round(j/realrepeats*256), :));
        plot3(x_all(end, PG2repeat_all == commonrepeat(j)), y_all(end, PG1repeat_all == commonrepeat(j)), z1_all(end, PG1repeat_all == commonrepeat(j)),...
            's', 'Color', linecolors(round(j/realrepeats*256), :), 'MarkerSize', markersize, 'MarkerFaceColor', linecolors(round(j/realrepeats*256), :));
    end
    axis equal;
    xlabel('X position (mm)');
    ylabel('Y position (mm)');
    zlabel('Z position (mm)');
    title(['N = ' num2str(realrepeats)]);
    set(p3axish1, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3axish2));
    set(p3axish2, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3axish1));
    
    % for publication start
    singlecolor = [0.85 0.85 0.85];
    avgcolor = [0 0 0];
    figure;
    for j = 1:realrepeats
        plot3(x_all(:, PG2repeat_all == commonrepeat(j)), y_all(:, PG1repeat_all == commonrepeat(j)), z1_all(:, PG1repeat_all == commonrepeat(j)),...
            '-', 'Color', singlecolor, 'LineWidth', 1);
        hold(gca, 'on');
        plot3(x_all(1, PG2repeat_all == commonrepeat(j)), y_all(1, PG1repeat_all == commonrepeat(j)), z1_all(1, PG1repeat_all == commonrepeat(j)),...
            'o', 'Color', singlecolor, 'MarkerSize', markersize, 'MarkerFaceColor', singlecolor);
        plot3(x_all(end, PG2repeat_all == commonrepeat(j)), y_all(end, PG1repeat_all == commonrepeat(j)), z1_all(end, PG1repeat_all == commonrepeat(j)),...
            's', 'Color', singlecolor, 'MarkerSize', markersize, 'MarkerFaceColor', singlecolor);
        x_all_temp(:, j) = x_all(:, PG2repeat_all == commonrepeat(j));
        y_all_temp(:, j) = y_all(:, PG1repeat_all == commonrepeat(j));
        z1_all_temp(:, j) = z1_all(:, PG1repeat_all == commonrepeat(j));
    end
    x_all_mean = mean(x_all_temp, 2);
    y_all_mean = mean(y_all_temp, 2);
    z1_all_mean = mean(z1_all_temp, 2);
    
    plot3(x_all_mean, y_all_mean, z1_all_mean,...
        '-', 'Color', avgcolor, 'LineWidth', 2);
    plot3(x_all_mean(1), y_all_mean(1), z1_all_mean(1),...
            'o', 'Color', avgcolor, 'MarkerSize', markersize, 'MarkerFaceColor', avgcolor);
        plot3(x_all_mean(end), y_all_mean(end), z1_all_mean(end),...
            's', 'Color', avgcolor, 'MarkerSize', markersize, 'MarkerFaceColor', avgcolor);
    axis equal;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title(['N = ' num2str(realrepeats)]);
    % for publication stop
    
    %SD plot
    figure;
    plot(t, std(x_all, 0, 2), '-r', 'LineWidth', 2);
    hold on;
    plot(t, std(y_all, 0, 2), '-g', 'LineWidth', 2);
    plot(t, std(z1_all, 0, 2), '-b', 'LineWidth', 2);
    common_xyzp(:, :, 1) = common_xyz(:, :, 1)-repmat(mean(common_xyz(:, :, 1), 2), 1, realrepeats);
    common_xyzp(:, :, 2) = common_xyz(:, :, 2)-repmat(mean(common_xyz(:, :, 2), 2), 1, realrepeats);
    common_xyzp(:, :, 3) = common_xyz(:, :, 3)-repmat(mean(common_xyz(:, :, 3), 2), 1, realrepeats);
    plot(t, mean(sqrt(sum(common_xyzp.^2, 3)), 2), '-k', 'LineWidth', 2);
    box off;
    xlabel('Time after light on (ms)');
    ylabel('SD (mm)');
    xlim([t(1)-1000/FrameRate t(end)+1000/FrameRate]);
    legend('X', 'Y', 'Z', 'All');
    legend boxoff;
    
    %speed plot
    flow = 5;
    fhigh = 15;
    p = 0.05;
    if get(handles.LForelimb_radiobutton, 'Value')
        % speed plot with baseline
        if FramesBL > 1
            speed_all = nan(Frames+FramesBL, realrepeats);
            speed = nan(Frames+FramesBL, 1);
            speedplot = zeros(1, realrepeats);
            t_with_baseline = (-FramesBL:Frames-1)/FrameRate*1000;
            figure;
            for j = 1:realrepeats
                distance_baseline = sqrt(diff(x_baseline_all(:, PG2repeat_all == commonrepeat(j))).^2+diff(y_baseline_all(:, PG1repeat_all == commonrepeat(j))).^2+diff(z1_baseline_all(:, PG1repeat_all == commonrepeat(j))).^2);
                distance = sqrt(diff(x_all(:, PG2repeat_all == commonrepeat(j))).^2+diff(y_all(:, PG1repeat_all == commonrepeat(j))).^2+diff(z1_all(:, PG1repeat_all == commonrepeat(j))).^2);
                speed(2:end) = [distance_baseline; distance]/(1/FrameRate);
                speedplot(j) = plot(t_with_baseline, speed, '-', 'Color', singlecolor, 'LineWidth', 1);
                hold on;
                set(speedplot(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, commonrepeat(j)));
                speed_all(:, j) = speed;
            end
            speedmean = mean(speed_all, 2);
            plot(t_with_baseline, speedmean, '-k', 'LineWidth', 2);
            yl = ylim;
            plot([0 0], yl, '--r');
            title(['N = ' num2str(realrepeats) ' | speed']);
            box off;
            xlabel('Time (ms)');
            ylabel('Speed (mm/s)');
            xlim([t_with_baseline(1) t_with_baseline(end)+1000/FrameRate]);
        end
        
        figure;
        subplot(1, 2, 1);
        speed_all = zeros(size(x_all, 1), realrepeats);
        speed = zeros(size(x_all, 1), 1)/0;
        speedplot = zeros(1, realrepeats);
        start2ref_all = zeros(1, realrepeats);
        end2ref_all = zeros(1, realrepeats);
        distance_all = zeros(1, realrepeats);
        distanceshort_all = zeros(1, realrepeats);
        frequency_all = zeros(1, realrepeats);
        power_all = zeros(1, realrepeats);
        for j = 1:realrepeats
            distance = sqrt(diff(x_all(:, PG2repeat_all == commonrepeat(j))).^2+diff(y_all(:, PG1repeat_all == commonrepeat(j))).^2+diff(z1_all(:, PG1repeat_all == commonrepeat(j))).^2);
            speed(2:end) = distance/(1/FrameRate);
            speedplot(j) = plot(t, speed, '-', 'Color', linecolors(round(j/realrepeats*256), :), 'LineWidth', 2);
            hold on;
            set(speedplot(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, commonrepeat(j)));
            speed_all(:, j) = speed;
            
            start2ref = sqrt(x_all(1, PG2repeat_all == commonrepeat(j))^2+y_all(1, PG1repeat_all == commonrepeat(j))^2+z1_all(1, PG1repeat_all == commonrepeat(j))^2);
            start2ref_all(j) = start2ref;
            end2ref = sqrt(x_all(end, PG2repeat_all == commonrepeat(j))^2+y_all(end, PG1repeat_all == commonrepeat(j))^2+z1_all(end, PG1repeat_all == commonrepeat(j))^2);
            end2ref_all(j) = end2ref;
            distance_all(j) = sum(distance);
            distanceshort = sqrt((x_all(end, PG2repeat_all == commonrepeat(j))-x_all(1, PG2repeat_all == commonrepeat(j)))^2+...
                (y_all(end, PG1repeat_all == commonrepeat(j))-y_all(1, PG1repeat_all == commonrepeat(j)))^2+...
                (z1_all(end, PG1repeat_all == commonrepeat(j))-z1_all(1, PG1repeat_all == commonrepeat(j)))^2);
            distanceshort_all(j) = distanceshort;
            
            [f, MX] = FFT_trajectory(x_all(:, PG2repeat_all == commonrepeat(j))-mean(x_all(:, PG2repeat_all == commonrepeat(j))), Fn);
            id = find(MX == max(MX));
            id = id(1);
            fmean = f(id);
            fpower = sum(MX(f >= flow & f <= fhigh));
            [f, MX] = FFT_trajectory(y_all(:, PG1repeat_all == commonrepeat(j))-mean(y_all(:, PG1repeat_all == commonrepeat(j))), Fn);
            id = find(MX == max(MX));
            id = id(1);
            fmean = f(id)+fmean;
            fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
            [f, MX] = FFT_trajectory(z1_all(:, PG1repeat_all == commonrepeat(j))-mean(z1_all(:, PG1repeat_all == commonrepeat(j))), Fn);
            id = find(MX == max(MX));
            id = id(1);
            fmean = f(id)+fmean;
            fmean = fmean/3;
            fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
            fpower = fpower/3;
            frequency_all(j) = fmean;
            power_all(j) = fpower;
        end
        start2ref_mean = mean(start2ref_all);
        end2ref_mean = mean(end2ref_all);
        distance_mean = mean(distance_all);
        distanceshort_mean = mean(distanceshort_all);
        straightnessindex = distanceshort_mean/distance_mean;
        frequency_mean = mean(frequency_all);
        power_mean = mean(power_all);
        speedmean = mean(speed_all, 2);
        peakspeed = max(speedmean);
        for j = 3:size(speed_all, 1)
            H = ttest(speed_all(j, :)-speed_all(2, :), 0, 'alpha', p, 'tail', 'right');
            if H == 1
                delay = (j-2)/FrameRate*1000;
                break;
            end
        end
        if H == 0
            delay = inf;
        end
        disp(['Start2Ref = ' num2str(start2ref_mean, '%.2f') ' mm']);
        disp(['End2Ref = ' num2str(end2ref_mean, '%.2f') ' mm']);
        disp(['Total Distance = ' num2str(distance_mean, '%.2f') ' mm']);
        disp(['Linear Distance = ' num2str(distanceshort_mean, '%.2f') ' mm']);
        disp(['Straightness Index = ' num2str(straightnessindex, '%.2f')]);
        disp(['Peakspeed = ' num2str(peakspeed, '%.2f') ' mm/s']);
        disp(['Delay = ' num2str(delay, '%.2f') ' ms']);
        disp(['Mean Frequency = ' num2str(frequency_mean, '%.2f') ' Hz']);
        disp(['Mean Power = ' num2str(power_mean, '%.2f') ' mm^{2}']);
        plot(t, speedmean, '-k', 'LineWidth', 4);
        title(['N = ' num2str(realrepeats) ' | speed']);
        box off;
        if size(speed_all, 2) == 1
            tstart = t(speedmean == max(speedmean));
            tstart = tstart(1);
            f = fittype('a*exp(-((x-b)/c)^2)+d', 'coefficients', {'a', 'b', 'c', 'd'},...
                'independent', 'x');
            options = fitoptions('Method', 'NonlinearLeastSquares',...
                'Robust', 'on',...
                'StartPoint', [max(speedmean) tstart 150 0],...
                'Lower', [0 0 0 0],...
                'Upper', [1500 500 inf 1500]);
            [curve, gof, ~] = fit(t(2:end)', speedmean(2:end), f, options);
            rsquare = gof.rsquare
            hfit = plot(curve, ':k');
            set(hfit, 'LineWidth', 2);
        end
        xlabel('Time after light on (ms)');
        ylabel('Speed (mm/s)');
        xlim([t(1) t(end)+1000/FrameRate]);
    elseif get(handles.Jaw_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        figure;
        subplot(1, 2, 1);
        speed_all = zeros(size(x_all));
        speed = zeros(size(x_all, 1), 1)/0;
        speedplot = zeros(1, size(x_all, 2));
        if get(handles.RForelimb_radiobutton, 'Value')
            start2ref_all = zeros(1, size(x_all, 2));
            end2ref_all = zeros(1, size(x_all, 2));
        end
        distance_all = zeros(1, size(x_all, 2));
        distanceshort_all = zeros(1, size(x_all, 2));
        frequency_all = zeros(1, size(x_all, 2));
        power_all = zeros(1, size(x_all, 2));
        for j = 1:length(speedplot)
            distance = sqrt(diff(x_all(:, j)).^2+diff(z2_all(:, j)).^2);
            speed(2:end) = distance/(1/FrameRate);
            speedplot(j) = plot(t, speed, '-', 'Color', linecolors(round(j/length(speedplot)*256), :), 'LineWidth', 2);
            hold on;
            set(speedplot(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
            speed_all(:, j) = speed;
            
            if get(handles.RForelimb_radiobutton, 'Value')
                start2ref = sqrt(x_all(1, j)^2+z2_all(1, j)^2);
                start2ref_all(j) = start2ref;
                end2ref = sqrt(x_all(end, j)^2+z2_all(end, j)^2);
                end2ref_all(j) = end2ref;
            end
            distance_all(j) = sum(distance);
            distanceshort = sqrt((x_all(end, j)-x_all(1, j))^2+(z2_all(end, j)-z2_all(1, j))^2);
            distanceshort_all(j) = distanceshort;
            
            [f, MX] = FFT_trajectory(z2_all(:, j)-mean(z2_all(:, j)), Fn);
            id = find(MX == max(MX));
            id = id(1);
            fmean = f(id);
            fpower = sum(MX(f >= flow & f <= fhigh));
            frequency_all(j) = fmean;
            power_all(j) = fpower;
        end
        if get(handles.RForelimb_radiobutton, 'Value')
            start2ref_mean = mean(start2ref_all);
            end2ref_mean = mean(end2ref_all);
        end
        distance_mean = mean(distance_all);
        distanceshort_mean = mean(distanceshort_all);
        straightnessindex = distanceshort_mean/distance_mean;
        frequency_mean = mean(frequency_all);
        power_mean = mean(power_all);
        speedmean = mean(speed_all, 2);
        peakspeed = max(speedmean);
        for j = 3:size(speed_all, 1)
            H = ttest(speed_all(j, :)-speed_all(2, :), 0, 'alpha', p, 'tail', 'right');
            if H == 1
                delay = (j-2)/FrameRate*1000;
                break;
            end
        end
        if H == 0
            delay = inf;
        end
        if get(handles.RForelimb_radiobutton, 'Value')
            disp(['Start2Ref = ' num2str(start2ref_mean, '%.2f') ' mm']);
            disp(['End2Ref = ' num2str(end2ref_mean, '%.2f') ' mm']);
        end
        disp(['Total Distance = ' num2str(distance_mean, '%.2f') ' mm']);
        disp(['Linear Distance = ' num2str(distanceshort_mean, '%.2f') ' mm']);
        disp(['Straightness Index = ' num2str(straightnessindex, '%.2f')]);
        disp(['Peakspeed = ' num2str(peakspeed, '%.2f') ' mm/s']);
        disp(['Delay = ' num2str(delay, '%.2f') ' ms']);
        disp(['Mean Frequency = ' num2str(frequency_mean, '%.2f') ' Hz']);
        disp(['Mean Power = ' num2str(power_mean, '%.2f') ' mm^{2}']);
        plot(t, speedmean, '-k', 'LineWidth', 4);
        if get(handles.Jaw_radiobutton, 'Value')
            title(['Jaw: N = ' num2str(length(speedplot)) ' | speed']);
        else
            title(['Right Forelimb: N = ' num2str(length(speedplot)) ' | speed']);
        end
        box off;
        if size(speed_all, 2) == 1
            tstart = t(speedmean == max(speedmean));
            tstart = tstart(1);
            f = fittype('a*exp(-((x-b)/c)^2)+d', 'coefficients', {'a', 'b', 'c', 'd'},...
                'independent', 'x');
            options = fitoptions('Method', 'NonlinearLeastSquares',...
                'Robust', 'on',...
                'StartPoint', [max(speedmean) tstart 150 0],...
                'Lower', [0 0 0 0],...
                'Upper', [1500 500 inf 1500]);
            [curve, gof, ~] = fit(t(2:end)', speedmean(2:end), f, options);
            rsquare = gof.rsquare
            hfit = plot(curve, ':k');
            set(hfit, 'LineWidth', 2);
        end
        xlabel('Time after light on (ms)');
        ylabel('Speed (mm/s)');
        xlim([t(1) t(end)+1000/FrameRate]);
    end
    
    %acceleration plot
    if get(handles.LForelimb_radiobutton, 'Value')
        subplot(1, 2, 2);
        acceleration_all = zeros(size(speed_all))/0;
        acceleration_all(2:end, :) = diff(speed_all)/(1/FrameRate);
        accelerationplot = zeros(1, realrepeats);
        for j = 1:realrepeats
            accelerationplot(j) = plot(t, acceleration_all(:, j), '-', 'Color', linecolors(round(j/realrepeats*256), :), 'LineWidth', 2);
            hold on;
            set(accelerationplot(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, commonrepeat(j)));
        end
        plot(t, mean(acceleration_all, 2), '-k', 'LineWidth', 4);
        xlabel('Time after light on (ms)');
        ylabel('Acceleration (mm/s^{2})');
        title(['N = ' num2str(realrepeats) ' | acceleration']);
        xlim([t(1)+1000/FrameRate t(end)+1000/FrameRate]);
        box off;
    elseif get(handles.Jaw_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        subplot(1, 2, 2);
        acceleration_all = zeros(size(speed_all))/0;
        acceleration_all(2:end, :) = diff(speed_all)/(1/FrameRate);
        accelerationplot = zeros(1, size(speed_all, 2));
        for j = 1:size(speed_all, 2)
            accelerationplot(j) = plot(t, acceleration_all(:, j), '-', 'Color', linecolors(round(j/size(speed_all, 2)*256), :), 'LineWidth', 2);
            hold on;
            set(accelerationplot(j), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, PG2repeat_all(j)));
        end
        plot(t, mean(acceleration_all, 2), '-k', 'LineWidth', 4);
        xlabel('Time after light on (ms)');
        ylabel('Acceleration (mm/s^{2})');
        if get(handles.Jaw_radiobutton, 'Value')
            title(['Jaw: N = ' num2str(size(speed_all, 2)) ' | acceleration']);
        else
            title(['Right Forelimb: N = ' num2str(size(speed_all, 2)) ' | acceleration']);
        end
        xlim([t(1)+1000/FrameRate t(end)+1000/FrameRate]);
        box off;
    end
end
if get(handles.LForelimb_radiobutton, 'Value')
    result.meanx = mean(common_xyz(:, :, 1), 2);
    result.meany = mean(common_xyz(:, :, 2), 2);
    result.meanz = mean(common_xyz(:, :, 3), 2);
    result.x_pixel = x_all_pixel;
    result.y_pixel = y_all_pixel;
    result.z1_pixel = z1_all_pixel;
    result.z2_pixel = z2_all_pixel;
    result.speedmean = speedmean;
    result.accelerationmean = mean(acceleration_all, 2);
    result.start2ref = start2ref_mean;
    result.end2ref = end2ref_mean;
    result.distance = distance_mean;
    result.distanceshort = distanceshort_mean;
    result.straightnessindex = straightnessindex;
    result.peakspeed = peakspeed;
    result.delay = delay;
    result.frequency = frequency_mean;
    result.power = power_mean;
elseif get(handles.Jaw_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    result.meanx = mean(x_all, 2);
    result.meanz = mean(z2_all, 2);
    result.speedmean = speedmean;
    result.accelerationmean = mean(acceleration_all, 2);
    if get(handles.RForelimb_radiobutton, 'Value')
        result.start2ref = start2ref_mean;
        result.end2ref = end2ref_mean;
        result.x_pixel = x_all_pixel;
        result.y_pixel = y_all_pixel;
        result.z1_pixel = z1_all_pixel;
        result.z2_pixel = z2_all_pixel;
    end
    result.distance = distance_mean;
    result.distanceshort = distanceshort_mean;
    result.straightnessindex = straightnessindex;
    result.peakspeed = peakspeed;
    result.delay = delay;
    result.frequency = frequency_mean;
    result.power = power_mean;
end
[file, path] = uiputfile('*.mat','Save');
if file ~= 0
	save([path file], 'result');
end

function ShowRepeatNumber(object, event, j)
disp(j);

function FramesBL_text_Callback(hObject, eventdata, handles)

function FramesBL_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FramesAL_text_Callback(hObject, eventdata, handles)

function FramesAL_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PG1Y0_text_Callback(hObject, eventdata, handles)

function PG1Y0_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PG1Z0_text_Callback(hObject, eventdata, handles)

function PG1Z0_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PG2X0_text_Callback(hObject, eventdata, handles)

function PG2X0_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PG2Z0_text_Callback(hObject, eventdata, handles)

function PG2Z0_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DeleteTrajectory_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = get(handles.Spot_list, 'Value');
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
for i = 1:length(ExpNumber)
    try
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat']);
    catch
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\exp' num2str(ExpNumber(i), '%03d') '.mat']);
        SP.Order = StimulusParameter.Order;
    end
    for j = 1:length(SpotNumber)
        for k = 1:length(RepeatNumber)
            TrialNumber = find(SP.Order == SpotNumber(j));
            TrialNumber = TrialNumber(RepeatNumber(k));
            if get(handles.LForelimb_radiobutton, 'Value')
                filename1 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_trajectory.mat'];
                filename2 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_trajectory.mat'];
            elseif get(handles.RForelimb_radiobutton, 'Value')
                filename1 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_rightpaw_trajectory.mat'];
                filename2 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_rightpaw_trajectory.mat'];
            elseif get(handles.Jaw_radiobutton, 'Value')
                filename2 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_jaw_trajectory.mat'];
            end
            if ~get(handles.Jaw_radiobutton, 'Value')
                try
                    delete(filename1);
                end
            end
            try
                delete(filename2);
            end
        end
    end
end

function ProcessTJ_button_Callback(hObject, eventdata, handles)
tic;
disp('Processing Started !');
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = get(handles.Spot_list, 'Value');
RepeatNumber = get(handles.Repeat_list, 'Value');
spot_string = get(handles.Spot_list, 'String');
if length(SpotNumber) ~= size(spot_string, 1)
    errordlg('You must select all the Spots !', 'ERROR');
    return;
end
PG1Y0 = str2double(get(handles.PG1Y0_text, 'String'));
PG1Z0 = str2double(get(handles.PG1Z0_text, 'String'));
PG2X0 = str2double(get(handles.PG2X0_text, 'String'));
PG2Z0 = str2double(get(handles.PG2Z0_text, 'String'));
if isnan(PG1Y0) || isnan(PG1Z0) || isnan(PG2X0) || isnan(PG2Z0)
    try
        load([DataDir '\ReferencePoint.mat']);
        PG1Y0 = Rf.PG1Y0;
        PG1Z0 = Rf.PG1Z0;
        PG2X0 = Rf.PG2X0;
        PG2Z0 = Rf.PG2Z0;
        set(handles.PG1Y0_text, 'String', num2str(PG1Y0));
        set(handles.PG1Z0_text, 'String', num2str(PG1Z0));
        set(handles.PG2X0_text, 'String', num2str(PG2X0));
        set(handles.PG2Z0_text, 'String', num2str(PG2Z0));
    catch
        errordlg('You need to set the Reference Point first!', 'Error');
        return;
    end
end
FrameRate = 100;
resolution = [512 640];
CP1 = load([DataDir '\Camera Calibration\Camera Calibration (PG1)\cameraParams.mat']);
CP1 = CP1.cameraParams;
CP2 = load([DataDir '\Camera Calibration\Camera Calibration (PG2)\cameraParams.mat']);
CP2 = CP2.cameraParams;
inch2mmcst = 25.4;
PG1REF = undistortPoints([PG1Y0 resolution(1)-PG1Z0], CP1);
PG1REF = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), PG1REF);
PG2REF = undistortPoints([PG2X0 resolution(1)-PG2Z0], CP1);
PG2REF = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), PG2REF);
distance_all = zeros(length(SpotNumber), length(ExpNumber)*length(RepeatNumber))/0;
distanceshort_all = zeros(size(distance_all))/0;
flow = 5;
fhigh = 15;
Fs = FrameRate;                             % sampling frequency
Fn = Fs/2;                                  % Nyquist frequency
frequency_all = zeros(size(distance_all))/0;
power_all = zeros(size(distance_all))/0;
pangle_all = zeros(size(distance_all))/0;
nangle_all = zeros(size(distance_all))/0;
endpoints3_all = zeros([size(distance_all) 3])/0;
startpoints3_all = zeros([size(distance_all) 3])/0;
startpoints_all = zeros([size(distance_all) 4])/0;
endpoints_all = zeros([size(distance_all) 4])/0;
end2ref_all = zeros(size(distance_all))/0;
speed_all = cell(1, length(SpotNumber));
trajectory_all = [];
if ~get(handles.Parfor_check, 'Value')
    workbar(0, 'Computing...', 'Progress');
end
for i = 1:length(ExpNumber)
    try
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat']);
        Frames = ceil(FrameRate*SP.Duration);
    catch
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\exp' num2str(ExpNumber(i), '%03d') '.mat']);
        SP.Order = StimulusParameter.Order;
        Frames = ceil(FrameRate*StimulusParameter.Duration);
    end
    if i == 1
        trajectory_all = zeros(length(ExpNumber)*length(RepeatNumber), Frames+1, length(SpotNumber), 4)/0;
    end
    tic;
    if get(handles.Parfor_check, 'Value')
        parfor j = 1:length(SpotNumber)
            startpoints_all_sliced = startpoints_all(j, :, :);
            endpoints_all_sliced = endpoints_all(j, :, :);
            trajectory_all_sliced = trajectory_all(:, :, j, :);
            pangle_all_sliced = pangle_all(j, :);
            nangle_all_sliced = nangle_all(j, :);
            distance_all_sliced = distance_all(j, :);
            distanceshort_all_sliced = distanceshort_all(j, :);
            endpoints3_all_sliced = endpoints3_all(j, :, :);
            startpoints3_all_sliced = startpoints3_all(j, :, :);
            end2ref_all_sliced = end2ref_all(j, :);
            frequency_all_sliced = frequency_all(j, :);
            power_all_sliced = power_all(j, :);
            for k = 1:length(RepeatNumber)
                trueid = (i-1)*length(RepeatNumber)+k;
                TrialNumber = find(SP.Order == SpotNumber(j));
                try
                    TrialNumber = TrialNumber(RepeatNumber(k));
                catch
                    continue;
                end
                y = [];
                z1 = [];
                x = [];
                z2 = [];
                v = struct();
                try
                    v.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
                    if get(handles.RForelimb_radiobutton, 'Value')
                        traj = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                    else
                        traj = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_trajectory.mat']);
                    end
                    temp = find(v.data(:, 1) == 1);
                    y = traj.trajectory(:, 1);
                    y = y([temp(1)-1; temp]);
                    z1 = traj.trajectory(:, 2);
                    z1 = z1([temp(1)-1; temp]);
                    if any(isnan(y)) || any(isnan(z1))
                        disp(['Exp' num2str(ExpNumber(i), '%03d') ' Spot' num2str(SpotNumber(j)) ' Repeat' num2str(RepeatNumber(k)) ' has a problem']);
                        continue;
                    end
                    if numel(y) < Frames+1
                        if numel(y) == Frames
                            frameid = v.data(:, 2);
                            frameid = frameid([temp(1)-1; temp]);
                            frameid = frameid-frameid(1)+1;
                            lostframe = setdiff(1:Frames+1, frameid);
                            if numel(lostframe) ~= 1
                                continue;
                            end
                            if lostframe == Frames+1
                                y = [y; y(end)];
                                z1 = [z1; z1(end)];
                            else
                                yq = interp1(frameid, y, lostframe);
                                y = [y(1:lostframe-1); yq; y(lostframe:end)];
                                z1q = interp1(frameid, z1, lostframe);
                                z1 = [z1(1:lostframe-1); z1q; z1(lostframe:end)];
                            end
                        else
                            continue;
                        end
                    elseif numel(y) > Frames+1
                        y = y(1:Frames+1);
                        z1 = z1(1:Frames+1);
                    end
                    temp = [y resolution(1)-z1];
                    temp = undistortPoints(temp, CP1);
                    temp = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), temp);
                    y = (temp(:, 1)-PG1REF(1))*inch2mmcst;
                    z1 = (temp(:, 2)-PG1REF(2))*inch2mmcst;
                    startpoints_all_sliced(1, trueid, 2) = y(1);
                    endpoints_all_sliced(1, trueid, 2) = y(end);
                    trajectory_all_sliced(trueid, :, 1, 2) = y;
                    startpoints_all_sliced(1, trueid, 3) = z1(1);
                    endpoints_all_sliced(1, trueid, 3) = z1(end);
                    trajectory_all_sliced(trueid, :, 1, 3) = z1;
                    pangle = 0;
                    nangle = 0;
                    ysmooth = smooth(y, 5);
                    z1smooth = smooth(z1, 5);
                    for kk = 1:length(y)-2
                        vector1 = [ysmooth(kk+1)-ysmooth(kk) z1smooth(kk+1)-z1smooth(kk)];
                        vector2 = [ysmooth(kk+2)-ysmooth(kk+1) z1smooth(kk+2)-z1smooth(kk+1)];
                        if sqrt((y(kk+2)-y(kk+1))^2+(z1(kk+2)-z1(kk+1))^2) < 1
                            continue;
                        end
                        if min(abs(vector1)) == 0 || min(abs(vector2)) == 0
                            continue;
                        end
                        rotate_angle = acosd(dot(vector1,vector2)/(norm(vector1)*norm(vector2)));
                        if vector1(1)*vector2(2)-vector1(2)*vector2(1) < 0
                            nangle = nangle+rotate_angle;
                        elseif vector1(1)*vector2(2)-vector1(2)*vector2(1) > 0
                            pangle = pangle+rotate_angle;
                        end
                    end
                    pangle_all_sliced(1, trueid) = pangle;
                    nangle_all_sliced(1, trueid) = nangle;
                end
                try
                    v.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
                    if get(handles.LForelimb_radiobutton, 'Value')
                        traj = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_trajectory.mat']);
                    elseif get(handles.RForelimb_radiobutton, 'Value')
                        traj = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                    elseif get(handles.Jaw_radiobutton, 'Value')
                        traj = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_jaw_trajectory.mat']);
                    end
                    temp = find(v.data(:, 1) == 1);
                    x = traj.trajectory(:, 1);
                    x = x([temp(1)-1; temp]);
                    z2 = traj.trajectory(:, 2);
                    z2 = z2([temp(1)-1; temp]);
                    if any(isnan(x)) || any(isnan(z2))
                        disp(['Exp' num2str(ExpNumber(i), '%03d') ' Spot' num2str(SpotNumber(j)) ' Repeat' num2str(RepeatNumber(k)) ' has a problem']);
                        continue;
                    end
                    if numel(x) < Frames+1
                        if numel(x) == Frames
                            frameid = v.data(:, 2);
                            frameid = frameid([temp(1)-1; temp]);
                            frameid = frameid-frameid(1)+1;
                            lostframe = setdiff(1:Frames+1, frameid);
                            if numel(lostframe) ~= 1
                                continue;
                            end
                            if lostframe == Frames+1
                                x = [x; x(end)];
                                z2 = [z2; z2(end)];
                            else
                                xq = interp1(frameid, x, lostframe);
                                x = [x(1:lostframe-1); xq; x(lostframe:end)];
                                z2q = interp1(frameid, z2, lostframe);
                                z2 = [z2(1:lostframe-1); z2q; z2(lostframe:end)];
                            end
                        else
                            continue;
                        end
                    elseif numel(x) > Frames+1
                        x = x(1:Frames+1);
                        z2 = z2(1:Frames+1);
                    end
                    temp = [x resolution(1)-z2];
                    temp = undistortPoints(temp, CP2);
                    temp = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), temp);
                    x = (temp(:, 1)-PG2REF(1))*inch2mmcst;
                    z2 = (temp(:, 2)-PG2REF(2))*inch2mmcst;
%                     if get(handles.Jaw_radiobutton, 'Value')
%                         if max(abs(z2-z2(1))) > 2 || max(abs(x-x(1))) > 1  %%%%%%%%%%%%%%%
%                             continue;
%                         end
%                     end
                    
%                     if get(handles.Jaw_radiobutton, 'Value')
%                         if max(abs(diff(z2))) > 1 || max(abs(diff(x))) > 1  %%%%%%%%%%%%%%%
%                             continue;
%                         end
%                     end
                    startpoints_all_sliced(1, trueid, 1) = x(1);
                    endpoints_all_sliced(1, trueid, 1) = x(end);
                    trajectory_all_sliced(trueid, :, 1, 1) = x;
                    startpoints_all_sliced(1, trueid, 4) = z2(1);
                    endpoints_all_sliced(1, trueid, 4) = z2(end);
                    trajectory_all_sliced(trueid, :, 1, 4) = z2;
                end
                if ~isempty(x)
                    if get(handles.LForelimb_radiobutton, 'Value') && ~isempty(y)
                        speed = zeros(1, length(x))/0;
                        speed(2:end) = sqrt(diff(x).^2+diff(y).^2+diff(z1).^2)/(1/FrameRate);
                        speed_all{j} = [speed_all{j}; speed];
                        distance = sum(sqrt(diff(x).^2+diff(y).^2+diff(z1).^2));
                        distance_all_sliced(1, trueid) = distance;
                        distanceshort_all_sliced(1, trueid) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z1(end)-z1(1))^2);
                        endpoints3_temp = [x(end) y(end) z1(end)];
                        startpoints3_temp = [x(1) y(1) z1(1)];
                        for kk = 1:3
                            endpoints3_all_sliced(1, trueid, kk) = endpoints3_temp(kk);
                            startpoints3_all_sliced(1, trueid, kk) = startpoints3_temp(kk);
                        end
                        end2ref_all_sliced(1, trueid) = sqrt(x(end)^2+y(end)^2+z1(end)^2);
                        [f, MX] = FFT_trajectory(x-mean(x), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id);
                        fpower = sum(MX(f >= flow & f <= fhigh));
                        [f, MX] = FFT_trajectory(y-mean(y), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id)+fmean;
                        fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
                        [f, MX] = FFT_trajectory(z1-mean(z1), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id)+fmean;
                        fmean = fmean/3;
                        fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
                        fpower = fpower/3;
                        frequency_all_sliced(1, trueid) = fmean;
                        power_all_sliced(1, trueid) = fpower;
                    elseif get(handles.Jaw_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
                        speed = zeros(1, length(x))/0;
                        speed(2:end) = sqrt(diff(x).^2+diff(z2).^2)/(1/FrameRate);
                        speed_all{j} = [speed_all{j}; speed];
                        distance = sum(sqrt(diff(x).^2+diff(z2).^2));
                        distance_all_sliced(1, trueid) = distance;
                        distanceshort_all_sliced(1, trueid) = sqrt((x(end)-x(1))^2+(z2(end)-z2(1))^2);
                        if get(handles.RForelimb_radiobutton, 'Value') && ~isempty(y)
                            endpoints3_temp = [x(end) y(end) z1(end)];
                            startpoints3_temp = [x(1) y(1) z1(1)];
                            for kk = 1:3
                                endpoints3_all_sliced(1, trueid, kk) = endpoints3_temp(kk);
                                startpoints3_all_sliced(1, trueid, kk) = startpoints3_temp(kk);
                            end
                            end2ref_all_sliced(1, trueid) = sqrt(x(end)^2+y(end)^2+z1(end)^2);
                        end
                        [f, MX] = FFT_trajectory(z2-mean(z2), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id);
                        fpower = sum(MX(f >= flow & f <= fhigh));
                        frequency_all_sliced(1, trueid) = fmean;
                        power_all_sliced(1, trueid) = fpower;
                    end
                end
            end
            startpoints_all(j, :, :) = startpoints_all_sliced;
            endpoints_all(j, :, :) = endpoints_all_sliced;
            trajectory_all(:, :, j, :) = trajectory_all_sliced;
            pangle_all(j, :) = pangle_all_sliced;
            nangle_all(j, :) = nangle_all_sliced;
            distance_all(j, :) = distance_all_sliced;
            distanceshort_all(j, :) = distanceshort_all_sliced;
            endpoints3_all(j, :, :) = endpoints3_all_sliced;
            startpoints3_all(j, :, :) = startpoints3_all_sliced;
            end2ref_all(j, :) = end2ref_all_sliced;
            frequency_all(j, :) = frequency_all_sliced;
            power_all(j, :) = power_all_sliced;
        end
    else
        for j = 1:length(SpotNumber)
            for k = 1:length(RepeatNumber)
                trueid = (i-1)*length(RepeatNumber)+k;
                TrialNumber = find(SP.Order == SpotNumber(j));
                try
                    TrialNumber = TrialNumber(RepeatNumber(k));
                catch
                    continue;
                end
                y = [];
                z1 = [];
                x = [];
                z2 = [];
                try
                    v.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
                    if get(handles.RForelimb_radiobutton, 'Value')
                        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                    else
                        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_trajectory.mat']);
                    end
                    temp = find(v.data(:, 1) == 1);
                    y = trajectory(:, 1);
                    y = y([temp(1)-1; temp]);
                    z1 = trajectory(:, 2);
                    z1 = z1([temp(1)-1; temp]);
                    if any(isnan(y)) || any(isnan(z1))
                        disp(['Exp' num2str(ExpNumber(i), '%03d') ' Spot' num2str(SpotNumber(j)) ' Repeat' num2str(RepeatNumber(k)) ' has a problem']);
                        continue;
                    end
                    if numel(y) < Frames+1
                        if numel(y) == Frames
                            frameid = v.data(:, 2);
                            frameid = frameid([temp(1)-1; temp]);
                            frameid = frameid-frameid(1)+1;
                            lostframe = setdiff(1:Frames+1, frameid);
                            if numel(lostframe) ~= 1
                                continue;
                            end
                            if lostframe == Frames+1
                                y = [y; y(end)];
                                z1 = [z1; z1(end)];
                            else
                                yq = interp1(frameid, y, lostframe);
                                y = [y(1:lostframe-1); yq; y(lostframe:end)];
                                z1q = interp1(frameid, z1, lostframe);
                                z1 = [z1(1:lostframe-1); z1q; z1(lostframe:end)];
                            end
                        else
                            continue;
                        end
                    elseif numel(y) > Frames+1
                        y = y(1:Frames+1);
                        z1 = z1(1:Frames+1);
                    end
                    temp = [y resolution(1)-z1];
                    temp = undistortPoints(temp, CP1);
                    temp = pointsToWorld(CP1, CP1.RotationMatrices(:, :, end), CP1.TranslationVectors(end, :), temp);
                    y = (temp(:, 1)-PG1REF(1))*inch2mmcst;
                    z1 = (temp(:, 2)-PG1REF(2))*inch2mmcst;
                    startpoints_all(j, trueid, 2) = y(1);
                    startpoints_all(j, trueid, 3) = z1(1);
                    endpoints_all(j, trueid, 2) = y(end);
                    endpoints_all(j, trueid, 3) = z1(end);
                    trajectory_all(trueid, :, j, 2) = y;
                    trajectory_all(trueid, :, j, 3) = z1;
                    pangle = 0;
                    nangle = 0;
                    ysmooth = smooth(y, 5);
                    z1smooth = smooth(z1, 5);
                    for kk = 1:length(y)-2
                        vector1 = [ysmooth(kk+1)-ysmooth(kk) z1smooth(kk+1)-z1smooth(kk)];
                        vector2 = [ysmooth(kk+2)-ysmooth(kk+1) z1smooth(kk+2)-z1smooth(kk+1)];
                        if sqrt((y(kk+2)-y(kk+1))^2+(z1(kk+2)-z1(kk+1))^2) < 1
                            continue;
                        end
                        if min(abs(vector1)) == 0 || min(abs(vector2)) == 0
                            continue;
                        end
                        rotate_angle = acosd(dot(vector1,vector2)/(norm(vector1)*norm(vector2)));
                        if vector1(1)*vector2(2)-vector1(2)*vector2(1) < 0
                            nangle = nangle+rotate_angle;
                        elseif vector1(1)*vector2(2)-vector1(2)*vector2(1) > 0
                            pangle = pangle+rotate_angle;
                        end
                    end
                    pangle_all(j, trueid) = pangle;
                    nangle_all(j, trueid) = nangle;
                end
                try
                    v.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
                    if get(handles.LForelimb_radiobutton, 'Value')
                        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_trajectory.mat']);
                    elseif get(handles.RForelimb_radiobutton, 'Value')
                        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_rightpaw_trajectory.mat']);
                    elseif get(handles.Jaw_radiobutton, 'Value')
                        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_jaw_trajectory.mat']);
                    end
                    temp = find(v.data(:, 1) == 1);
                    x = trajectory(:, 1);
                    x = x([temp(1)-1; temp]);
                    z2 = trajectory(:, 2);
                    z2 = z2([temp(1)-1; temp]);
                    if any(isnan(x)) || any(isnan(z2))
                        disp(['Exp' num2str(ExpNumber(i), '%03d') ' Spot' num2str(SpotNumber(j)) ' Repeat' num2str(RepeatNumber(k)) ' has a problem']);
                        continue;
                    end
                    if numel(x) < Frames+1
                        if numel(x) == Frames
                            frameid = v.data(:, 2);
                            frameid = frameid([temp(1)-1; temp]);
                            frameid = frameid-frameid(1)+1;
                            lostframe = setdiff(1:Frames+1, frameid);
                            if numel(lostframe) ~= 1
                                continue;
                            end
                            if lostframe == Frames+1
                                x = [x; x(end)];
                                z2 = [z2; z2(end)];
                            else
                                xq = interp1(frameid, x, lostframe);
                                x = [x(1:lostframe-1); xq; x(lostframe:end)];
                                z2q = interp1(frameid, z2, lostframe);
                                z2 = [z2(1:lostframe-1); z2q; z2(lostframe:end)];
                            end
                        else
                            continue;
                        end
                    elseif numel(x) > Frames+1
                        x = x(1:Frames+1);
                        z2 = z2(1:Frames+1);
                    end
                    temp = [x resolution(1)-z2];
                    temp = undistortPoints(temp, CP2);
                    temp = pointsToWorld(CP2, CP2.RotationMatrices(:, :, end), CP2.TranslationVectors(end, :), temp);
                    x = (temp(:, 1)-PG2REF(1))*inch2mmcst;
                    z2 = (temp(:, 2)-PG2REF(2))*inch2mmcst;
%                     if get(handles.Jaw_radiobutton, 'Value')
%                         if max(abs(z2-z2(1))) > 2 || max(abs(x-x(1))) > 1  %%%%%%%%%%%%%%%
%                             continue;
%                         end
%                     end
                    
%                     if get(handles.Jaw_radiobutton, 'Value')
%                         if max(abs(diff(z2))) > 1 || max(abs(diff(x))) > 1  %%%%%%%%%%%%%%%
%                             continue;
%                         end
%                     end
                    startpoints_all(j, trueid, 1) = x(1);
                    endpoints_all(j, trueid, 1) = x(end);
                    startpoints_all(j, trueid, 4) = z2(1);
                    endpoints_all(j, trueid, 4) = z2(end);
                    trajectory_all(trueid, :, j, 1) = x;
                    trajectory_all(trueid, :, j, 4) = z2;
                end
                if ~isempty(x)
                    if get(handles.LForelimb_radiobutton, 'Value') && ~isempty(y)
                        speed = zeros(1, length(x))/0;
                        speed(2:end) = sqrt(diff(x).^2+diff(y).^2+diff(z1).^2)/(1/FrameRate);
                        speed_all{j} = [speed_all{j}; speed];
                        distance = sum(sqrt(diff(x).^2+diff(y).^2+diff(z1).^2));
                        distance_all(j, trueid) = distance;
                        distanceshort_all(j, trueid) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z1(end)-z1(1))^2);
                        endpoints3_all(j, trueid, 1) = x(end);
                        endpoints3_all(j, trueid, 2) = y(end);
                        endpoints3_all(j, trueid, 3) = z1(end);
                        startpoints3_all(j, trueid, 1) = x(1);
                        startpoints3_all(j, trueid, 2) = y(1);
                        startpoints3_all(j, trueid, 3) = z1(1);
                        end2ref_all(j, (i-1)*length(RepeatNumber)+k) = sqrt(x(end)^2+y(end)^2+z1(end)^2);
                        [f, MX] = FFT_trajectory(x-mean(x), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id);
                        fpower = sum(MX(f >= flow & f <= fhigh));
                        [f, MX] = FFT_trajectory(y-mean(y), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id)+fmean;
                        fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
                        [f, MX] = FFT_trajectory(z1-mean(z1), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id)+fmean;
                        fmean = fmean/3;
                        fpower = sum(MX(f >= flow & f <= fhigh))+fpower;
                        fpower = fpower/3;
                        frequency_all(j, trueid) = fmean;
                        power_all(j, trueid) = fpower;
                    elseif get(handles.Jaw_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
                        speed = zeros(1, length(x))/0;
                        speed(2:end) = sqrt(diff(x).^2+diff(z2).^2)/(1/FrameRate);
                        speed_all{j} = [speed_all{j}; speed];
                        distance = sum(sqrt(diff(x).^2+diff(z2).^2));
                        distance_all(j, trueid) = distance;
                        distanceshort_all(j, trueid) = sqrt((x(end)-x(1))^2+(z2(end)-z2(1))^2);
                        if get(handles.RForelimb_radiobutton, 'Value') && ~isempty(y)
                            endpoints3_all(j, trueid, 1) = x(end);
                            endpoints3_all(j, trueid, 2) = y(end);
                            endpoints3_all(j, trueid, 3) = z1(end);
                            startpoints3_all(j, trueid, 1) = x(1);
                            startpoints3_all(j, trueid, 2) = y(1);
                            startpoints3_all(j, trueid, 3) = z1(1);
                            end2ref_all(j, (i-1)*length(RepeatNumber)+k) = sqrt(x(end)^2+y(end)^2+z1(end)^2);
                        end
                        [f, MX] = FFT_trajectory(z2-mean(z2), Fn);
                        id = find(MX == max(MX));
                        id = id(1);
                        fmean = f(id);
                        fpower = sum(MX(f >= flow & f <= fhigh));
                        frequency_all(j, trueid) = fmean;
                        power_all(j, trueid) = fpower;
                    end
                end
            end
            workbar(((i-1)*length(SpotNumber)+j)/length(ExpNumber)/length(SpotNumber), ['Exp' num2str(i) ' : Spot' num2str(j)], 'Progress');
        end
    end
    toc;
end
data.startpoints_all = startpoints_all;
data.endpoints_all = endpoints_all;
data.endpoints3_all = endpoints3_all;
data.startpoints3_all = startpoints3_all;
data.end2ref_all = end2ref_all;
data.trajectory_all = trajectory_all;
data.framerate = FrameRate;
data.speed_all = speed_all;
data.distance_all = distance_all;
data.distanceshort_all = distanceshort_all;
data.frequency_all = frequency_all;
data.power_all = power_all;
data.pangle_all = pangle_all;
data.nangle_all = nangle_all;
data.rows = StimulusParameter.Ysites;
data.columns = StimulusParameter.Xsites;

% uncomment below if you stimulate the left hemisphere and track the right forelimb
% mapping = reshape(1:128, 16, 8);
% mapping = fliplr(mapping);
% mapping = mapping(:);
% data.startpoints_all = data.startpoints_all(mapping, :, :);
% data.startpoints_all(:, :, 1) = -data.startpoints_all(:, :, 1);
% data.endpoints_all = data.endpoints_all(mapping, :, :);
% data.endpoints_all(:, :, 1) = -data.endpoints_all(:, :, 1);
% data.startpoints3_all = data.startpoints3_all(mapping, :, :);
% data.startpoints3_all(:, :, 1) = -data.startpoints3_all(:, :, 1);
% data.endpoints3_all = data.endpoints3_all(mapping, :, :);
% data.endpoints3_all(:, :, 1) = -data.endpoints3_all(:, :, 1);
% data.end2ref_all = data.end2ref_all(mapping, :);
% data.trajectory_all = data.trajectory_all(:, :, mapping, :);
% data.trajectory_all(:, :, :, 1) = -data.trajectory_all(:, :, :, 1);
% data.speed_all = data.speed_all(mapping);
% data.distance_all = data.distance_all(mapping, :);
% data.distanceshort_all = data.distanceshort_all(mapping, :);
% data.frequency_all = data.frequency_all(mapping, :);
% data.power_all = data.power_all(mapping, :);
% data.pangle_all = data.pangle_all(mapping, :);
% data.nangle_all = data.nangle_all(mapping, :);

set(handles.ProcessTJ_button, 'UserData', data);
% [file,path] = uiputfile('*.mat','Save');
% if file ~= 0
%     save([path file], 'data');
% end
disp('Processing Finished !');
toc;
msgbox('Processing Finished !');

function GenerateMaps_button_Callback(hObject, eventdata, handles)
scalefactor = 50;
data = get(handles.ProcessTJ_button, 'UserData');
if isempty(data)
    data = evalin('base', 'data');
end
%1D maps
colorbarim = ones(data.rows*scalefactor, 1*scalefactor, 3);
colorbarim(:, :, 1) = repmat(linspace(0.67, 0, data.rows*scalefactor).', 1, 1*scalefactor);
colorbarim(:, :, 3) = repmat(linspace(0, 1, 1*scalefactor), data.rows*scalefactor, 1);
im = zeros(data.rows, data.columns, 3)/0;
im(:, :, 2) = 1;
title_text = {'dx map' 'x end point map' 'dy map' 'y end point map' 'dz map' 'z end point map'};
for i = 1:3
    if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
        [d_mean, pvalue, end_mean, end_std] = takecareNaN(data.startpoints_all(:, :, i), data.endpoints_all(:, :, i));
    elseif get(handles.Jaw_radiobutton, 'Value')
        if i == 1
            [d_mean, pvalue, end_mean, end_std] = takecareNaN(data.startpoints_all(:, :, i), data.endpoints_all(:, :, i));
        elseif i == 2
            continue;
        elseif i == 3
            [d_mean, pvalue, end_mean, end_std] = takecareNaN(data.startpoints_all(:, :, 4), data.endpoints_all(:, :, 4));
        end
    end
    d_absmax = max(abs(d_mean));
    im(:, :, 1) = mat2gray(reshape(d_mean, data.rows, data.columns), [-d_absmax d_absmax])*0.67;
    im(:, :, 3) = 1-reshape(pvalue, data.rows, data.columns);
    figure;
    subplot(1, 2, 1);
    imshow(imresize(hsv2rgb(im), scalefactor, 'nearest'));
    title(title_text{i*2-1});
    cbaxish = subplot(1, 2, 2);
    imshow(hsv2rgb(colorbarim))
    axis on;
    box off;
    set(cbaxish, 'XTick', [0.5 1*scalefactor+0.5], 'XTickLabel', {'1' '0'})
    set(cbaxish, 'YTick', [0.5 0.5+data.rows*scalefactor/2 data.rows*scalefactor+0.5], 'YTickLabel', {num2str(d_absmax, '%.1f') '0' num2str(-d_absmax, '%.1f')})
    set(cbaxish, 'YAxisLocation', 'right')
    set(cbaxish, 'Position', [0.31 cbaxish.Position(2:4)]);
    
    im(:, :, 1) = mat2gray(reshape(end_mean, data.rows, data.columns))*0.67;
    im(:, :, 3) = 1-mat2gray(reshape(end_std, data.rows, data.columns));
    figure;
    subplot(1, 2, 1);
    imshow(imresize(hsv2rgb(im), scalefactor, 'nearest'));
    title(title_text{i*2});
    cbaxish = subplot(1, 2, 2);
    imshow(hsv2rgb(colorbarim))
    axis on;
    box off;
    set(cbaxish, 'XTick', [0.5 1*scalefactor+0.5], 'XTickLabel', {num2str(max(end_std), '%.1f') num2str(min(end_std), '%.1f')})
    set(cbaxish, 'YTick', [0.5 0.5+data.rows*scalefactor/2 data.rows*scalefactor+0.5], 'YTickLabel', {num2str(max(end_mean), '%.1f') num2str(max(end_mean)/2+min(end_mean)/2, '%.1f') num2str(min(end_mean), '%.1f')})
    set(cbaxish, 'YAxisLocation', 'right')
    set(cbaxish, 'Position', [0.31 cbaxish.Position(2:4)]);
end

%2D maps
d_mean = zeros(size(data.startpoints_all, 1), 4)/0;
for i = 1:4
    [d_mean(:, i), ~, ~, ~] = takecareNaN(data.startpoints_all(:, :, i), data.endpoints_all(:, :, i)); 
end
x = d_mean(:, 1);
y = d_mean(:, 2);
z1 = d_mean(:, 3);
z2 = d_mean(:, 4);
x = x/max(abs([x; z2]));
z2 = z2/max(abs([x; z2]));
y = y/max(abs([y; z1]));
z1 = z1/max(abs([y; z1]));
interval = 1;
[X, Y] = meshgrid(interval/2:interval:interval/2+interval*(data.columns-1), interval/2+interval*(data.rows-1):-interval:interval/2);
figure;
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    subplot(1, 2, 1);
end
hold on;
for i = 1:length(x)
    quiver(X(i), Y(i), x(i), z2(i), 'Color', 'black', 'LineWidth', 1, 'MaxHeadSize', 1);
end
axis image;
box on;
set(gca, 'XTick', [], 'YTick', [], 'XLim', [0 interval*data.columns], 'YLim', [0 interval*data.rows]);
xlabel('X');
ylabel('Z');
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    subplot(1, 2, 2);
    hold on;
    for i = 1:length(y)
        quiver(X(i), Y(i), y(i), z1(i), 'Color', 'black', 'LineWidth', 1, 'MaxHeadSize', 1);
    end
    axis image;
    box on;
    set(gca, 'XTick', [], 'YTick', [], 'XLim', [0 interval*data.columns], 'YLim', [0 interval*data.rows]);
    xlabel('Y');
    ylabel('Z');
end

%Scatter map
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    startpoints3 = zeros(size(data.endpoints3_all, 1), 4);
    endpoints3 = zeros(size(data.endpoints3_all, 1), 4);
    for i = 1:size(data.endpoints3_all, 1)
        x_row = data.startpoints3_all(i, :, 1);
        x_row(isnan(x_row)) = [];
        y_row = data.startpoints3_all(i, :, 2);
        y_row(isnan(y_row)) = [];
        z_row = data.startpoints3_all(i, :, 3);
        z_row(isnan(z_row)) = [];
        startpoints3(i, 1) = mean(x_row);
        startpoints3(i, 2) = mean(y_row);
        startpoints3(i, 3) = mean(z_row);
        
        x_row = data.endpoints3_all(i, :, 1);
        x_row(isnan(x_row)) = [];
        y_row = data.endpoints3_all(i, :, 2);
        y_row(isnan(y_row)) = [];
        z_row = data.endpoints3_all(i, :, 3);
        z_row(isnan(z_row)) = [];
        endpoints3(i, 1) = mean(x_row);
        endpoints3(i, 2) = mean(y_row);
        endpoints3(i, 3) = mean(z_row);
        endpoints3(i, 4) = mean(sqrt((x_row-mean(x_row)).^2+(y_row-mean(y_row)).^2+(z_row-mean(z_row)).^2));
    end
    figure;
    imshow(imresize(reshape(endpoints3(:, 4), data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Scatter map');
    colormap parula;
    colorbar;
    %3D plot
    temp = repmat(linspace(0, 1, data.columns), data.rows, 1);
    facecolor_all(:, 1) = temp(:);
    temp = repmat(linspace(0, 1, data.rows)', 1, data.columns);
    facecolor_all(:, 2) = temp(:);
    facecolor_all(:, 3) = 0;
    figure;
    subplot(1, 3, 1);
    imshow(imresize(reshape(facecolor_all, data.rows, data.columns, 3), scalefactor, 'nearest'), []);
    p3h1 = subplot(1, 3, 2);
    for i = 1:size(data.endpoints3_all, 1)
        plot3(endpoints3(i, 1)-startpoints3(i, 1), endpoints3(i, 2)-startpoints3(i, 2), endpoints3(i, 3)-startpoints3(i, 3), 'o', 'MarkerFaceColor',facecolor_all(i, :), 'MarkerEdgeColor', facecolor_all(i, :));
        hold on;
    end
    slocationh = plot3(0, 0, 0, 'o', 'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerSize', 12);
    axis equal;
    xlabel('dX position (mm)');
    ylabel('dY position (mm)');
    zlabel('dZ position (mm)');
    title('Distribution of end points in space');
    legend(slocationh, 'Rest Location');
    legend boxoff;
    p3h2 = subplot(1, 3, 3);
    for i = 1:size(data.endpoints3_all, 1)
        plot3(endpoints3(i, 1), endpoints3(i, 2), endpoints3(i, 3), 'o', 'MarkerFaceColor',facecolor_all(i, :), 'MarkerEdgeColor', facecolor_all(i, :));
        hold on;
    end
    for i = 1:3
        x0_all = data.startpoints_all(:, :, 1);
        x0_all(isnan(x0_all)) = [];
        y0_all = data.startpoints_all(:, :, 2);
        y0_all(isnan(y0_all)) = [];
        z0_all = data.startpoints_all(:, :, 3);
        z0_all(isnan(z0_all)) = [];
    end
    slocationh = plot3(mean(x0_all(:)), mean(y0_all(:)), mean(z0_all(:)), 'o', 'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerSize', 12);
    axis equal;
    xlabel('X position (mm)');
    ylabel('Y position (mm)');
    zlabel('Z position (mm)');
    title('Distribution of end points in space');
    legend(slocationh, 'Rest Location');
    legend boxoff;
    % colormap(gray)
    % colorbar('YTick', [0 0.5 1], 'YTickLabel', {num2str(min(endpoints3(:, 4)), '%.1f'), num2str((min(endpoints3(:, 4))+max(endpoints3(:, 4))/2), '%.1f'), num2str(max(endpoints3(:, 4)), '%.1f')});
    set(p3h1, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3h2));
    set(p3h2, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3h1));
    
    figure;
    p3th1 = subplot(1, 2, 1);
    p3ths = zeros(1, size(data.trajectory_all, 3));
    for i = 1:size(data.trajectory_all, 3)
        temp1 = [];
        temp2 = [];
        temp3 = [];
        for j = 1:size(data.trajectory_all, 1)
            if ~isnan(data.trajectory_all(j, 1, i, 1)) && ~isnan(data.trajectory_all(j, 1, i, 2)) && ~isnan(data.trajectory_all(j, 1, i, 3))
                temp1 = [temp1; data.trajectory_all(j, :, i, 1)-data.trajectory_all(j, 1, i, 1)];
                temp2 = [temp2; data.trajectory_all(j, :, i, 2)-data.trajectory_all(j, 1, i, 2)];
                temp3 = [temp3; data.trajectory_all(j, :, i, 3)-data.trajectory_all(j, 1, i, 3)];
            end
        end
        p3ths(i) = plot3(mean(temp1), mean(temp2), mean(temp3), '-', 'Color',facecolor_all(i, :), 'LineWidth', 2);
        set(p3ths(i), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, i));
        hold on;
        plot3(mean(temp1(:, end)), mean(temp2(:, end)), mean(temp3(:, end)), 's', 'Color', facecolor_all(i, :), 'MarkerSize', 8, 'MarkerFaceColor', facecolor_all(i, :));
    end
    xlabel('dX (mm)');
    ylabel('dY (mm)');
    zlabel('dZ (mm)');
    p3th2 = subplot(1, 2, 2);
    p3th = zeros(1, size(data.trajectory_all, 3));
    for i = 1:size(data.trajectory_all, 3)
        temp1 = [];
        temp2 = [];
        temp3 = [];
        for j = 1:size(data.trajectory_all, 1)
            if ~isnan(data.trajectory_all(j, 1, i, 1)) && ~isnan(data.trajectory_all(j, 1, i, 2)) && ~isnan(data.trajectory_all(j, 1, i, 3))
                temp1 = [temp1; data.trajectory_all(j, :, i, 1)];
                temp2 = [temp2; data.trajectory_all(j, :, i, 2)];
                temp3 = [temp3; data.trajectory_all(j, :, i, 3)];
            end
        end
        p3th(i) = plot3(mean(temp1), mean(temp2), mean(temp3), '-', 'Color',facecolor_all(i, :), 'LineWidth', 2);
        set(p3th(i), 'ButtonDownFcn', @(object, event)ShowRepeatNumber(object, event, i));
        hold on;
        plot3(mean(temp1(:, 1)), mean(temp2(:, 1)), mean(temp3(:, 1)), 'o', 'Color', facecolor_all(i, :), 'MarkerSize', 8, 'MarkerFaceColor', facecolor_all(i, :));
        plot3(mean(temp1(:, end)), mean(temp2(:, end)), mean(temp3(:, end)), 's', 'Color', facecolor_all(i, :), 'MarkerSize', 8, 'MarkerFaceColor', facecolor_all(i, :));
    end
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    set(p3th1, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3th2));
    set(p3th2, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3th1));
end

%distance maps
distance = zeros(size(data.distance_all, 1), 1);
distanceshort = distance;
for i = 1:size(data.distance_all, 1)
    distance_row = data.distance_all(i, :);
    distance_row(isnan(distance_row)) = [];
    distance(i) = mean(distance_row);
    distanceshort_row = data.distanceshort_all(i, :);
    distanceshort_row(isnan(distanceshort_row)) = [];
    distanceshort(i) = mean(distanceshort_row);
end
figure;
subplot(1, 3, 1);
imshow(imresize(reshape(distanceshort, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Linear distance map');
colormap(viridis);
colorbar;
subplot(1, 3, 2);
imshow(imresize(reshape(distance, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Travel distance map');
colormap(viridis);
colorbar;
subplot(1, 3, 3);
imshow(imresize(reshape(distanceshort./distance, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Straightness index map');
colormap(viridis);
colorbar;

%end2ref map
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    end2ref = zeros(size(data.end2ref_all, 1), 1);
    for i = 1:size(data.end2ref_all, 1)
        end2ref_row = data.end2ref_all(i, :);
        end2ref_row(isnan(end2ref_row)) = [];
        end2ref(i) = mean(end2ref_row);
    end
    figure;
    imshow(imresize(reshape(end2ref, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('End2Ref map');
    colormap(viridis);
    colorbar;
end

%fft maps
frequency = zeros(size(data.frequency_all, 1), 1);
pw = frequency;
for i = 1:size(data.frequency_all, 1)
    frequency_row = data.frequency_all(i, :);
    frequency_row(isnan(frequency_row)) = [];
    frequency(i) = mean(frequency_row);
    pw_row = data.power_all(i, :);
    pw_row(isnan(pw_row)) = [];
    pw(i) = mean(pw_row);
end
figure;
subplot(1, 2, 1);
imshow(imresize(reshape(frequency, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Frequency map');
colormap parula;
colorbar;
subplot(1, 2, 2);
imshow(imresize(reshape(pw, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Rhythmicity map');
colormap parula;
colorbar;

%rotation angle maps
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    pangle = zeros(size(data.pangle_all, 1), 1);
    nangle = pangle;
    sd_pangle = pangle;
    sd_nangle = nangle;
    for i = 1:size(data.pangle_all, 1)
        pangle_row = data.pangle_all(i, :);
        pangle_row(isnan(pangle_row)) = [];
        pangle(i) = mean(pangle_row);
        sd_pangle(i) = std(pangle_row);
        nangle_row = data.nangle_all(i, :);
        nangle_row(isnan(nangle_row)) = [];
        nangle(i) = mean(nangle_row);
        sd_nangle(i) = std(nangle_row);
    end
    figure;
    subplot(1, 4, 1);
    imshow(imresize(reshape(pangle, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Forward rotation angle map');
    colormap(gca, viridis(256, 'viridis'));
    colorbar;
    subplot(1, 4, 2);
    imshow(imresize(reshape(nangle, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Backward rotation angle map');
    colormap(gca, viridis(256, 'viridis'));
    colorbar;
    subplot(1, 4, 3);
    imshow(imresize(reshape(mat2gray(sd_pangle), data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Normalized SD map of forward rotation angle');
    colormap(gca, 'gray');
    colorbar;
    subplot(1, 4, 4);
    imshow(imresize(reshape(mat2gray(sd_nangle), data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Normalized SD map of backward rotation angle');
    colormap(gca, 'gray');
    colorbar;
end

%peak speed map
speed_all = data.speed_all;
peakspeed = zeros(1, length(speed_all));
p = 0.05;
FrameRate = data.framerate;
delay = zeros(1, length(speed_all));
delay_rmANOVA = zeros(1, length(speed_all));
parfor i = 1:length(speed_all)
    speed = speed_all{i};
    peakspeed(i) = max(mean(speed_all{i}));
    H = 0;
    for j = 3:size(speed, 2)
%         H = ttest(speed(:, j)-speed(:, 2), 0, 'alpha', p, 'tail', 'right');
        [~, H] = signrank(speed(:, j), speed(:, 2), 'alpha', p, 'tail', 'right');
        if H == 1
            delay(i) = (j-2)/FrameRate*1000;
            break;
        end
    end
    if H == 0
        delay(i) = NaN;
    end
    
    speed(:, 1) = [];
    t = array2table(speed);
    rm = fitrm(t, ['speed1-speed' num2str(size(speed, 2)) '~ 1']);
    tbl = multcompare(rm, 'Time');
    tbl = table2array(tbl);
    tbl(tbl(:, 1)~=1, :) = [];
    id = find(tbl(:, 5) < p, 1);
    if isempty(id)
        delay_rmANOVA(i) = NaN;
    else
        delay_rmANOVA(i) = id/FrameRate*1000;
    end
end
figure;
imshow(imresize(reshape(peakspeed, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Peak speed map');
colormap(gca, 'parula');
colorbar;

%delay map
figure;
if any(~isnan(delay))
    subplot(1, 2, 1);
    alpha = imresize(reshape(delay, data.rows, data.columns), scalefactor, 'nearest');
    alpha(~isnan(alpha)) = 1;
    alpha(isnan(alpha)) = 0;
    imshow(zeros([size(alpha) 3]), []);
    hold on;
    h = imshow(imresize(reshape(delay, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Delay map (Sign Rank Test)');
    colormap(gca, 'parula');
    colorbar;
    set(h, 'AlphaData', alpha);
end

if any(~isnan(delay_rmANOVA))
    subplot(1, 2, 2);
    alpha = imresize(reshape(delay_rmANOVA, data.rows, data.columns), scalefactor, 'nearest');
    alpha(~isnan(alpha)) = 1;
    alpha(isnan(alpha)) = 0;
    imshow(zeros([size(alpha) 3]), []);
    hold on;
    h = imshow(imresize(reshape(delay_rmANOVA, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Delay map (rmANOVA)');
    colormap(gca, 'parula');
    colorbar;
    set(h, 'AlphaData', alpha);
end

data.peakspeed = peakspeed;
data.delay = delay;
[file,path] = uiputfile('*.mat','Save');
if file ~= 0
    save([path file], 'data');
end

%discreteness index map
if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
    trajectory = data.trajectory_all;
    t = (0:size(trajectory, 2)-1)/FrameRate*1000;
    t(1) = [];
    rsquare = zeros(size(trajectory, 3), 1);
    f = fittype('a*exp(-((x-b)/c)^2)+d', 'coefficients', {'a', 'b', 'c', 'd'},...
        'independent', 'x');
    workbar(0, 'Computing...', 'Progress');
    for i = 1:size(trajectory, 3)
        R2 = [];
        for j = 1:size(trajectory, 1)
            if ~isnan(trajectory(j, 1, i, 1)) && ~isnan(trajectory(j, 1, i, 2)) && ~isnan(trajectory(j, 1, i, 3))
                distance = sqrt(diff(trajectory(j, :, i, 1)).^2+diff(trajectory(j, :, i, 2)).^2+diff(trajectory(j, :, i, 3)).^2);
                speed = distance/(1/FrameRate);
                tstart = t(speed == max(speed));
                tstart = tstart(1);
                options = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Robust', 'on',...
                    'StartPoint', [max(speed) tstart 150 0],...
                    'Lower', [0 0 0 0],...
                    'Upper', [inf 500 inf inf]);
                [~, gof, ~] = fit(t', speed', f, options);
                R2 = [R2 gof.rsquare];
            end
        end
        if length(R2) > 1
            rsquare(i) = mean(R2);
        else
            rsquare(i) = R2;
        end
        workbar(i/size(trajectory,3), 'Computing...', 'Progress');
    end
    figure;
    imshow(imresize(reshape(rsquare, data.rows, data.columns), scalefactor, 'nearest'), []);
    title('Discreteness index map');
    colormap(gca, 'parula');
    colorbar;
end
    

function SaveRP_button_Callback(hObject, eventdata, handles)
PG1Y0 = str2double(get(handles.PG1Y0_text, 'String'));
if ~isnan(PG1Y0)
    Rf.PG1Y0 = PG1Y0;
else
    errordlg('PG1 Y0, PG1 Z0, PG2 X0, PG2 Z0 all need to have a value !', 'ERROR');
    return;
end
PG1Z0 = str2double(get(handles.PG1Z0_text, 'String'));
if ~isnan(PG1Z0)
    Rf.PG1Z0 = PG1Z0;
else
    errordlg('PG1 Y0, PG1 Z0, PG2 X0, PG2 Z0 all need to have a value !', 'ERROR');
    return;
end
PG2X0 = str2double(get(handles.PG2X0_text, 'String'));
if ~isnan(PG2X0)
    Rf.PG2X0 = PG2X0;
else
    errordlg('PG1 Y0, PG1 Z0, PG2 X0, PG2 Z0 all need to have a value !', 'ERROR');
    return;
end
PG2Z0 = str2double(get(handles.PG2Z0_text, 'String'));
if ~isnan(PG2Z0)
    Rf.PG2Z0 = PG2Z0;
else
    errordlg('PG1 Y0, PG1 Z0, PG2 X0, PG2 Z0 all need to have a value !', 'ERROR');
    return;
end
DataDir = get(handles.dir_text, 'String');
if exist([DataDir '\ReferencePoint.mat'], 'file')
    user_response = modaldlg;
    switch lower(user_response)
        case 'no'
            % take no action
            return;
        case 'yes'
    end
end
save([DataDir '\ReferencePoint.mat'], 'Rf');
disp('Done');

function MarkerThreshold_text_Callback(hObject, eventdata, handles)

function MarkerThreshold_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Autoprocessing_button_Callback(hObject, eventdata, handles)
tic;
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = get(handles.Spot_list, 'Value');
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
FramesBL = str2double(get(handles.FramesBL_text, 'String'));
FramesAL = str2double(get(handles.FramesAL_text, 'String'));
v.FrameRange = [FramesBL FramesAL];
v.threshold = str2double(get(handles.MarkerThreshold_text, 'String'));
v.threshold = v.threshold/100;
v.hue(1) = str2double(get(handles.huelow_text, 'String'));
v.hue(2) = str2double(get(handles.huehigh_text, 'String'));
filtersize = 9;
v.hfilter = fspecial('gaussian', filtersize, filtersize/5);
v.se = strel('disk', 1);
v.hblob = vision.BlobAnalysis('AreaOutputPort', true, 'CentroidOutputPort', true);
v.method = 1;
v.width = 20;
v.pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
v.thresholdcheck = get(handles.Thresholding_check, 'Value');

try
    load([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\VRDLPParameter.mat']);
catch
    load([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\exp' num2str(ExpNumber(1), '%03d') '.mat']);
    SP.Order = StimulusParameter.Order;
end
TrialNumber = find(SP.Order == SpotNumber(1));
TrialNumber = TrialNumber(RepeatNumber(1));
vd.data = load([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
lightframeid = find(vd.data(:, 1) == 1);
skipframes = lightframeid(1)-v.FrameRange(1)-1;
videoinfo = info(vd.FileReader);
v.FrameRate = videoinfo.VideoFrameRate;
v.Height = videoinfo.VideoSize(2);
v.Width = videoinfo.VideoSize(1);
for i = 1:skipframes
    step(vd.FileReader);
end
im = step(vd.FileReader);
im_temp = rgb2hsv(im);
im_temp = imfilter(im_temp, v.hfilter, 'same');
im_v = im_temp(:, :, 3);
im_h = im_temp(:, :, 1);
BW_h = im_h > v.hue(1) & im_h < v.hue(2);
imtool(im_h);
imtool(im_v);
figure;
imshow(im, []);
h = drawfreehand;
BW1 = createMask(h);
v.BW1 = BW1;
close;
ROI = im_v(BW1&BW_h);
ROI = sort(ROI);
BWmarker = im_v >= ROI(floor(v.threshold*numel(ROI))+1);
BWmarker = BWmarker&BW1&BW_h;
BWmarker = imerode(BWmarker, v.se);
BWmarker = imdilate(BWmarker, v.se);
if get(handles.Jaw_radiobutton, 'Value')
    [row, column] = find(BWmarker == 1);
    center = [mean(column(row == max(row))) max(row)];
elseif get(handles.LForelimb_radiobutton, 'Value')
    [area, centroid] = step(v.hblob, BWmarker);
    center = centroid(area == max(area), :);
end
simage = insertShape(im, 'FilledCircle', [center(1, :) 5], 'Color', [0 1 0], 'Opacity', 1);
if vd.data(skipframes+i, 1)
    simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
end
release(vd.FileReader);
if ~get(handles.Jaw_radiobutton, 'Value')
    vd.data = load([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
    vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(1), '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
    lightframeid = find(vd.data(:, 1) == 1);
    skipframes = lightframeid(1)-v.FrameRange(1)-1;
    for i = 1:skipframes
        step(vd.FileReader);
    end
    im = step(vd.FileReader);
    im_temp = rgb2hsv(im);
    im_temp = imfilter(im_temp, v.hfilter, 'same');
    im_v = im_temp(:, :, 3);
    im_h = im_temp(:, :, 1);
    imtool(im_h);
    imtool(im_v);
    figure;
    imshow(im, []);
    h = drawfreehand;
    BW2 = createMask(h);
    v.BW2 = BW2;
    close;
    release(vd.FileReader);
end
if get(handles.Play_checkbox, 'Value')
    figure;
    subplot(1, 2, 1);
    him(1) = imshow(BWmarker, []);
    subplot(1, 2, 2);
    him(2) = imshow(simage, []);
    v.him = him;
end

disp('Tracking Started !');
if ~get(handles.Parfor_check, 'Value')
    n = 1;
    workbar(0, 'Computing Ongoing...', 'Progress');
end
for i = 1:length(ExpNumber)
    try
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat']);
    catch
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\exp' num2str(ExpNumber(i), '%03d') '.mat']);
        SP.Order = StimulusParameter.Order;
    end
    if get(handles.Parfor_check, 'Value')
        parfor j = 1:length(SpotNumber)
            for k = 1:length(RepeatNumber)
                vd = struct();
                try
                    TrialNumber = find(SP.Order == SpotNumber(j));
                    TrialNumber = TrialNumber(RepeatNumber(k));
                catch
                    continue;
                end
                try
                    vd.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
                    vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
                    vd.FileName = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber)];
                    vd.thresholdcheck = v.thresholdcheck;
                    tracking_state = auto_tracking_parfor(v, vd, 1, handles);
                    if tracking_state ~= 1 && ~v.thresholdcheck
                        vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
                        vd.thresholdcheck = 1;
                        auto_tracking_parfor(v, vd, 1, handles);
                    end
                end
                if ~get(handles.Jaw_radiobutton, 'Value')
                    try
                        vd.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
                        vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
                        vd.FileName = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber)];
                        vd.thresholdcheck = v.thresholdcheck;
                        tracking_state = auto_tracking_parfor(v, vd, 2, handles);
                        if tracking_state ~= 1 && ~v.thresholdcheck
                            vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
                            vd.thresholdcheck = 1;
                            auto_tracking_parfor(v, vd, 2, handles);
                        end
                    end
                end
            end
        end
    else
        for j = 1:length(SpotNumber)
            for k = 1:length(RepeatNumber)
                try
                    TrialNumber = find(SP.Order == SpotNumber(j));
                    TrialNumber = TrialNumber(RepeatNumber(k));
                catch
                    workbar(n/length(RepeatNumber)/length(SpotNumber)/length(ExpNumber), 'Computing Ongoing...', 'Progress'); 
                    n = n+1;
                    continue;
                end
                try
                    vd.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2gpio' num2str(TrialNumber) '.csv']);
                    vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
                    vd.FileName = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber)];
                    vd.thresholdcheck = v.thresholdcheck;
                    tracking_state = auto_tracking(v, vd, 1, handles);
                    if tracking_state ~= 1 && ~v.thresholdcheck
                        vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '.avi']);
                        vd.thresholdcheck = 1;
                        auto_tracking(v, vd, 1, handles);
                    end
                end
                if ~get(handles.Jaw_radiobutton, 'Value')
                    try
                        vd.data = load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1gpio' num2str(TrialNumber) '.csv']);
                        vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
                        vd.FileName = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber)];
                        vd.thresholdcheck = v.thresholdcheck;
                        tracking_state = auto_tracking(v, vd, 2, handles);
                        if tracking_state ~= 1 && ~v.thresholdcheck
                            vd.FileReader = vision.VideoFileReader([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '.avi']);
                            vd.thresholdcheck = 1;
                            auto_tracking(v, vd, 2, handles);
                        end
                    end
                end
                workbar(n/length(RepeatNumber)/length(SpotNumber)/length(ExpNumber), 'Computing Ongoing...', 'Progress'); 
                n = n+1;
            end
        end
    end
end
if get(handles.Parfor_check, 'Value')
    msgbox('Tracking Finished !');
end
disp('Tracking Finished !');
toc;

function tracking_state = auto_tracking_parfor(v, vd, callerid, handles)
tracking_state = 1;
lightframeid = find(vd.data(:, 1) == 1);
Frames = v.FrameRange(1)+numel(lightframeid)+v.FrameRange(2);
skipframes = lightframeid(1)-v.FrameRange(1)-1;
sim_all = single(zeros(v.Height, v.Width, 3, Frames)/0);
trajectory = zeros(Frames, 2)/0;
for i = 1:skipframes
    step(vd.FileReader);
end
for i = 1:Frames
    if i == 1 || vd.thresholdcheck
        im = step(vd.FileReader);
        im_temp = rgb2hsv(im);
        im_temp = imfilter(im_temp, v.hfilter, 'same');
        im_v = im_temp(:, :, 3);
        im_h = im_temp(:, :, 1);
        BW_h = im_h > v.hue(1) & im_h < v.hue(2);
        if callerid == 1
            BW = v.BW1;
        elseif callerid == 2
            BW = v.BW2;
        end
        ROI = im_v(BW&BW_h);
        ROI = sort(ROI);
        BWmarker = im_v >= ROI(floor(v.threshold*numel(ROI))+1);
        BWmarker = BWmarker&BW&BW_h;
        BWmarker = imerode(BWmarker, v.se);
        BWmarker = imdilate(BWmarker, v.se);
        if get(handles.Jaw_radiobutton, 'Value')
            [row, column] = find(BWmarker == 1);
            center = [mean(column(row == max(row))) max(row)];
        elseif get(handles.LForelimb_radiobutton, 'Value')
            [area, centroid] = step(v.hblob, BWmarker);
            center = centroid(area == max(area), :);
        end
        trajectory(i, :) = center(1, :);
        simage = insertShape(im, 'FilledCircle', [center(1, :) 5], 'Color', [0 1 0], 'Opacity', 1);
        if vd.data(skipframes+i, 1)
            simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
        end
        sim_all(:, :, :, i) = simage;
    end
    if ~vd.thresholdcheck
        if i == 1
            rect = [center(1, 1)-v.width/2 center(1, 2)-v.width/2 v.width v.width];
            switch v.method
                case 1
                    points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 2
                    points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 3
                    points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 4
                    points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 5
                    points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 6
                    points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(rect));
            end
            points = points.Location;
            initialize(v.pointTracker, points, im);
            oldPoints = points;
            
            x = rect(1, 1);
            y = rect(1, 2);
            w = rect(1, 3);
            h = rect(1, 4);
            bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
        else
            % get the next frame
            im = step(vd.FileReader);
            
            % Track the points. Note that some points may be lost.
            [points, isFound] = step(v.pointTracker, im);
            visiblePoints = points(isFound, :);
            oldInliers = oldPoints(isFound, :);
            
            try % need at least 2 points
                
                % Estimate the geometric transformation between the old points
                % and the new points and eliminate outliers
                [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                    oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
                
                % Apply the transformation to the bounding box
                [bboxPolygon(1:2:end), bboxPolygon(2:2:end)] ...
                    = transformPointsForward(xform, bboxPolygon(1:2:end), bboxPolygon(2:2:end));
                
                BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(im, 1), size(im, 2));
                [~, centroid] = step(v.hblob, BW);
                trajectory(i, :) = centroid;
                simage = insertShape(im, 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
                if vd.data(skipframes+i, 1)
                    simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
                end
                sim_all(:, :, :, i) = simage;
                
                % Reset the points
                oldPoints = visiblePoints;
                setPoints(v.pointTracker, oldPoints);
            catch
                tracking_state = 0;
                break;
            end
        end
    end
end
release(v.pointTracker);
release(vd.FileReader);
trajectory_all = zeros(size(vd.data, 1), 2)/0;
trajectory_all(1:size(trajectory, 1), :) = trajectory;
trajectory_all = circshift(trajectory_all, skipframes);
trajectory = trajectory_all;
FileName = vd.FileName;
% if exist([FileName '_trajectory.mat'], 'file')
%     user_response = modaldlg;
%     switch lower(user_response)
%     case 'no'
%         % take no action
%     case 'yes'  
%         save([FileName '_trajectory.mat'], 'trajectory');
%     end
% else
%     save([FileName '_trajectory.mat'], 'trajectory');
% end
if get(handles.LForelimb_radiobutton, 'Value')
    save([FileName '_trajectory.mat'], 'trajectory');
elseif get(handles.Jaw_radiobutton, 'Value')
    save([FileName '_Jaw_trajectory.mat'], 'trajectory');
end

% if exist([FileName '_annotated.mp4'], 'file')
%     user_response = modaldlg;
%     switch lower(user_response)
%     case 'no'
%         % take no action
%         return;
%     case 'yes'  
%     end
% end
if get(handles.LForelimb_radiobutton, 'Value')
    vidObj = VideoWriter([FileName '_annotated'], 'MPEG-4');
elseif get(handles.Jaw_radiobutton, 'Value')
    vidObj = VideoWriter([FileName '_Jaw_annotated'], 'MPEG-4');
end
set(vidObj, 'FrameRate', v.FrameRate, 'Quality', 75);
open(vidObj);
imvideo = sim_all;
writeVideo(vidObj, imvideo);
close(vidObj);
clear imvideo;

function tracking_state = auto_tracking(v, vd, callerid, handles)
tracking_state = 1;
lightframeid = find(vd.data(:, 1) == 1);
Frames = v.FrameRange(1)+numel(lightframeid)+v.FrameRange(2);
skipframes = lightframeid(1)-v.FrameRange(1)-1;
sim_all = single(zeros(v.Height, v.Width, 3, Frames)/0);
trajectory = zeros(Frames, 2)/0;
for i = 1:skipframes
    step(vd.FileReader);
end
for i = 1:Frames
    if i == 1 || vd.thresholdcheck
        im = step(vd.FileReader);
        im_temp = rgb2hsv(im);
        im_temp = imfilter(im_temp, v.hfilter, 'same');
        im_v = im_temp(:, :, 3);
        im_h = im_temp(:, :, 1);
        BW_h = im_h > v.hue(1) & im_h < v.hue(2);
        if callerid == 1
            BW = v.BW1;
        elseif callerid == 2
            BW = v.BW2;
        end
        ROI = im_v(BW&BW_h);
        ROI = sort(ROI);
        BWmarker = im_v >= ROI(floor(v.threshold*numel(ROI))+1);
        BWmarker = BWmarker&BW&BW_h;
        BWmarker = imerode(BWmarker, v.se);
        BWmarker = imdilate(BWmarker, v.se);
        if get(handles.Play_checkbox, 'Value')
            set(v.him(1), 'CData', BWmarker);
            drawnow;
        end
        if get(handles.Jaw_radiobutton, 'Value')
            [row, column] = find(BWmarker == 1);
            center = [mean(column(row == max(row))) max(row)];
        elseif get(handles.LForelimb_radiobutton, 'Value')
            [area, centroid] = step(v.hblob, BWmarker);
            center = centroid(area == max(area), :);
        end
        trajectory(i, :) = center(1, :);
        simage = insertShape(im, 'FilledCircle', [center(1, :) 5], 'Color', [0 1 0], 'Opacity', 1);
        if vd.data(skipframes+i, 1)
            simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
        end
        sim_all(:, :, :, i) = simage;
        if get(handles.Play_checkbox, 'Value')
            set(v.him(2), 'CData', simage);
            drawnow;
            %             pause(0.1);
        end
    end
    if ~vd.thresholdcheck
        if i == 1
            rect = [center(1, 1)-v.width/2 center(1, 2)-v.width/2 v.width v.width];
            switch v.method
                case 1
                    points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 2
                    points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 3
                    points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 4
                    points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 5
                    points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(rect));
                case 6
                    points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(rect));
            end
            points = points.Location;
            initialize(v.pointTracker, points, im);
            oldPoints = points;
            
            x = rect(1, 1);
            y = rect(1, 2);
            w = rect(1, 3);
            h = rect(1, 4);
            bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
        else
            % get the next frame
            im = step(vd.FileReader);
            
            % Track the points. Note that some points may be lost.
            [points, isFound] = step(v.pointTracker, im);
            visiblePoints = points(isFound, :);
            oldInliers = oldPoints(isFound, :);
            
            try % need at least 2 points
                
                % Estimate the geometric transformation between the old points
                % and the new points and eliminate outliers
                [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                    oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
                
                % Apply the transformation to the bounding box
                [bboxPolygon(1:2:end), bboxPolygon(2:2:end)] ...
                    = transformPointsForward(xform, bboxPolygon(1:2:end), bboxPolygon(2:2:end));
                
                BW = poly2mask(bboxPolygon(1:2:end), bboxPolygon(2:2:end), size(im, 1), size(im, 2));
                [~, centroid] = step(v.hblob, BW);
                trajectory(i, :) = centroid;
                simage = insertShape(im, 'FilledCircle', [centroid 5], 'Color', [0 1 0], 'Opacity', 1);
                if vd.data(skipframes+i, 1)
                    simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
                end
                sim_all(:, :, :, i) = simage;
                % Display the annotated video frame using the video player object
                if get(handles.Play_checkbox, 'Value')
                    simage = insertShape(simage, 'Polygon', bboxPolygon);
                    set(v.him(2), 'CData', simage);
                    drawnow;
                end
                
%                 release(v.pointTracker);
%                 rect = [centroid(1, 1)-v.width/2 centroid(1, 2)-v.width/2 v.width v.width];
%                 switch v.method
%                     case 1
%                         points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 2
%                         points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 3
%                         points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 4
%                         points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 5
%                         points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 6
%                         points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                 end
%                 points = points.Location;
%                 initialize(v.pointTracker, points, im);
%                 oldPoints = points;
%                 
%                 x = rect(1, 1);
%                 y = rect(1, 2);
%                 w = rect(1, 3);
%                 h = rect(1, 4);
%                 bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
                % Reset the points
                oldPoints = visiblePoints;
                setPoints(v.pointTracker, oldPoints);
            catch
%                 release(v.pointTracker);
%                 im_temp = rgb2hsv(im);
%                 im_temp = imfilter(im_temp, v.hfilter, 'same');
%                 im_v = im_temp(:, :, 3);
%                 im_h = im_temp(:, :, 1);
%                 BW_h = im_h > v.hue(1) & im_h < v.hue(2);
%                 if callerid == 1
%                     BW = BW1;
%                 elseif callerid == 2
%                     BW = BW2;
%                 end
%                 ROI = im_v(BW&BW_h);
%                 ROI = sort(ROI);
%                 BWmarker = im_v >= ROI(floor(v.threshold*numel(ROI))+1);
%                 BWmarker = BWmarker&BW&BW_h;
%                 BWmarker = imerode(BWmarker, v.se);
%                 BWmarker = imdilate(BWmarker, v.se);
%                 if get(handles.Play_checkbox, 'Value')
%                     set(him(1), 'CData', BWmarker);
%                     drawnow;
%                 end
%                 if get(handles.Jaw_radiobutton, 'Value')
%                     [row, column] = find(BWmarker == 1);
%                     center = [mean(column(row == max(row))) max(row)];
%                 elseif get(handles.LForelimb_radiobutton, 'Value')
%                     [area, centroid] = step(v.hblob, BWmarker);
%                     center = centroid(area == max(area), :);
%                 end
%                 trajectory(i, :) = center(1, :);
%                 simage = insertShape(im, 'FilledCircle', [center(1, :) 5], 'Color', [0 1 0], 'Opacity', 1);
%                 if v.data(skipframes+i, 1)
%                     simage = insertShape(simage, 'FilledRectangle', [0 0 100 100], 'Color', [1 1 1], 'Opacity', 1);
%                 end
%                 sim_all(:, :, :, i) = simage;
%                 if get(handles.Play_checkbox, 'Value')
%                     set(him(2), 'CData', simage);
%                     drawnow;
%                     %             pause(0.1);
%                 end
%                 
%                 rect = [center(1, 1)-v.width/2 center(1, 2)-v.width/2 v.width v.width];
%                 switch v.method
%                     case 1
%                         points = detectMinEigenFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 2
%                         points = detectHarrisFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 3
%                         points = detectFASTFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 4
%                         points = detectBRISKFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 5
%                         points = detectSURFFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                     case 6
%                         points = detectMSERFeatures(rgb2gray(im), 'ROI', ceil(rect));
%                 end
%                 points = points.Location;
%                 initialize(v.pointTracker, points, im);
%                 oldPoints = points;
%                 
%                 x = rect(1, 1);
%                 y = rect(1, 2);
%                 w = rect(1, 3);
%                 h = rect(1, 4);
%                 bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
                tracking_state = 0;
                break;
            end
        end
    end
end
release(v.pointTracker);
release(vd.FileReader);
trajectory_all = zeros(size(vd.data, 1), 2)/0;
trajectory_all(1:size(trajectory, 1), :) = trajectory;
trajectory_all = circshift(trajectory_all, skipframes);
trajectory = trajectory_all;
FileName = vd.FileName;
% if exist([FileName '_trajectory.mat'], 'file')
%     user_response = modaldlg;
%     switch lower(user_response)
%     case 'no'
%         % take no action
%     case 'yes'  
%         save([FileName '_trajectory.mat'], 'trajectory');
%     end
% else
%     save([FileName '_trajectory.mat'], 'trajectory');
% end
if get(handles.LForelimb_radiobutton, 'Value')
    save([FileName '_trajectory.mat'], 'trajectory');
elseif get(handles.Jaw_radiobutton, 'Value')
    save([FileName '_Jaw_trajectory.mat'], 'trajectory');
end

% if exist([FileName '_annotated.mp4'], 'file')
%     user_response = modaldlg;
%     switch lower(user_response)
%     case 'no'
%         % take no action
%         return;
%     case 'yes'  
%     end
% end
if get(handles.LForelimb_radiobutton, 'Value')
    vidObj = VideoWriter([FileName '_annotated'], 'MPEG-4');
elseif get(handles.Jaw_radiobutton, 'Value')
    vidObj = VideoWriter([FileName '_Jaw_annotated'], 'MPEG-4');
end
set(vidObj, 'FrameRate', v.FrameRate, 'Quality', 75);
open(vidObj);
imvideo = sim_all;
writeVideo(vidObj, imvideo);
close(vidObj);
clear imvideo;

function CheckVideo_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
if length(ExpNumber) ~= 1
    errordlg('You can only choose one exp !', 'ERROR');
    return;
end
SpotNumber = get(handles.Spot_list, 'Value');
if length(SpotNumber) ~= 1
    errordlg('You can only choose one spot !', 'ERROR');
    return;
end
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
try
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\VRDLPParameter.mat']);
catch
    load([DataDir '\exp' num2str(ExpNumber, '%03d') '\exp' num2str(ExpNumber, '%03d') '.mat']);
    SP.Order = StimulusParameter.Order;
end
parametertext = ['Laser Knob: ' num2str(StimulusParameter.LaserKnob, '%.2f')...
    '; Pulse Width: ' num2str(StimulusParameter.PulseWidth) ' ms; Frequency: ' num2str(StimulusParameter.Frequency)...
    ' Hz; Duration: ' num2str(StimulusParameter.Duration) ' s'];
for i = 1:length(RepeatNumber)
    TrialNumber = find(SP.Order == SpotNumber);
    TrialNumber = TrialNumber(RepeatNumber(i));
    if ~get(handles.Jaw_radiobutton, 'Value')
        try
            mh1 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber) '_annotated.mp4']);
            text1 = ['exp' num2str(ExpNumber, '%03d') '\PG1camera' num2str(TrialNumber) ' | ' parametertext];
            mh1.Parent.Name = text1;
            mh1.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
        catch
            errordlg(['Cannot play exp' num2str(ExpNumber, '%03d') ': PG1 trial' num2str(TrialNumber, '%03d')], 'ERROR');
        end
    end
    try
        if get(handles.LForelimb_radiobutton, 'Value') || get(handles.RForelimb_radiobutton, 'Value')
            mh2 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) '_annotated.mp4']);
        elseif get(handles.Jaw_radiobutton, 'Value')
            mh2 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) '_Jaw_annotated.mp4']);
        end
        text2 = ['exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) ' | ' parametertext];
        mh2.Parent.Name = text2;
        mh2.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
    catch
        try
            if get(handles.Jaw_radiobutton, 'Value')
                mh2 = implay([DataDir '\exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) '_annotated.mp4']);
                text2 = ['exp' num2str(ExpNumber, '%03d') '\PG2camera' num2str(TrialNumber) ' | ' parametertext];
                mh2.Parent.Name = text2;
                mh2.Parent.Position = [412.2000  189.0000  748.0000  570.4000];
            end
        catch
            errordlg(['Cannot play exp' num2str(ExpNumber, '%03d') ': PG2 trial' num2str(TrialNumber, '%03d')], 'ERROR');
        end
    end
end
% if length(RepeatNumber) == 1
%     if RepeatNumber ~= size(get(handles.Repeat_list, 'String'), 1);
%         set(handles.Repeat_list, 'Value', RepeatNumber+1);
%     else
%         set(handles.Spot_list, 'Value', SpotNumber+1);
%         set(handles.Repeat_list, 'Value', 1);
%     end
% end

function Deletecheckvideo_button_Callback(hObject, eventdata, handles)
DataDir = get(handles.dir_text, 'String');
ExpNumber = get(handles.Exp_list, 'Value');
SpotNumber = get(handles.Spot_list, 'Value');
RepeatNumber = get(handles.Repeat_list, 'Value');
if isempty(RepeatNumber)
    errordlg('You forgot to choose repeat !', 'ERROR');
    return;
end
for i = 1:length(ExpNumber)
    try
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\VRDLPParameter.mat']);
    catch
        load([DataDir '\exp' num2str(ExpNumber(i), '%03d') '\exp' num2str(ExpNumber(i), '%03d') '.mat']);
        SP.Order = StimulusParameter.Order;
    end
    for j = 1:length(SpotNumber)
        for k = 1:length(RepeatNumber)
            TrialNumber = find(SP.Order == SpotNumber(j));
            TrialNumber = TrialNumber(RepeatNumber(k));
            if get(handles.LForelimb_radiobutton, 'Value')
                filename1 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG1camera' num2str(TrialNumber) '_annotated.mp4'];
                filename2 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_annotated.mp4'];
            elseif get(handles.Jaw_radiobutton, 'Value')
                filename2 = [DataDir '\exp' num2str(ExpNumber(i), '%03d') '\PG2camera' num2str(TrialNumber) '_Jaw_annotated.mp4'];
            end
            if ~get(handles.Jaw_radiobutton, 'Value')
                try
                    delete(filename1);
                end
            end
            try
                delete(filename2);
            end
        end
    end
end

function huelow_text_Callback(hObject, eventdata, handles)

function huelow_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function huehigh_text_Callback(hObject, eventdata, handles)

function huehigh_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Thresholding_check_Callback(hObject, eventdata, handles)

function Play_checkbox_Callback(hObject, eventdata, handles)

function Parfor_check_Callback(hObject, eventdata, handles)
