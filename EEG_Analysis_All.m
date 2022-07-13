function varargout = EEG_Analysis_All(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EEG_Analysis_All_OpeningFcn, ...
                   'gui_OutputFcn',  @EEG_Analysis_All_OutputFcn, ...
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

function EEG_Analysis_All_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = EEG_Analysis_All_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function dir_button_Callback(hObject, eventdata, handles)
DataDir = uigetdir('C:\', 'Select a animal folder to analysis');
if DataDir ~= 0
    handles.DataDir = DataDir;
    set(handles.dir_text, 'String', DataDir);
    RefreshLists(hObject, eventdata, handles);
end

function RefreshLists(hObject, eventdata, handles)
handles.DataDir = get(handles.dir_text, 'String');
LastExpNumber = FindLastExpNumber(handles.DataDir);
ExpNumber = str2double(LastExpNumber(8:10));
set(handles.Exp_list, 'Value', []);
set(handles.Exp_list, 'String', num2str([1:ExpNumber]'));

function dir_text_Callback(hObject, eventdata, handles)

function dir_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Close_All_button_Callback(hObject, eventdata, handles)
set(handles.figure1, 'HandleVisibility', 'off');
close all;
set(handles.figure1, 'HandleVisibility', 'on');

function RefreshLists_Callback(hObject, eventdata, handles)
RefreshLists(hObject, eventdata, handles);

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
SampleRate = 1000;
Gain = 1000;
wduration = str2double(get(handles.WDuration_text, 'String'));
wduration = round(wduration*SampleRate);
dinterval = str2double(get(handles.DInterval_text, 'String'));
dinterval = round(dinterval/1000*SampleRate);
winfun = hanning(wduration);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(wduration);          % Next highest power of 2 greater than length(x).
try
    load([DataDir '\eeg_exp' num2str(ExpNumber, '%03d') '.mat']);
    seeg = y(:, 2)/Gain;
    figure('Color', [1 1 1]);
    subplot(2, 3, 1);
    hp = plot(x, seeg, '-k');
    hold on;
    yrange = max(seeg)-min(seeg);
    ylim = [min(seeg)-0.05*yrange max(seeg)+0.05*yrange];
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
    xtime = x(1+floor(wduration/2):dinterval:length(seeg)-wduration+1+floor(wduration/2));
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
    errordlg(['Cannot analyze eeg_exp' num2str(ExpNumber, '%03d')], 'ERROR');
end

function myplot(hObject, eventdata, ha)
hf = figure('Color', [1 1 1]);
copyobj(ha, hf);
set(gca, 'Position', [0.13 0.11 0.775 0.815]);

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
        colorbar('YTick', [0.5 32.5 64.5], 'YTickLabel', {'0', '180', '360'});
    case 3
        colormap hsv;
        colorbar('YLim', [0.5 32.5], 'YTick', [0.5 16.5 32.5], 'YTickLabel', {'-90', '0', '90'});
    case 4
        colormap gray;
        hcb = colorbar;
        set(hcb, 'YTick', [0.5 32.5 64.5], 'YTickLabel', {num2str(max(UserData{3}(:))) num2str(max(UserData{3}(:))/2+min(UserData{3}(:))/2) num2str(min(UserData{3}(:)))});
    case 5
        clf;
        copyobj(UserData{3}, gcf);
        set(gca, 'Position', [0.13 0.11 0.775 0.815]);
        colormap jet;
        hcb = colorbar;
        title(hcb, UserData{4});
end
title(UserData{2});

function Exp_list_Callback(hObject, eventdata, handles)

function Exp_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
