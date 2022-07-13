function varargout = Rotary_Encoder(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Rotary_Encoder_OpeningFcn, ...
                   'gui_OutputFcn',  @Rotary_Encoder_OutputFcn, ...
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

function Rotary_Encoder_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = Rotary_Encoder_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Start_Callback(hObject, eventdata, handles)
global lasttotalsteps
global lasttrig_in
global lasttimestamp
global laststeps
global x
global y
sencoder = daq.createSession('ni');
sencoder.IsContinuous = true;
ch1 = addCounterInputChannel(sencoder, 'Dev1', 0, 'Position');
ch1.EncoderType = 'X4';
counterNBits = 32;
signedThreshold = 2^(counterNBits-1);

ch2 = addAnalogInputChannel(sencoder, 'Dev1', 0, 'Voltage');
% sencoder.DurationInSeconds = Inf;
sencoder.Rate = 100;
sencoder.NotifyWhenDataAvailableExceeds = 1;
% sencoder.addTriggerConnection('External', 'Dev1/PFI0', 'StartTrigger');
% sencoder.ExternalTriggerTimeout = 10;
% sencoder.Connections(1).TriggerCondition = 'RisingEdge';

sttl = daq.createSession('ni');
addDigitalChannel(sttl, 'Dev1', 'Port0/Line0:1', 'OutputOnly');

fid = fopen('log.bin', 'wb'); % change file name in here
daxes = handles.display_axes;
cla(daxes);
showtime = 10;
x = zeros(showtime*sencoder.Rate, 1)/0;
y = zeros(showtime*sencoder.Rate, 2)/0;
hp(1) = plot(daxes, x, y(:, 1), '-k'); 
hold(daxes, 'on');
hp(2) = plot(daxes, x, y(:, 2), '-r');
box(daxes, 'off');
xlabel(daxes, 'Time (s)');
legend(daxes, 'Speed', 'Trigger In', 'Location', 'NorthEastOutside');
legend(daxes, 'boxoff');

constantv(1) = 2^counterNBits;
constantv(2) = signedThreshold;
constantv(3) = sencoder.NotifyWhenDataAvailableExceeds;
lh = addlistener(sencoder, 'DataAvailable', @(src, event)myfunc(src, event, constantv, sttl, fid, hp, daxes));

lasttotalsteps = [];
lasttrig_in = [];
lasttimestamp = [];
laststeps = [];

handles.sencoder = sencoder;
handles.sttl = sttl;
handles.lh = lh;
handles.fid = fid;
guidata(hObject, handles);

startBackground(sencoder);
% wait(sencoder, 60*60*10);

function myfunc(src, event, constantv, sttl, fid, hp, daxes)
global lasttotalsteps
global lasttrig_in
global lasttimestamp
global laststeps
global x
global y
tic
totalsteps = event.Data(:, 1);
trig_in = event.Data(:, 2);
timestamps = event.TimeStamps;
totalsteps(totalsteps > constantv(2)) = totalsteps(totalsteps > constantv(2))-constantv(1);
if isempty(lasttotalsteps)
    steps = [totalsteps(1); diff(totalsteps)];
else
    steps = [totalsteps(1)-lasttotalsteps; diff(totalsteps)];
end
data = [timestamps steps trig_in].';
fwrite(fid, data, 'double');

x(1:constantv(3), :) = [];
x = [x; timestamps];
y(1:constantv(3), :) = [];
y = [y; [steps trig_in]];
set(hp(1), 'XData', x, 'YData', y(:, 1));
set(hp(2), 'XData', x, 'YData', y(:, 2));
if ~isempty(lasttotalsteps)
    set(daxes, 'XLim', [min(x) max(x)]);
end
     
if sum(steps) > 20   %running threshold
    outputSingleScan(sttl, [1 0]);
    outputSingleScan(sttl, [0 0]);
end

if sum(steps) < -20   %running threshold
    outputSingleScan(sttl, [0 1]);
    outputSingleScan(sttl, [0 0]);
end

lasttotalsteps = totalsteps(end);
lasttrig_in = trig_in(end);
lasttimestamp = timestamps(end);
laststeps = steps(end);
toc

function Stop_Callback(hObject, eventdata, handles)
fclose(handles.fid);
delete(handles.sencoder);
delete(handles.sttl);
delete(handles.lh);
clear handles.fid handles.sencoder handles.sttl handles.lh;
disp('Experiment Over');
