function varargout = AddFrame_random(varargin)
%*------------------------------------------------*
% Developmental Plasticity of Sensory System Lab
% University of Science and Technology of China
% Author : Xu An, Aug. 2011
% anxu@mail.ustc.edu.cn
%*------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AddFrame_random_OpeningFcn, ...
                   'gui_OutputFcn',  @AddFrame_random_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function AddFrame_random_OpeningFcn(hObject, eventdata, handles, varargin)
global TagDir
set(handles.DataDir,'String',TagDir)
RefreshExp(hObject, eventdata, handles)
handles.output = hObject;
guidata(hObject, handles);
function varargout = AddFrame_random_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
function figure1_CreateFcn(hObject, eventdata, handles)
function ExpText_Callback(hObject, eventdata, handles)
function ExpText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DataDir_Callback(hObject, eventdata, handles)	
function DataDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ConA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ConB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Check1_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check2_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check3_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check4_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check5_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check6_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)
function Check7_Callback(hObject, eventdata, handles)
    RefreshImage(hObject, eventdata, handles)

function DirButton_Callback(hObject, eventdata, handles)
	DataDir = uigetdir('','');
	if DataDir ~= 0
		set(handles.DataDir,'String',DataDir);
	end	
	RefreshExp(hObject, eventdata, handles)
function RefreshExp(hObject, eventdata, handles)
	try
        set(handles.ConA,'Value',[]);
        set(handles.ConB,'Value',[]);
        set(handles.ConA_rand,'Value',[]);
        set(handles.ConB_rand,'Value',[]);
		LastExpFolder = FindLastExpFolder(get(handles.DataDir,'String'));
		ExpNumber = str2double(LastExpFolder(4:6));
		set(handles.ConA, 'String',num2str([1:ExpNumber]'))
		set(handles.ConB, 'String',num2str([1:ExpNumber]'))
        set(handles.ConA_rand, 'String',num2str([1:ExpNumber]'))
		set(handles.ConB_rand, 'String',num2str([1:ExpNumber]'))
	end
		
	
function PlotButton_Callback(hObject, eventdata, handles)
    global stack_diff
    global NumData
    [stack_diff NumData] = AddFunction(hObject, eventdata, handles);
    colormap(gray)
    for i = 1:7
        eval(['axes (handles.Frame' num2str(i) ')']);
        imagesc (reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3)), ...
            [mean2(reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3)))-3*std2(reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3)))...
            mean2(reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3)))+3*std2(reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3)))])
        title(['Frame' num2str(i)])
        axis off
        axis image
        colorbar
    end
    RefreshImage(hObject, eventdata, handles)

    
function RefreshImage(hObject, eventdata, handles)
    global stack_diff
    global NumData
    global FrameImage
    global FrameImage_show
    FrameImage = zeros(NumData.ROI(4),NumData.ROI(3));
    for i = 1:7
        FrameNum(i) = eval(['get(handles.Check' num2str(i) ',''Value'')']);        
        FrameImage = FrameImage + FrameNum(i)*reshape(stack_diff(:,i),NumData.ROI(4),NumData.ROI(3));
    end    
    try
        FrameImage = FrameImage./(sum(FrameNum(:)));
        FrameImage_show = FrameImage;
        if get(handles.filter_text,'Value')==1
            FrameImage_show = bandpass_disk_filter(FrameImage,3);
        end
    end
    axes (handles.Image);
    try
        h=imagesc (FrameImage_show, [mean2(FrameImage_show(:))-3*std2(FrameImage_show(:)) mean2(FrameImage_show(:))+3*std2(FrameImage_show(:))]);
    catch
        h=imagesc (FrameImage_show);
    end
    axis off
    axis image
    colorbar
	set(h, 'ButtonDownFcn', @Image_ButtonDownFcn)
	
	
function SaveImage_Callback(hObject, eventdata, handles)
    global FrameImage
    DataDir = get(handles.DataDir, 'String');
    curpwd = pwd;
    cd(DataDir)
    [file,path] = uiputfile('*.mat','Save Frame');
    if file ~= 0
		save([path file], 'FrameImage');
	end
    cd(curpwd)


function AtoB_Callback(hObject, eventdata, handles)
	if get(handles.AtoB,'Value')
		AValue = get(handles.ConA, 'Value');
		set(handles.ConB, 'Value', AValue);
	end
	
function ConA_Callback(hObject, eventdata, handles)
	AtoB_Callback(hObject, eventdata, handles)

function ConB_Callback(hObject, eventdata, handles)
	AtoB_Callback(hObject, eventdata, handles)

function ESDButton_Callback(hObject, eventdata, handles)
    [stack_diff NumData] = AddFunction(hObject, eventdata, handles);
	[ind_comp,W,D,M1,M2] = ESD_Core(stack_diff',NumData.ROI(4),NumData.ROI(3),5,5,0);
    if get(handles.AtoB,'Value')||get(handles.randAtorandB,'Value')
        TitleText = ['Exp ' num2str(get(handles.ConA,'Value')) '  ' num2str(get(handles.ConA_rand,'Value'))];
    else
        TitleText = ['ConA:  ' num2str(get(handles.ConA,'Value')) '  ' num2str(get(handles.ConA_rand,'Value')) '     '  'ConB:  ' num2str(get(handles.ConB,'Value')) '  ' num2str(get(handles.ConB_rand,'Value'))];
    end
    ESDView(NumData,TitleText,stack_diff,ind_comp);
	
function Image_ButtonDownFcn(hObject, eventdata, handles)
	global FrameImage_show
	figure
	try
        imagesc (FrameImage_show, [mean2(FrameImage_show(:))-3*std2(FrameImage_show(:)) mean2(FrameImage_show(:))+3*std2(FrameImage_show(:))]);
    catch
        imagesc (FrameImage_show)
    end
    axis off
    axis image
    colorbar
	colormap(gray)
    uicontrol('Style', 'PushButton', 'String', 'Save', 'Position', [100,50,100,30],...
			'CallBack', @SaveImage);

function AddMat_Callback(hObject, eventdata, handles)
    AddMat;


function [stack_diff NumData]=AddFunction(hObject, eventdata, handles)
    DataDir = get(handles.DataDir,'String');
	ConAText = get(handles.ConA,'Value');
	ConBText = get(handles.ConB,'Value');
    ConA_randText = get(handles.ConA_rand,'Value');
	ConB_randText = get(handles.ConB_rand,'Value');
	stack1 = 0;
	stack2 = 0;
    stack3 = 0;
    stack4 = 0;
    stack_diff1 = 0;
    stack_diff2 = 0;
	stack_diff = 0;
    if ConAText~=0
	    load ([DataDir '\esd\exp' num2str(ConAText(1),'%03u') '\NumData.mat'])
    else
        load ([DataDir '\esd\exp' num2str(ConA_randText(1),'%03u') '\NumData.mat'])
    end
    if get(handles.minus,'value')
	    for i = 1:length(ConAText)
		    if get(handles.ConAodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack1;
            elseif get(handles.ConAeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack1;
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack2;
            end
        end
        for i = 1:length(ConA_randText)
		    if get(handles.ConA_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack1;
            elseif get(handles.ConA_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack2;
            elseif get(handles.ConA_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack3.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack3;
            elseif get(handles.ConA_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack4.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack1;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack2;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack3.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack3;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\ffa_stack4.mat'])
			    stack_diff1 = stack_diff1 + ffa_stack4;
            end
        end
	    for i = 1:length(ConBText)
		    if get(handles.ConBodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack1;
		    elseif get(handles.ConBeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack1;
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack2;
            end
        end
        for i = 1:length(ConB_randText)
		    if get(handles.ConB_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack1;
            elseif get(handles.ConB_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack2;
            elseif get(handles.ConB_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack3.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack3;
            elseif get(handles.ConB_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack4.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack1.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack1;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack2.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack2;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack3.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack3;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\ffa_stack4.mat'])
			    stack_diff2 = stack_diff2 - ffa_stack4;
            end
        end
    elseif get(handles.amplitude,'value')
        for i = 1:length(ConAText)
		    if get(handles.ConAodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff1 = stack_diff1 + norm_stack1;
            elseif get(handles.ConAeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff1 = stack_diff1 + norm_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff1 = stack_diff1 + norm_stack1;
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff1 = stack_diff1 + norm_stack2;
            end
        end
        for i = 1:length(ConA_randText)
		    if get(handles.ConA_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff1 = stack_diff1 + norm_stack1;
            elseif get(handles.ConA_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff1 = stack_diff1 + norm_stack2;
            elseif get(handles.ConA_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack3.mat'])
			    stack_diff1 = stack_diff1 + norm_stack3;
            elseif get(handles.ConA_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack4.mat'])
			    stack_diff1 = stack_diff1 + norm_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff1 = stack_diff1 + norm_stack1;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff1 = stack_diff1 + norm_stack2;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack3.mat'])
			    stack_diff1 = stack_diff1 + norm_stack3;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\norm_stack4.mat'])
			    stack_diff1 = stack_diff1 + norm_stack4;
            end
        end
	    for i = 1:length(ConBText)
		    if get(handles.ConBodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff2 = stack_diff2 - norm_stack1;
		    elseif get(handles.ConBeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff2 = stack_diff2 - norm_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff2 = stack_diff2 - norm_stack1;
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff2 = stack_diff2 - norm_stack2;
            end
        end
        for i = 1:length(ConB_randText)
		    if get(handles.ConB_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff2 = stack_diff2 - norm_stack1;
            elseif get(handles.ConB_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff2 = stack_diff2 - norm_stack2;
            elseif get(handles.ConB_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack3.mat'])
			    stack_diff2 = stack_diff2 - norm_stack3;
            elseif get(handles.ConB_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack4.mat'])
			    stack_diff2 = stack_diff2 - norm_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack1.mat'])
			    stack_diff2 = stack_diff2 - norm_stack1;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack2.mat'])
			    stack_diff2 = stack_diff2 - norm_stack2;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack3.mat'])
			    stack_diff2 = stack_diff2 - norm_stack3;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\norm_stack4.mat'])
			    stack_diff2 = stack_diff2 - norm_stack4;
            end
        end
    else
        for i = 1:length(ConAText)
		    if get(handles.ConAodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff1 = stack_diff1 + divide_stack1;
            elseif get(handles.ConAeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff1 = stack_diff1 + divide_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff1 = stack_diff1 + divide_stack1;
                load ([DataDir '\esd\exp' num2str(ConAText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff1 = stack_diff1 + divide_stack2;
            end
        end
        for i = 1:length(ConA_randText)
		    if get(handles.ConA_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff1 = stack_diff1 + divide_stack1;
            elseif get(handles.ConA_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff1 = stack_diff1 + divide_stack2;
            elseif get(handles.ConA_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack3.mat'])
			    stack_diff1 = stack_diff1 + divide_stack3;
            elseif get(handles.ConA_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack4.mat'])
			    stack_diff1 = stack_diff1 + divide_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff1 = stack_diff1 + divide_stack1;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff1 = stack_diff1 + divide_stack2;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack3.mat'])
			    stack_diff1 = stack_diff1 + divide_stack3;
                load ([DataDir '\esd\exp' num2str(ConA_randText(i),'%03u') '\divide_stack4.mat'])
			    stack_diff1 = stack_diff1 + divide_stack4;
            end
        end
	    for i = 1:length(ConBText)
		    if get(handles.ConBodd,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff2 = stack_diff2 - divide_stack1;
		    elseif get(handles.ConBeven,'Value')
			    load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff2 = stack_diff2 - divide_stack2;
            else
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff2 = stack_diff2 - divide_stack1;
                load ([DataDir '\esd\exp' num2str(ConBText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff2 = stack_diff2 - divide_stack2;
            end
        end
        for i = 1:length(ConB_randText)
		    if get(handles.ConB_rand1,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff2 = stack_diff2 - divide_stack1;
            elseif get(handles.ConB_rand2,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff2 = stack_diff2 - divide_stack2;
            elseif get(handles.ConB_rand3,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack3.mat'])
			    stack_diff2 = stack_diff2 - divide_stack3;
            elseif get(handles.ConB_rand4,'Value')
			    load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack4.mat'])
			    stack_diff2 = stack_diff2 - divide_stack4;
            else
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack1.mat'])
			    stack_diff2 = stack_diff2 - divide_stack1;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack2.mat'])
			    stack_diff2 = stack_diff2 - divide_stack2;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack3.mat'])
			    stack_diff2 = stack_diff2 - divide_stack3;
                load ([DataDir '\esd\exp' num2str(ConB_randText(i),'%03u') '\divide_stack4.mat'])
			    stack_diff2 = stack_diff2 - divide_stack4;
            end
        end
    end
    if get(handles.ConAoddeven,'Value')||get(handles.ConA_all,'Value')
        stack_diff1=stack_diff1./(2*length(ConAText)+4*length(ConA_randText));
    else
        stack_diff1=stack_diff1./(length(ConAText)+length(ConA_randText));
    end
    if isempty(ConAText) && isempty(ConA_randText)
        stack_diff1 = 0;
    end
    if get(handles.ConBoddeven,'Value')||get(handles.ConB_all,'Value')
        stack_diff2=stack_diff2./(2*length(ConBText)+4*length(ConB_randText));
    else
        stack_diff2=stack_diff2./(length(ConBText)+length(ConB_randText));
    end
    if isempty(ConBText) && isempty(ConB_randText)
        stack_diff2 = 0;
    end
    stack_diff=stack_diff1+stack_diff2;
    
    
function SaveImage(hObject, eventdata, handles)
	global FrameImage
	[file,path] = uiputfile('*.mat','Save Frame');
	if file ~= 0
        save([path file], 'FrameImage');
    end
function ConA_rand_Callback(hObject, eventdata, handles)
randAtorandB_Callback(hObject, eventdata, handles)
function ConA_rand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ConB_rand_Callback(hObject, eventdata, handles)
randAtorandB_Callback(hObject, eventdata, handles)
function ConB_rand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function randAtorandB_Callback(hObject, eventdata, handles)
    if get(handles.randAtorandB,'Value')
        AValue = get(handles.ConA_rand, 'Value');
		set(handles.ConB_rand, 'Value', AValue);
    end
function filter_text_Callback(hObject, eventdata, handles)
