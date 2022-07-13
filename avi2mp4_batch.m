%%
thr = 1;
quality = 100;
curpwd = pwd;
if ~exist('directoryname', 'var')
    directoryname = uigetdir('C:\', 'Pick a Directory');
else
    if directoryname ~= 0
        directoryname = uigetdir(directoryname, 'Pick a Directory');
    else
        directoryname = uigetdir('C:\', 'Pick a Directory');
    end
end
cd(directoryname);
DirList = dir;
% workbar(0, 'Converting Ongoing...', 'Progress');
total = length(DirList);
parfor n = 3:total
    try
        if strcmp(DirList(n).name(end-3:end), '.avi')
            readerobj = VideoReader(DirList(n).name);
            FrameRate = readerobj.FrameRate;
            vidObj = VideoWriter(DirList(n).name(1:end-4), 'MPEG-4');
            set(vidObj, 'FrameRate', FrameRate, 'Quality', quality);
            open(vidObj);
            while hasFrame(readerobj)
                vidFrame = readFrame(readerobj);
                if thr ~= 1
                    vidFrame = rgb2hsv(vidFrame);
                    intensity = vidFrame(:, :, 3);
                    intensity = mat2gray(intensity, [0 thr]);
                    vidFrame(:, :, 3) = intensity;
                    vidFrame = hsv2rgb(vidFrame);
                end
                writeVideo(vidObj, vidFrame);
            end
            close(vidObj);
            disp([DirList(n).name(1:end-4) '.mp4']);
        end
    end
    %     workbar((n-2)/(total-2), DirList(n).name, 'Progress');
end
clear readerobj;
cd(curpwd);
msgbox('Done');