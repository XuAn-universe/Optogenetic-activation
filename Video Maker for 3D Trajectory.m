%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
cycletime = 5;
framerate = 60;
[filename, pathname] = uiputfile('*.*', 'Save Video as');
if filename == 0
    return;
end
vidObj = VideoWriter([pathname filename], 'MPEG-4');
set(vidObj, 'FrameRate', framerate, 'Quality', 75);
open(vidObj);
ha = gca;
figure;
set(gcf, 'Position', get(0, 'ScreenSize'))
c = copyobj(ha, gcf);
set(c, 'Position', 'default');
axis tight;
set(gcf, 'Color', [1 1 1]);
v = get(gca, 'View');
for i = linspace(0, 360, cycletime*framerate)
    set(gca, 'View', [i, v(2)]);
    f = getframe(gcf);
    writeVideo(vidObj, f);
end
close(vidObj);
disp('Done!');