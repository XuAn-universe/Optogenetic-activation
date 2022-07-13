%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
FrameRate = 20;
nframe = 150;
curpwd = pwd;
try
    cd(foldername);
end
foldername = uigetdir('', 'Select a folder to save');
if foldername == 0
    return;
else
    cd(foldername);
end
vidObj = VideoWriter('Supplementary Video 2 2nd', 'MPEG-4');
set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
open(vidObj);
for i = 1:nframe
    imvideo(:, :, :, i) = [PG2camera5_LeftForelimb_annotated(i).cdata PG1camera5_LeftForelimb_annotated(i).cdata];
end
writeVideo(vidObj, imvideo);
close(vidObj);
clear imvideo;
cd(curpwd)
msgbox('Done !');