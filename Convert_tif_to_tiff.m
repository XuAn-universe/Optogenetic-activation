%%
curpwd = pwd;
if exist('DataDir', 'var') ~= 1
    DataDir = uigetdir('C:\', 'Select an image folder to process');
else
    DataDir = uigetdir(DataDir, 'Select an image folder to process');
end
if isempty(DataDir)
    cd(curpwd);
    return;
else
    cd(DataDir);
end
DirList = dir;
n = 0;
image_name = [];
for i = 3:size(DirList)
    if strncmp(DirList(i).name(end-3:end), '.tif', 4)
        n = n+1;
        image_name{n} = DirList(i).name;
    end
end
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(image_name)
    im = imread([DataDir '\' image_name{i}]);
    if exist([DataDir '\TIFF']) ~= 7
        mkdir([DataDir '\TIFF']);
    end
    imwrite(im, [DataDir '\TIFF\' image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none');
    workbar(i/numel(image_name), [num2str(i) '/' num2str(numel(image_name))], 'Progress'); 
end
cd(curpwd);