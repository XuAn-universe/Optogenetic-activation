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
    if strncmp(DirList(i).name(end-3:end), '.jp2', 4)
        n = n+1;
        image_name{n} = DirList(i).name;
    end
end
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(image_name)
    im = imread([DataDir '\' image_name{i}]);
    for j = 1:size(im, 3)
        if exist([DataDir '\Channel' num2str(j)]) ~= 7
            mkdir([DataDir '\Channel' num2str(j)]);
        end
        cd([DataDir '\Channel' num2str(j)]);
        imwrite(im(:,:,j), [image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none');
    end
    workbar(i/numel(image_name), [num2str(i) '/' num2str(numel(image_name))], 'Progress'); 
end
cd(curpwd);