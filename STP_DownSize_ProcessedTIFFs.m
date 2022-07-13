%%
downsize_factor = 0.1;
curpwd = pwd;
if exist('DataDir', 'var') ~= 1
    DataDir = uigetdir('C:\', 'Pick an image folder to process');
else
    if ~DataDir
        DataDir = uigetdir('C:\', 'Pick an image folder to process');
    else
        DataDir = uigetdir(DataDir, 'Pick an image folder to process');
    end
end
if ~DataDir
    cd(curpwd);
    return;
else
    cd(DataDir);
end
DirList = dir;
n = 0;
image_name = [];
for i = 3:size(DirList)
    if numel(DirList(i).name) < 15
        continue;
    end
    if strncmp(DirList(i).name(1:15), 'StitchedImage_Z', 15)
        n = n+1;
        image_name{n} = DirList(i).name;
    end
end

if exist('SaveDir', 'var') ~= 1
    SaveDir = uigetdir('C:\', 'Pick a folder to save image stack');
else
    if ~SaveDir
        SaveDir = uigetdir('C:\', 'Pick a folder to save image stack');
    else
        SaveDir = uigetdir(SaveDir, 'Pick a folder to save image stack');
    end
end
if ~SaveDir
    cd(curpwd);
    return;
else
    cd(SaveDir);
end
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(image_name)
    im = imread([DataDir '\' image_name{i}]);
    im = imresize(im, downsize_factor);
    if i == 1
        imwrite(im, 'Whole Brain.tiff', 'tiff', 'Compression', 'none');
    else
        imwrite(im, 'Whole Brain.tiff', 'tiff', 'Compression', 'none', 'WriteMode', 'append');
    end
    workbar(i/numel(image_name), [num2str(i) '/' num2str(numel(image_name))], 'Progress'); 
end
cd(curpwd);