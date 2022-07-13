%%
n = 5;
DataDir = uigetdir('C:\', 'Select a folder with videos');
OutputDir = uigetdir('C:\', 'Pick a directory to rename the videos');
curpwd = pwd;
cd(DataDir);
videos = dir('*.avi');
cd(curpwd);
for i = 1:length(videos)
    vname = videos(i).name;
    movefile([DataDir '\' vname], [OutputDir '\' vname(1:9) num2str(str2double(vname(10:end-4))+n, '%d') '.avi']);
end

%%
n = 5;
DataDir = uigetdir('C:\', 'Select a folder with CSVs');
OutputDir = uigetdir('C:\', 'Pick a directory to rename the CSVs');
curpwd = pwd;
cd(DataDir);
csvs = dir('*.csv');
cd(curpwd);
for i = 1:length(csvs)
    cname = csvs(i).name;
    movefile([DataDir '\' cname], [OutputDir '\' cname(1:7) num2str(str2double(cname(8:end-4))+n, '%d') '.csv']);
end

%%
DataDir = uigetdir('C:\', 'Select a folder with tif');
curpwd = pwd;
cd(DataDir);
tifs = dir('*.tif');
for i = 1:length(tifs)
    imname = tifs(i).name;
    movefile([DataDir '\' imname], [DataDir '\' imname(1:end-4) '.tiff']);
end
cd(curpwd);

%%
DataDir = uigetdir('C:\', 'Select a folder with tiff');
curpwd = pwd;
cd(DataDir);
tiffs = dir('*.tiff');
for i = 1:length(tiffs)
    imname = tiffs(i).name;
    movefile([DataDir '\' imname], [DataDir '\' imname(1:end-5) '.tif']);
end
cd(curpwd);