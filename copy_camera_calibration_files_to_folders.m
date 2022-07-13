%%
DataFolder = 'F:\Videos for head fixed stimulation';
CCFolder = 'C:\Users\Xu An\Desktop\Camera Calibration';
DirList = dir(DataFolder);
for n = 3:size(DirList)
    ExpName = DirList(n).name;
    if str2double(ExpName(5:6)) >= 17
        copyfile(CCFolder, [DataFolder '\' ExpName '\Camera Calibration']);
    end
end
disp('Done !');