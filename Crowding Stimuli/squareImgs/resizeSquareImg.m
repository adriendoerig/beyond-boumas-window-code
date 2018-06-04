% get imgs in left/right folder and make them newSize*newSize

% directory management
motherShip = fileparts(which(mfilename)); % The program directory
cd(motherShip) % go there just in case we are far away
addpath(genpath(motherShip)); % add the folder and subfolders to path
imgsPathL = [motherShip, '/left/NoZoom'];
imgsPathR = [motherShip, '/right/NoZoom'];
new_imgsPathL = [motherShip, '/left/new'];
new_imgsPathR = [motherShip, '/right/new'];


newSize = 100;


cd(imgsPathL)
files = dir('*.jpg');
names = {files.name};
for j = 1:length(names)
    cd(imgsPathL)
    img = imread(names{j});
    img = imresize(img, [newSize,newSize]);
    cd(new_imgsPathL)
    imwrite(img,[names{j}(1:end-4), '.png'])
end

cd(imgsPathR)
files = dir('*.jpg');
names = {files.name};
for j = 1:length(names)
    cd(imgsPathR)
    img = imread(names{j});
    img = imresize(img, [newSize,newSize]);
    cd(new_imgsPathR)
    imwrite(img,[names{j}(1:end-4), '.png'])
end
