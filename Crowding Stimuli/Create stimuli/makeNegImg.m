function [ img ] = makeNegImg( input, varargin )
% MAKENEGIMG returns the negative of an input img.
% use an image name to convert a single image or a folder name to convert
% all images in that folder. If you give a file, give it WITH extension.
% optional arguments: you can make square images and zoom in on the square.
% the first optional argument is the side of the square in pixels, the
% second is how much you want to zoom in (e.g. 2x). Use 1 if you don't want
% to zoom.

if exist(input) == 2 % input is a file
    img = rgb2gray(imread(input));
    img = imcomplement(img);
    figure(1)
    imshow(img)
    imwrite(img,[input(end-4:end), '.jpg'])

elseif exist(input) == 7 % input is a folder
    files = dir([input, '/*.jpg']);
    names = {files.name};
    
    for i = 1:length(names)
        imName = [input, '/', names{i}];
        img = rgb2gray(imread(imName));
        img = imcomplement(img);
        if varargin{1}
            img = makeSquareImg(img,varargin{1},varargin{2});
        end
        figure(1)
        imshow(img)
        drawnow
        imwrite(img,[names{i}(1:end-4), '.jpg'])
    end
    
else
    error('enter a file or folder name')
end

end

