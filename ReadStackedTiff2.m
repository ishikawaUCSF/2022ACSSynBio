function images = ReadStackedTiff2(filename)

% 07/18/2017
% Read stacked tiff file in the folder
% Input filename and get image
% Updated on 12/13/2021
% 


info = imfinfo(filename);
num_images = numel(info);
imStack1 = zeros(info(1,1).Height, info(1,1).Width, num_images);
for n01 = 1:num_images
    imStack1(:,:,n01) = im2double(imread(filename, n01, 'Info', info));
end

images = imStack1;