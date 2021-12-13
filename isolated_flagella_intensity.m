% Detect isolated Chlamydomonas flagella and cell body and measure flagella length

% 05/30/2018
% Input image is stacked tiff file.
% Detect flagella by intensity.
% Cannot analyze crossing flagella.
% Remove flagella which connected to image boarder.
% Output flagella length and intensity as well as parameters.
% gauss3filter.m: https://jp.mathworks.com/matlabcentral/fileexchange/39013-efficient-three-dimensional-3d-gaussian-smoothing-using-convolution-via-frequency-domain


%% Read TIFF file

clear
close all

fname1 = uigetfile('*.tif');
Image1 = ReadStackedTiff2(fname1);
sizeImg1 = size(Image1);

%% parameters

radius2 = 3; % analysis radius of flagella for measuring intensity
thresh_flagella = 2.0; % thershold range for flagella (large number is more severe)
pixel_size = 0.091; % pixel size of the input image (ignore z size)

%% Show projected image

maxImg1 = max(Image1, [], 3);
figure, imshow(maxImg1, []);
% gaussImg1 = gauss3filter(Image1, [5 5 5]); 

%% Segment flagella

Image2 = gauss3filter(Image1);
mean_Image2 = mean(Image2(:));
std_Image2 = std(Image2(:));
thresh2 = mean_Image2+ std_Image2 * thresh_flagella;
level1 = graythresh(Image2);
ImageBW3 = false(sizeImg1);
ImageBW3(Image2>thresh2) = 1; % make binary images
%ImageBW3b = imclearborder(ImageBW3);
ImageBW4 = bwareaopen(ImageBW3, 500);
max_ImageBW4 = max(ImageBW4, [], 3);
sizeProImg1 = size(max_ImageBW4);
CC_BW4 = bwconncomp(ImageBW4, 26);


[LabeledImg3, ~] = bwlabeln(ImageBW4,26);
LabeledImg4 = max(LabeledImg3,[],3);
RGB2 = label2rgb(LabeledImg4, 'jet', 'k', 'shuffle');
figure, imshow(RGB2, []);
   
%% Measure flagella length and calculate total intensity

[xgrid2, ygrid2, zgrid2] = meshgrid(-radius2:radius2);
ball = (sqrt(xgrid2.^2 + ygrid2.^2 + zgrid2.^2) <= radius2);

sizeFlage1 = CC_BW4.NumObjects;
ImageData1 = cell(sizeFlage1,1);
Data2 = zeros(sizeFlage1,5);

for n07 = 1:sizeFlage1
    skelBW1 = false(sizeImg1);
    skelBW1(CC_BW4.PixelIdxList{1,n07}) = 1;
    skelBW3 = Skeleton3D(skelBW1);
    %[~,node,~] = Skel2Graph3D(skelBW3,1);
    %if length(node) == 2
        L1 = MeasureLine3D2(skelBW3);
        Data2(n07,1) = sum(L1(:,5)) * pixel_size;
        dilatedSkelBW1 = imdilate(skelBW3,ball);
        FlageRegion1 = Image1 .* dilatedSkelBW1;
        Data2(n07,2) = sum(FlageRegion1(:));
        Data2(n07,3:5) = L1(1,1:3);
        ImageData1{n07,1} = dilatedSkelBW1;
    %end
end


% dilatedSkelBW_xy = max(dilatedSkelBW1, [], 3);
% imsho(dilatedSkelBW_xy);
% figure, isosurface(dilatedSkelBW1);


%% Visualize image

ImageBW5 = zeros(sizeImg1);
for n08 = 1:sizeFlage1
    ImageBW5(ImageData1{n08,1}) = n08;
end

%figure, isosurface(ImageBW5);

LabeledImg5 = max(ImageBW5,[],3);
RGB3 = label2rgb(LabeledImg5, 'jet', 'k', 'shuffle');
figure, imshow(RGB3, []);
hold on 

for n10 = 1:sizeFlage1
    if Data2(n10,3) >0
        plot(Data2(n10,4), Data2(n10,3), 'Marker', 'o', 'MarkerFaceColor', 'w');
        text(Data2(n10,4)+6, Data2(n10,3)+6, num2str(n10), 'color','white');
    end
end
hold off


%% Extract data


writefile = fopen('flagella_intensity.txt','a+');
for n09 = 1:size(Data2, 1)
        if Data2(n09,1) >0
            fprintf(writefile, '%s	%d	%.4f	%.4f	%.2f	%d\n', fname, n09, Data2(n09,1), Data2(n09,2), thresh_flagella, radius2);
        end
end
fclose(writefile);

fname2 = strrep(fname, '.tif', '.fig');
saveas(gcf, fname2);

