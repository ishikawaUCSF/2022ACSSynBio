

% Detect Chlamydomonas flagella and cell body and measure flagella length
% and intensity.(08/23/2017)
% Output flagellar length and intensity as Data2
% Updated on 9/18/2017
%   slightly changed the method for detecting flagella
%   choose the flagella by operator for analysis and display
% Updated on 11/13/2017
%   This script works on Matlab2011a
% Need functions "ReadStackedTiff1.m", "Direction26.m", "gauss3filter.m",
% "MeasureLine3D2.m", "PeakPoint3.m"
% gauss3filter.m: https://jp.mathworks.com/matlabcentral/fileexchange/39013-efficient-three-dimensional-3d-gaussian-smoothing-using-convolution-via-frequency-domain

%% Read TIFF file

clear
close all

fname1 = uigetfile('*.tif');
Image1 = ReadStackedTiff2(fname1);
sizeImg1 = size(Image1);

%% parameters

formatSpec1 = '%s%d%d%f%f%f%f%d%d%f%d%C';
T1 = readtable('flagella_intensity1.txt','Delimiter','\t','Format',formatSpec1);

NameIdx1 = find(ismember(T1.Var1, fname1));
thresh_cell = double(table2array(T1(NameIdx1(1),6))); % thershold range for cell body (large number is more severe)
thresh_flagella = double(table2array(T1(NameIdx1(1),7))); % thershold range for flagella (large number is more severe)
radius1 = double(table2array(T1(NameIdx1(1),8))); % use for detecting flagella
radius2 = double(table2array(T1(NameIdx1(1),9))); % analysis radius of flagella for measuring intensity
matrix_correction = double(table2array(T1(NameIdx1(1),10))); % change the orientation matrix when detecting flagella
measure_length = double(table2array(T1(NameIdx1(1),11)));
TIP = double(table2array(T1(NameIdx1,3)));

pixel_size = 0.091; % pixel size of the input image (ignore z size)

%% Show projected image

maxImg1 = max(Image1, [], 3);
%figure, imshow(maxImg1, []);
gaussImg1 = gauss3filter(Image1, [7 7 5]); 

%% Detect cell bodies 

mean_gaussImg1 = mean(gaussImg1(:));
std_gaussImg1 = std(gaussImg1(:));
thresh1 = mean_gaussImg1 + std_gaussImg1 * thresh_cell;
ImageBW1 = false(sizeImg1);
ImageBW1(gaussImg1>thresh1) = 1; % make binary images
ImageBW1 = imfill(ImageBW1,'holes');
[xgrid3, ygrid3, zgrid3] = meshgrid(-3:3);
ball3 = (sqrt(xgrid3.^2 + ygrid3.^2 + zgrid3.^2) <= 3);
ImageBW1b = imerode(ImageBW1, ball3);
% CC_BW1 = bwconncomp(ImageBW1, 26);
% figure, isosurface(ImageBW1);

[LabeledImg1, ~] = bwlabeln(ImageBW1b,26);
LabeledImg2 = max(LabeledImg1,[],3);
RGB1 = label2rgb(LabeledImg2, 'jet', 'k', 'shuffle');
%figure, imshow(RGB1, []);


%% Remove cell body from image for detecting flagella

% Dilate cell body binary image
ImageBW2 = bwareaopen(ImageBW1, 1000);
[xgrid1, ygrid1, zgrid1] = meshgrid(-9:9);
structure1 = (sqrt(xgrid1.^2 + ygrid1.^2 + zgrid1.^2) <= 9);
% structure1 = strel('disk', 9);
ImageBW3 = imdilate(ImageBW2, structure1);

% Subtract cell body from filtered image
gaussImg1b = gauss3filter(Image1);
Image2 = gaussImg1b;
Image2(ImageBW3) = 0;
mean_Image2 = mean(Image2(:));
std_Image2 = std(Image2(:));
thresh2 = mean_Image2 + std_Image2 * thresh_flagella;
ImageBW4 = false(sizeImg1);
ImageBW4(Image2>thresh2) = 1; % make binary images
max_ImageBW4 = max(ImageBW4, [], 3);
sizeProImg1 = size(max_ImageBW4);

[LabeledImg3, ~] = bwlabeln(ImageBW4,26);
LabeledImg4 = max(LabeledImg3,[],3);
RGB2 = label2rgb(LabeledImg4, 'jet', 'k', 'shuffle');
%figure, imshow(RGB2, []);

%% Detect potential flagella tips

% Decide centroid of cell body in 2D
ImageXBW1 = zeros(sizeImg1 + radius1*2);
ImageXBW1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1) = ImageBW1b;
sizeImgX1 = size(ImageXBW1);
CC_XBW1 = bwconncomp(ImageXBW1, 26);
numPixels1 = cellfun(@numel,CC_XBW1.PixelIdxList);
[~,idx1] = max(numPixels1);
ImageCellBody1 = false(sizeImgX1);
ImageCellBody1(CC_XBW1.PixelIdxList{1, idx1}) = 1;
max_ImgCB1 = imfill(max(ImageCellBody1, [],3),'holes');
sizeProImgX1 = size(ImageXBW1(:,:,1));

%figure, imshow(max_ImgCB1, []);
%figure, isosurface(ImageCellBody1);


% Decide potential basal body region
RegionMaxCB1 = regionprops(max_ImgCB1, 'Area', 'centroid', 'MajorAxisLength', 'Orientation', 'Perimeter');
metric = 4*pi*RegionMaxCB1.Area/RegionMaxCB1.Perimeter^2;
Cent1 = round(RegionMaxCB1.Centroid);
if metric > 0.9
    Base1 = Cent1;
else
    angle1 = RegionMaxCB1.Orientation - 45;
    BW1 = false([1000 1000]);
    BW1(1:500,501:1000) = 1;
    BW1(501:1000,1:500) = 1;
    BW2 = imrotate(BW1, angle1);
    sizeBW2 = size(BW2);
    Cent2 = round(sizeBW2 ./2);
    x_start = Cent2(1,2)-Cent1(1,2)+1;
    y_start = Cent2(1,1)-Cent1(1,1)+1;
    BW3 = BW2(x_start:x_start+sizeProImgX1(1,1)-1, y_start:y_start+sizeProImgX1(1,2)-1);
    BW3(Cent1) = 0;
    CC_BW3 = bwconncomp(BW3);
    ImgMaxXBW4 = false(sizeProImg1 + radius1*2);
    ImgMaxXBW4(radius1+1:end-radius1, radius1+1:end-radius1) = max_ImageBW4;
    
    Probability1 = zeros(CC_BW3.NumObjects,1);
    for n00 = 1:CC_BW3.NumObjects
        %area1 = length(CC_BW3.PixelIdxList{1,n00});
        pixelNumber1 = sum(ImgMaxXBW4(CC_BW3.PixelIdxList{1,n00})); %/area1;
        Probability1(n00,1) = pixelNumber1;
    end
    [~,idx3] = max(Probability1);
    angle2 = RegionMaxCB1.Orientation * pi/180;
    x_base = round(sin(angle2) * RegionMaxCB1.MajorAxisLength * 0.45);
    y_base = round(cos(angle2) * RegionMaxCB1.MajorAxisLength * 0.45);
    base1a = [Cent1(1,1)+y_base, Cent1(1,2)-x_base];
    base2a = [Cent1(1,1)-y_base, Cent1(1,2)+x_base];
    base1b = sub2ind(sizeProImgX1, base1a(1,2), base1a(1,1));
    base2b = sub2ind(sizeProImgX1, base2a(1,2), base2a(1,1));
    if find(CC_BW3.PixelIdxList{1,idx3} == base1b)
        Base1 = base1a;
    else
        Base1 = base2a;
    end
end

% Detect potential flagella tips in 2D using the distance from cell body
ImageXBW2 = zeros(sizeImg1 + radius1*2);
ImageXBW2(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1) = ImageBW4;
CC_FLA1 = bwconncomp(ImageXBW2, 26);
TipPoint1 = zeros(CC_FLA1.NumObjects, 3);

for n01 = 1:CC_FLA1.NumObjects
    ImageFLA2 = false(sizeImgX1);
    ImageFLA2(CC_FLA1.PixelIdxList{1,n01}) = 1;
    max_ImgFLA2 = max(ImageFLA2, [],3);
    skeFla1 = bwmorph(max_ImgFLA2, 'skel',Inf);
    endpoints1 = find(bwmorph(skeFla1, 'endpoints'));
    [x1 y1] = ind2sub(sizeProImgX1, endpoints1);
    dist1 = sqrt((y1-Base1(1,1)).^2 + (x1-Base1(1,2)).^2);
    [~,idx2] = max(dist1);
    TipPoint1(n01,1:2) = [x1(idx2), y1(idx2)];   
end
    
ImageX2 = zeros(sizeImg1 + radius1*2);
ImageX2(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1) = gaussImg1b;
sizeTP1 = size(TipPoint1);
for n02 = 1:sizeTP1(1)
    z1 = PeakPoint3(ImageX2, TipPoint1(n02,1), TipPoint1(n02,2));
    TipPoint1(n02,3) = z1;
end



%% Detect flagella step by step from the base 

M1 = Direction26(radius1); %call matrix
%M2 = Direction26matrix2;
M2 = Direction26matrix3;
%M3 = M2 * 0.1;
M3 = M2 * matrix_correction;
M4 = ones(26);
M4 = M4 + M3;

%prompt = 'Input tip numbers which should be analyzed \n';
%TIP = input(prompt);

CellBodyIndex1 = CC_XBW1.PixelIdxList{1, idx1};
sizeTP2 = length(TIP);
Flagella1 = cell(sizeTP2,1);

for n03 = 1:sizeTP2
    Data1 = zeros(measure_length, 4); %change by the flagella size
    Data1(1,1:3) = TipPoint1(TIP(n03),:);
    Data1(1,5) = sub2ind(sizeImgX1, Data1(1,1), Data1(1,2), Data1(1,3));
    ImageTemp1 = ImageX2;
    for n04 = 1:length(Data1)
        Intensity1 = zeros(26,1);
        for n05 = 1:26
            %ImageTemp2 = false(sizeImgX1);
            ImageTemp2 = ImageTemp1(Data1(n04,1)-radius1:Data1(n04,1)+radius1, Data1(n04,2)-radius1:Data1(n04,2)+radius1, Data1(n04,3)-radius1:Data1(n04,3)+radius1);
            %ImageTemp2(Data1(n04,1)-radius1:Data1(n04,1)+radius1, Data1(n04,2)-radius1:Data1(n04,2)+radius1,Data1(n04,3)-radius1:Data1(n04,3)+radius1) = M1{n05,1};
            Intensity1(n05,1) = sum(ImageTemp2(M1{n05,1}>0))/M1{n05,3};
        end
        
        ImageTemp3 = false(sizeImgX1);
        
        if n04 > 1
            Intensity2 = Intensity1 .* M4(:,Data1(n04,4));
            [Max1,Index1] = max(Intensity2);
            Sub2 = Data1(n04,1:3) + M1{Index1,2};
            Ind2 = sub2ind(sizeImgX1, Sub2(1,1), Sub2(1,2), Sub2(1,3));
            if isempty(find(CellBodyIndex1 == Ind2, 1)) && isempty(find(Data1(:,5) == Ind2, 1))
                    if Max1 > ImageTemp1(Data1(n04,1),Data1(n04,2),Data1(n04,3))/5
                        Data1(n04+1,1:3) = Data1(n04,1:3) + M1{Index1,2};
                        Data1(n04+1,4) = Index1;
                        Data1(n04+1,5) = Ind2;
                        ImageTemp3(Data1(n04,1)-radius1:Data1(n04,1)+radius1, Data1(n04,2)-radius1:Data1(n04,2)+radius1,Data1(n04,3)-radius1:Data1(n04,3)+radius1) = M1{M1{Index1,4},1};
                        ImageTemp1(ImageTemp3) = 0;
                        ImageTemp1(Data1(n04,5)) = 0;
                    else
                        comment = ['flagellum ', num2str(n03), ' might have accidentally stopped'];
                        disp(comment);
                        break
                    end
            elseif find(CellBodyIndex1 == Ind2, 1) > 0
                comment = ['flagellum ', num2str(n03), ' reached to the cell body.'];
                disp(comment);
                break
            elseif find(Data1(:,5) == Ind2, 1) > 0
                comment = ['flagellum ', num2str(n03), ' encounted to own flagellum.'];
                disp(comment);
                break
            end
        else
            [Max1,Index1] = max(Intensity1);
            Data1(n04+1,1:3) = Data1(n04,1:3) + M1{Index1,2};
            Data1(n04+1,4) = Index1;
            Data1(n04+1,5) = sub2ind(sizeImgX1, Data1(n04+1,1), Data1(n04+1,2), Data1(n04+1,3));
            ImageTemp3(Data1(n04,1)-radius1:Data1(n04,1)+radius1, Data1(n04,2)-radius1:Data1(n04,2)+radius1,Data1(n04,3)-radius1:Data1(n04,3)+radius1) = M1{M1{Index1,4},1};
            ImageTemp1(ImageTemp3) = 0;
            ImageTemp1(Data1(n04,5)) = 0;
        end
    end
    W1 = Data1(:,5) > 0;
    Flagella1{n03,1} = Data1(W1,5);
    %Flagella1{n03,2} = comment;
end


%% Measure flagella length and calculate total intensity

[xgrid2, ygrid2, zgrid2] = meshgrid(-radius2:radius2);
ball = (sqrt(xgrid2.^2 + ygrid2.^2 + zgrid2.^2) <= radius2);

sizeFlage1 = length(Flagella1);
ImageData1 = cell(sizeFlage1,1);
Data2 = zeros(sizeFlage1,2);
Flagella2 = cell(sizeTP2,1);

for n06 = 1:sizeFlage1
    skelBW1 = false(sizeImgX1);
    skelBW1(Flagella1{n06,1}) = 1;
    %skelBW2 = skelBW1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1);
    skelBW3 = Skeleton3D(skelBW1);
    [~,node,~] = Skel2Graph3D(skelBW3,1);

    L1 = MeasureLine3D2(skelBW3);
    Data2(n06,1) = sum(L1(:,5)) * pixel_size; %Flagellar legth (um)
    dilatedSkelBW1 = imdilate(skelBW3,ball);
    dilatedSkelBW2 = dilatedSkelBW1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1);
    FlageRegion1 = Image1 .* dilatedSkelBW2;
    Data2(n06,2) = sum(FlageRegion1(:)); %Total intensity of the flagellum
    ImageData1{n06,1} = dilatedSkelBW2;
    lengthL1 = size(L1, 1);  
    if length(node) == 2
        for n10 = 1:lengthL1
            skelBW4 = false(sizeImgX1);
            skelBW4(L1(n10,1),L1(n10,2),L1(n10,3)) = 1;
            dilatedPoint1 = imdilate(skelBW4, ball);
            dilatedPoint2 = dilatedPoint1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1);
            flagePoint1 = Image1 .* dilatedPoint2;
            L1(n10,6) = sum(flagePoint1(:));

        end
Flagella2{n06,1} = L1(:,6);
figure, plot(L1(:,6));
    end
end



dilatedSkelBW_xy = max(dilatedSkelBW1, [], 3);
imsho(dilatedSkelBW_xy);
figure, isosurface(dilatedSkelBW1);


%% Visualize image

ImageBW5 = ImageCellBody1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1);
ImageBW6 = uint8(ImageCellBody1(radius1+1:end-radius1, radius1+1:end-radius1, radius1+1:end-radius1));
for n07 = 1:sizeFlage1
    ImageBW5(ImageData1{n07,1}) = 1;
    ImageBW6(ImageData1{n07,1}) = n07+1;
end
%figure, isosurface(ImageBW5);

LabeledImg5 = max(ImageBW6,[],3);
RGB3 = label2rgb(LabeledImg5, 'jet', 'k', 'shuffle');
figure, imshow(RGB3, []);
hold on 

plot(Base1(1,1)-radius1, Base1(1,2)-radius1, 'Marker', 'o', 'MarkerFaceColor', 'red');
for n08 = 1:sizeFlage1
    plot(TipPoint1(TIP(n08),2)-radius1, TipPoint1(TIP(n08),1)-radius1, 'Marker', 'o', 'MarkerFaceColor', 'w');
    text(TipPoint1(TIP(n08),2)-radius1+4, TipPoint1(TIP(n08),1)-radius1+4, num2str(TIP(n08)), 'color','white');
end
hold off

%% Extract data

%{
writefile = fopen('flagella_intensity.txt','a+');
for n09 = 1:size(Data2, 1)
    fprintf(writefile, '%s	%d	%d	%.4f	%.4f	%.2f	%.2f	%d	%d	%.3f	%d\n', fname1, n09, TIP(n09), Data2(n09,1), Data2(n09,2), thresh_cell, thresh_flagella, radius1, radius2, matrix_correction, measure_length);
end
fclose(writefile);

fname12 = strrep(fname1, '.tif', '.fig');
saveas(gcf, fname12);
%}