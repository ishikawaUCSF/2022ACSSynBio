function z = PeakPoint3(image, x, y)

% detect a peak point (08/18/2017)
% input image should be 3D image


sizeImage = size(image);

ImageX1 = false(sizeImage(1), sizeImage(2));
ImageX1(x-2:x+2, y-2:y+2) = 1;
z1 = sizeImage(3);
Index1 = zeros(25, z1);

for n1 = 1:z1
    imageXY = image(:,:,n1);
    Index1(:,n1) = imageXY(ImageX1);
end

Index2 = sum(Index1, 1);
dataD1 = smooth(Index2, 3, 'sgolay', 1); 
PeakList1 = mspeaks(1:z1, dataD1, 'denoising', false);
[~, maxInd1] = max(PeakList1(:,2), [], 1);
z2 = round(PeakList1(maxInd1, 1));

z = z2;
