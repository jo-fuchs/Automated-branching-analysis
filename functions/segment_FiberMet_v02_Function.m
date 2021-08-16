function Seg_Image = segment_FiberMet_v02_Function(I , thickness_pixel)
% thickness_pixel
disp('Current version');
% clc; clear all; close all; % I = (imread('0020_2.tif')); I = uint8(I); 
if length(size(I)) == 3
    I = rgb2gray(I);
end
%figure; imshow(I,[]);  title ('Raw Image');

%imagesc(I);
%B = fibermetric(I,7,'ObjectPolarity','dark');
 tic

%% Fiber Metric
% B = fibermetric(I);
B = fibermetric(I,'StructureSensitivity',thickness_pixel);


%figure;imshow(B); title('After Using FiberMetric');
% figure; subplot 121; imshow(I,[]); title('Raw Image'); subplot 122; imshow(B);  title('After Using FiberMetric');

%% Initial Segmentation
BW = B > 0.05;
% figure; subplot 121; imshow(B);title('After Using FiberMetric'); subplot 122; imshow(BW);  title('Initial Segmentation Result');
% imshow(BW); title('Binarized Image');

% uncon = max(size(regionprops(BW, 'Centroid')));    % How many 'un-connected' structures in the segmented image? 
% if uncon>1
%     BW = bwareafilt(BW,1);    % in case there are more than one structure in the segmented image, just keep the largest one!
% end

%% Adaptive/Manual setting of the Size-Threshold 
Sizes = regionprops(BW, 'area'); Sizes = struct2table(Sizes);
Sorted_Sizes = sortrows(Sizes); Sorted_Sizes = table2array(Sorted_Sizes);
%figure; plot(Sorted_Sizes);
P = round(Sorted_Sizes(end,:)/10);

% Manual setting of the Size-Threshold
% P = 500;

%% Size-Thresholding (removing smaller objects)

BW_RemovedSmallObjs = bwareaopen(BW, P);  % Remove objects smaller than P pixel in size
% figure; subplot 121; imshow(BW);title('Segmentation Result'); subplot 122; imshow(BW_RemovedSmallObjs); title(['After Removing Objects Smaller than ', num2str(P), 'Pixels']);

%% Image Filling
% BW_filledHoles = bwmorph(BW_RemovedSmallObjs, 'fill');
% figure; imshow(BW_filledHoles); title('After Filling Holes (bwmorph)');

BW_filledHoles = imfill(BW_RemovedSmallObjs, 'holes');

% To check if Image Filling has changed the object size(s) dramatically!
% (In such cases, Im-filling will be avoided)
Sizes = regionprops(BW_filledHoles, 'area'); Sizes = struct2table(Sizes); 
Sorted_Sizes_new = sortrows(Sizes); Sorted_Sizes_new = table2array(Sorted_Sizes_new);
% figure; imshow(BW_filledHoles); title('After Filling Holes (imfill)');

%% Choose between imfill OR Soma-Approximation
if Sorted_Sizes_new(end,:) > 1.2 * Sorted_Sizes(end,:)
    display('Image Filling is avoided! Some will be found instead!');
    
    Soma = I > 0.4 * max(max(I)); 
%     Soma = bwareafilt(Soma,1);    % in case there are more than one Soma, just keep the largest one!
%     figure; imshow(Soma);  title ('Soma');
    Final_seg_OR = BW_RemovedSmallObjs | Soma;

else
    display('Image Filling is done!');
    Final_seg_OR = BW_filledHoles;
%     figure; subplot 121; imshow(BW_RemovedSmallObjs); subplot 122; imshow(BW_filledHoles); title('After Filling Holes (imfill)');
end



%% Final Segmentation Result
Final_seg =  bwareafilt(Final_seg_OR,1);    % in case there are more than one object, just keep the largest one!
Seg_Image = Final_seg;
 toc

% figure; imshow(Final_seg); title('Final Segmentation');
% figure; subplot 121; imshow(I , []); title('Raw Image'); subplot 122; imshow(Seg_Image); title('Final Segmentation');

% figure; subplot 121; imshow(I , []); title('Raw Image'); subplot 122; imshow(BW_filledHoles); title('Final Segmentation');
% figure; subplot 121; imshow(I , []); title('Raw Image'); subplot 122; imshow(BW_open); title('Final Segmentation (After Thickening)');

%% Thickening
% BW_thicken = bwmorph(BW_filledHoles , 'thicken');
% figure; subplot 121; imshow(BW_filledHoles); title('Before Thickening'); 
% subplot 122; imshow(BW_thicken); title('After Thickening');

%% Bridging
% BW_bridge = bwmorph(BW_thicken , 'bridge');
% figure; subplot 121; imshow(BW_thicken); title('Before Bridging'); 
% subplot 122; imshow(BW_bridge); title('After Bridging');
% 
% figure; subplot 121; imshow(BW_filledHoles); title('Before Thickening/Bridging'); 
% subplot 122; imshow(BW_bridge); title('After Thickening/Bridging');