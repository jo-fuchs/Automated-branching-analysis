function Seg_Image = segment_FiberMet_v06_Function(I , thickness_pixel, gap_size_pixel, gap_bridge_check)
%% different method of bridging (connecting endpoints) + added optionality
disp('version5');
% thickness_pixel
% clc; clear all; close all; % I = (imread('0020_2.tif')); I = uint8(I); 
% I = BW;
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
[row,col,~]=size(BW);
% figure; subplot 121; imshow(B);title('After Using FiberMetric'); subplot 122; imshow(BW);  title('Initial Segmentation Result');
% imshow(BW); title('Binarized Image');

% uncon = max(size(regionprops(BW, 'Centroid')));    % How many 'un-connected' structures in the segmented image? 
% if uncon>1
%     BW = bwareafilt(BW,1);    % in case there are more than one structure in the segmented image, just keep the largest one!
% end


%% remove perfectly vertical or horizontal lines (stitching artifacts)
% find horizontal
Cleaned = BW;
toKeep = false(size(BW));
for k = 1:row
    CountConsequtive = 1;
    for j = 2:col
        if Cleaned(k,j) == 1 && Cleaned(k,j-1) == 0
            CountConsequtive = 1;
        elseif Cleaned(k,j) == 1 && Cleaned(k,j-1) == 1
            toKeep(k, j) = 1; 
            CountConsequtive = CountConsequtive + 1;
        elseif Cleaned(k,j) == 0 && Cleaned(k,j-1) == 1
            if CountConsequtive < 300
                for l = 0:CountConsequtive-1
                    toKeep(k,j-l) = 0;
                end
            end
            CountConsequtive = 0;
        else
            CountConsequtive = 0;
            
        end
    end
end

Cleaned = Cleaned & ~toKeep;


%%

toKeep = false(size(BW));

% Detect vertical
for k = 1:col
    CountConsequtive = 0;
    for j = 2:row
        if Cleaned(j,k) == 1 && Cleaned(j-1,k) == 0
            CountConsequtive = 1;
        elseif Cleaned(j,k) == 1 && Cleaned(j-1,k) == 1
            toKeep(j,k) = 1; 
            CountConsequtive = CountConsequtive + 1;
        elseif Cleaned(j,k) == 0 && Cleaned(j-1,k) == 1
            if CountConsequtive < 300
                for l = 0:CountConsequtive-1
                    toKeep(j-l,k) = 0;
                end
            end
            CountConsequtive = 0;
        else
            CountConsequtive = 0;
        end
    end
end

Cleaned = Cleaned & ~toKeep;

% Cleaned = imclearborder(Cleaned);  %remove remaining neurons connected to
% border (too dangerous, removes too many neurons)

BW = Cleaned;
BW_thin0 = bwmorph(BW , 'open', 1); % remove the remaining thin stitching artifacts



%% Bridge by finding endpoints of skeleton and connecting those below threshold


BW_thicken = bwmorph(BW_thin0 , 'thicken', round(thickness_pixel/2)); % reduce the number of points to look at
BW_thicken = bwmorph(BW_thicken , 'bridge', round(thickness_pixel/2)); % 0.5 µm

if (gap_bridge_check)
    % Bridge by connecting endpoints
    disp('briging gaps');
    disp(gap_size_pixel);
    BW_bridges = bridge_gaps(BW_thicken, gap_size_pixel);
else
    BW_bridges = BW_thicken;    
end

BW_bridges = bwmorph( BW_bridges, 'dilate', round(thickness_pixel/2)); % make connections thicker
BW_bridged = BW_thicken | BW_bridges;

BW_thin = ~bwmorph(~BW_bridged, 'open', 1); % fill tiny holes introduced by bwmorph-bridge
BW_thin = bwmorph(BW_thin , 'thin', 2);

% figure; imshow(imoverlay(BW_thin0, cumulativeLines,'green'));

%% Adaptive/Manual setting of the Size-Threshold 
Sizes = regionprops(BW_thin, 'area'); Sizes = struct2table(Sizes);
Sorted_Sizes = sortrows(Sizes); Sorted_Sizes = table2array(Sorted_Sizes);
%figure; plot(Sorted_Sizes);
P = round(Sorted_Sizes(end,:)/10);

% Manual setting of the Size-Threshold
if (P > 500)
    P = 500;
end


%% Size-Thresholding (removing smaller objects)

BW_RemovedSmallObjs = bwareaopen(BW_thin, P);  % Remove objects smaller than P pixel in size
% figure; subplot 121; imshow(BW);title('Segmentation Result'); subplot 122; imshow(BW_RemovedSmallObjs); title(['After Removing Objects Smaller than ', num2str(P), 'Pixels']);


%% Image Filling for Soma detection
% BW_filledHoles = bwmorph(BW_RemovedSmallObjs, 'fill');
% figure; imshow(BW_filledHoles); title('After Filling Holes (bwmorph)');

BW_filledHoles = imfill(BW_RemovedSmallObjs, 'holes');

% To check if Image Filling has changed the object size(s) dramatically!
% (In such cases, Im-filling will be avoided)
% Sizes = regionprops(BW_filledHoles, 'area'); Sizes = struct2table(Sizes); 
% Sorted_Sizes_new = sortrows(Sizes); Sorted_Sizes_new = table2array(Sorted_Sizes_new);
% % figure; imshow(BW_filledHoles); title('After Filling Holes (imfill)');

% find brightest parts in image & remove attached neurites
SomaInt = I > 0.4 * max(max(I));
SomaInt = bwmorph(SomaInt, 'erode', thickness_pixel);
SomaInt = bwmorph(SomaInt, 'dilate', thickness_pixel);

% find largest overlap of filled parts and high intensity areas > potential soma-seeds
Relevant_filling = SomaInt & BW_filledHoles;
Relevant_filling = bwareafilt(Relevant_filling,2); % keep only largest 2

%% To avoid over-filling in cases where dendrites branch much close around soma
% fill only parts in vicinity (here: 3 µm, could be a variable) of "intensity" soma
Fill_area = bwmorph(Relevant_filling , 'dilate', 3 * thickness_pixel);
To_fill = Fill_area & BW_RemovedSmallObjs;

Soma = imfill(To_fill, 'holes');
% -> this way overfilling is never larger than the "intensity" soma

%% merge Soma region with refined skeleton
Final_Seg = BW_RemovedSmallObjs | Soma;

%% Test if filling did increase the skeleton area at all
% original skeleton area
Sizes = regionprops(BW_RemovedSmallObjs, 'area'); Sizes = struct2table(Sizes); 
Sorted_Sizes_beforeFill = sortrows(Sizes); Sorted_Sizes_beforeFill = table2array(Sorted_Sizes_beforeFill);

% after filling skeleton area
Sizes = regionprops(Final_Seg, 'area'); Sizes = struct2table(Sizes); 
Sorted_Sizes_afterFill = sortrows(Sizes); Sorted_Sizes_afterFill = table2array(Sorted_Sizes_afterFill);


% If no Filling occured > resort to Thresholding
if Sorted_Sizes_afterFill(end,:) < 1.001 * Sorted_Sizes_beforeFill(end,:) % if the size change was not noticable
    disp('Image Filling not successful! Soma will be found by thresholding!');
    
    % find brightest parts of image and remove thin extensions
    SomaInt = I > 0.3 * max(max(I));
    SomaInt = bwmorph(SomaInt, 'erode', thickness_pixel+1);
    SomaInt = bwmorph(SomaInt, 'dilate', thickness_pixel+1);
    
    % Find the one with largest overlap to refined skeleton
    Overlap = BW_RemovedSmallObjs & SomaInt;
    SomaSeed = bwareafilt(Overlap,1);
    
    % Keep only filled parts in vicinity of Soma
    Fill_area = bwmorph(SomaSeed , 'dilate', 10 * thickness_pixel);
    Fill_area = imfill(Fill_area, 'holes');
    SomaInt = imfill(SomaInt, 'holes');
    % Fill_area = bwmorph(Fill_area, 'close', thickness_pixel);
    Soma = Fill_area & SomaInt;

    
    % Soma = bwareafilt(SomaInt,1);    % in case there are still more than one Soma, just keep the largest one!

    Final_Seg = BW_RemovedSmallObjs | Soma;

else
    disp('Soma found by filling!');
end



%% Final Segmentation Result
Seg_Image =  bwareafilt(Final_Seg,1);    % in case there are more than one object, just keep the largest one!

% Seg_Image = Final_Seg;
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