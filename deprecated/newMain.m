%% Modifications to Segmentation & Classification
% Not up to date anymore


%%
clc, clear, close all

%% load raw image and pre-process
 fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_A2A_shPTEN_028.tif'; 
 pix_size = 0.2261;
 
BW = imread(fnameRaw);
[row,col,dim]=size(BW);
filt_im_ope=zeros(row,col,dim);     % final filtered data using opening morphology
conn=8;
se = strel('square',1);

for b=1:dim       % b = Index of each band of data
    mask=BW(:,:,b);    % Each band of data as the input
    marker_ope = imopen(mask, se);
    filt_im_ope(:,:,b) = imreconstruct(marker_ope,mask,conn);   % filtered band "b" , using Radius size "R"
end

Smoothed_Image = filt_im_ope;
I = Smoothed_Image; 
thickness_pixel = 4; % 1µm in these images
 figure; imshow(I, [0, 100]);
%% compare different segmentation functions

Seg_Image2 = segment_FiberMet_v02_Function(Smoothed_Image , thickness_pixel);         % Newer version based on FiberMetric method (v5_1 and later)
figure; imshow(Seg_Image2); title('Current Segmentation');

Seg_Image3 = segment_FiberMet_v03_Function(Smoothed_Image , thickness_pixel);         % modified soma-detection & remove stitching artifacts & bridge gaps by dilation
figure; imshow(Seg_Image3); title('Modified Soma detection with some dilation');

Seg_Image5 = segment_FiberMet_v05_Function(Smoothed_Image , thickness_pixel);         % v04 with additional gap-bridging (connecting close endpoints)
figure; imshow(Seg_Image5); title('Modified Soma detection with gap brige');

% % Differences of Gap-filled version to other versions
% figure; imshow(Seg_Image5- Seg_Image4, []);
% figure; imshow(Seg_Image5- Seg_Image2, []);




%% New Classification method on Seg_Image5

%% Find cell Body
s = 50; S = 50;   
tol = 10; extentThr = 0.25;    % extentThr lower worked slightly more consistently for me

[cBody, sklI, center] = FindCellBody(255*uint8(Seg_Image5), tol, extentThr , s , S, thickness_pixel);

allEndPoints = find(bwmorph(sklI, 'endpoints'));
% endPoints = allEndPoints;

%% Find Axon & Neurites
% tic
MIN_LEN = 44; % (10 µm in these images)
[Neurites, newSkel, axon, endPoints] = FindSomaNeurites(sklI, cBody, allEndPoints, MIN_LEN); 

% visualize it
image = BW;
% Axon
image = imoverlay(image, imdilate(axon{1}, ones(4)), [1 1 0] );

% Neurites
for j = 1:size(Neurites)
        image = imoverlay(image, imdilate(Neurites{j}, ones(4)), [1 0 0] );
end 

figure; imshow(image,[0, 150]);

%% Find Branches on Axon
AxonBranches = cell(10,1);
AxonBranches{1} = axon;
k = 1;
while ~isempty(AxonBranches{k}) || k > 10
  disp(['Analysing branch order ', num2str(k)]);
  k = k+1;
  [AxonBranches{k}, newSkel, endPoints] = findNextOrderBranch(newSkel, AxonBranches{k-1}, endPoints, MIN_LEN);
end
 AxonBranches = AxonBranches(~cellfun('isempty',AxonBranches));
  
disp('axon branches detected');

% visualize it 

% Branches
for k = 2:size(AxonBranches)
    %total = size(AxonBranches);
    Branches = AxonBranches{k};
    for j = 1:size(Branches)
        image = imoverlay(image, imdilate(Branches{j}, ones(4)), [0 1-mod(k-1,4)/3 1-mod(k-1,3)/4] );
    end
end


imshow(image);

%% Find branches on dendrites
NeuriteBranches = cell(10,1);
NeuriteBranches{1} = Neurites;
k = 1;
while ~isempty(NeuriteBranches{k}) || k > 10
    disp(['Analysing branch order ', num2str(k)]);  
    k = k+1;
    [NeuriteBranches{k}, newSkel, endPoints] = findNextOrderBranch(newSkel, NeuriteBranches{k-1}, endPoints, MIN_LEN);
end

disp('dendrite branches detected');

NeuriteBranches = NeuriteBranches(~cellfun('isempty',NeuriteBranches));

% visualize it 
for k = 2:size(NeuriteBranches)
    %total = size(NeuriteBranches);
    Branches = NeuriteBranches{k};
    for j = 1:size(Branches)
         image = imoverlay(image, imdilate(Branches{j}, ones(4)), [mod(k-1,4)/3 0 1-mod(k-1,3)/4] );
    end
end
imshow(image);
% toc