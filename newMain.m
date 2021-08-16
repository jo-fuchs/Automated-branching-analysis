%% Modifications to Segmentation & Classification


%%
clc, clear, close all

%% load raw image and pre-process
 fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_C3DIV5_C2_004.tif'; % stitching artifact example
% fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_B2B_shPTEN_024.tif'; % Soma example
% fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_B2B_shPTEN_001.tif'; % Bridge filling example

% fnameRaw = 'Raw images/Original_data_DIV7_AA6_016.tif'; % stitching artifact & soma example

% fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_B2B_ctrl_020.tif'; % Bridge filling example
% fnameRaw = 'Raw images/Original_data_DIV5_WT-KO_shRNA_B3DIV5_B_006.tif'; % Bridge filling & Soma example
% fnameRaw = 'Raw images/Original_data_DIV7_AA6_015.tif'; % Bridge filling example of unwanted filling

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
 figure; imshow(I, [0, 250]);
%% compare different segmentation functions

Seg_Image2 = segment_FiberMet_v02_Function(Smoothed_Image , thickness_pixel);         % Newer version based on FiberMetric method (v5_1 and later)
figure; imshow(Seg_Image2); title('Current Segmentation');

Seg_Image4 = segment_FiberMet_v04_Function(Smoothed_Image , thickness_pixel);         % modified soma-detection & remove stitching artifacts & bridge gaps by dilation
figure; imshow(Seg_Image4); title('Modified Soma detection with some dilation');

Seg_Image5 = segment_FiberMet_v05_Function(Smoothed_Image , thickness_pixel);         % v04 with additional gap-bridging (connecting close endpoints)
figure; imshow(Seg_Image5); title('Modified Soma detection with gap brige');

% % Differences of Gap-filled version to other versions
% figure; imshow(Seg_Image5- Seg_Image4, []);
% figure; imshow(Seg_Image5- Seg_Image2, []);




%% New Classification method on Seg_Image5

%% Find cell Body
s = 50; S = 50;   
tol = 10; extentThr = 0.25;    % extentThr lower worked slightly more consistently for me

[cBody, sklI, center] = FindCellBody(255*uint8(Seg_Image5), tol, extentThr , s , S);

allEndPoints = find(bwmorph(sklI, 'endpoints'));
% endPoints = allEndPoints;

%% Find Axon & Neurites
% tic
MIN_LEN = 40;
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

%% Measure relevant parameters

    total_neurite_length = pix_size * nnz(sklI_noSoma); 
    soma_size = pix_size * pix_size * nnz(cBody);
    
    axon_length = pix_size * nnz(axon{1});
    prim_branch_points = size(AxonBranches{2});
    prim_branch_length = pix_size * nnz(cell2mat(AxonBranches{2}));
    sec_branch_points = size(AxonBranches{3});
    sec_branch_length = pix_size * nnz(cell2mat(AxonBranches{3}));
    tert_branch_points = size(AxonBranches{4});
    tert_branch_length = pix_size * nnz(cell2mat(AxonBranches{4}));

    axon_branch_points = 0;
    axon_branch_length = 0;
    for m = 2:size(AxonBranches)
        num = size(AxonBranches{m});
        axon_branch_points = axon_branch_points + num(1);
        len = nnz(cell2mat(AxonBranches{m}));
        axon_branch_length = axon_branch_length + len;
    end
    
    
    dendrite_num = size(Neurites);
    dendrite_length = pix_size * nnz(cell2mat(Neurites));
    
    dendrite_branch_points = 0;
    dendrite_branch_length = 0;
    for m = 2:size(NeuriteBranches)
        num = size(NeuriteBranches{m});
        dendrite_branch_points = dendrite_branch_points + num(1);
        len = nnz(cell2mat(NeuriteBranches{m}));
        dendrite_branch_length = dendrite_branch_length + len;
    end

%% Visualize points
imshow(sklI);
hold on;
[r, c] = ind2sub(size(noBody), NeuriteStartPoints); plot(c, r, 'ro')