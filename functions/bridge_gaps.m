function [cumulativeLines] = bridge_gaps(inputMask, gap_size_pixel)
%bridge_gaps will take a skeleton and connect endpoints closer than a minimal gap size
%   Detailed explanation goes here
% tic
BW_skel = bwmorph(inputMask, 'skel', Inf);
%figure; imshow(BW_skel);

BW_skel2 = bwareaopen(BW_skel, 12);  % Remove objects smaller than 12 pixel (3µm) in size


% find endpoints
points = bwmorph(BW_skel2, 'endpoints');
%figure; imshow(points);

endPointsInd = find(points);
[endPointRows, endPointColumns] = ind2sub(size(BW_skel), endPointsInd);
% [endPointRows, endPointColumns] = find(points);
% which objects do they belong to
[labeledImage, numberOfSegments] = bwlabel(BW_skel);

theLabels = zeros(size(endPointRows));
numberOfEndpoints = length(endPointRows);

for k = 1:numberOfEndpoints
    theLabels(k) = labeledImage(endPointsInd(k));
end


%%
% create plot to connect lines on (doesn't really work otherwise..)
image = false(size(BW_skel));
cumulativeLines = image;
tempfig = figure;
group = axes;
imshow(image);

remainingEndPoints = endPointsInd;

%for k = flip(2:numberOfEndpoints)
 while length(remainingEndPoints) > 1   
%remainingEndPoints = endPointsInd;

[thisRow, thisColumn] = ind2sub(size(BW_skel), remainingEndPoints(k));

% 
%     thisRow = endPointRows(k);
%     thisColumn = endPointColumns(k);
%     
    % Get the label number of this segment
    thisLabel = labeledImage(remainingEndPoints(k));
    
    % Get indexes of the other end points.
%     otherEndpointIndexes = setdiff(1:numberOfEndpoints, k);
    
   % remainingEndPoints = setdiff(remainingEndPoints, remainingEndPoints(k));
    
    
    
    % check if they are on the same segment
%     otherLabels = theLabels(otherEndpointIndexes);
%     onSameSegment = (otherLabels == thisLabel); % List of what segments are the same as this segment
    onSameSegment = (theLabels == thisLabel); % List of what segments are the same as this segment
    
    otherSegmentPoints = remainingEndPoints;
%     otherEndpointIndexes(onSameSegment) = []; % Remove if on the same segment
    otherSegmentPoints(onSameSegment) = [];
    
    [otherRows, otherCols] = ind2sub(size(BW_skel), otherSegmentPoints);

%     otherCols = endPointColumns(otherEndpointIndexes);
%     otherRows = endPointRows(otherEndpointIndexes);
%     
    % find nearest neighbor
    [Idx, D] = knnsearch([otherRows otherCols],[thisRow thisColumn] );
    
    % if the distance is short enough, connect endpoints
    if D < gap_size_pixel
        % Draw line from this endpoint to the other endpoint.
        linemask = imline(group,[thisColumn, thisRow ; otherCols(Idx), otherRows(Idx)]).createMask();
        cumulativeLines = cumulativeLines | linemask;

        % only connect to endpoints once?
%         remainingEndPoints(Idx) = [];
%         theLabels(Idx) = [];
%         k = k-1;
        remainingEndPoints(k) = [];
        theLabels(k) = [];
    else
         remainingEndPoints(k) = [];
         theLabels(k) = [];
    end
      k = k-1;
 %   disp(k);
end
close(tempfig);
% toc
end

