function [cBody, skel, center] = FindCellBody(I, tol, extentThr , s , S, thickness_pixel)

% sklI = bwmorph(I, 'skel', Inf);
skel = bwmorph(I,'thin',inf);

while true
    B = imerode(I~=0, ones(s));
    B = bwareafilt(B,1);    % in case there are more than one Soma, just keep the largest one!
    if nnz(B) < S
        s = s - 2;
    else
        rp = regionprops(B, 'Centroid');
        y = round(rp.Centroid(1,1));
        x = round(rp.Centroid(1,2));
        break
    end
end

if~exist('extentThr', 'var')
    extentThr = -1;
end

cBody = false(size(I));
cBodyOld = cBody;
cBody(x, y) = true;

while(sum(cBody(:)) ~= sum(cBodyOld(:)))
    cBodyOld = cBody;
    p = imdilate(cBody, strel('disk',thickness_pixel, 0)) - cBody;
    vp = find(p);
    vpVal = I(vp);
    meanSeg = mean(I(cBody));
    cBody(vp(vpVal > meanSeg - tol & vpVal < meanSeg + tol)) = true;
    cBody = bwareafilt(cBody,1);    % in case there are more than one Soma, just keep the largest one!
    rp = regionprops(cBody, 'Extent');
    if rp.Extent < extentThr
        break
    end
end

cBody = imreconstruct(cBody,imerode(cBody, strel('disk', thickness_pixel, 0))); % erode by µm distance instead of pixel value (0.5µm)

I = cBody & skel;
% figure; imshow(I);
% figure; imshow(cBody);  title('Soma is found!'); pause(0.5);
cBody = bwareafilt(cBody,1);    % Again, in case there are more than one Soma, just keep the largest one!
% figure; imshow(cBody);
rp = regionprops(cBody, 'Centroid');

% rp.Centroid
c = round(rp.Centroid);
[~, ind] = bwdist(I);
% [x, y] = ind2sub(size(I), ind(c(2), c(1)));
% center = [x, y];
center = double(ind(c(2), c(1)));
