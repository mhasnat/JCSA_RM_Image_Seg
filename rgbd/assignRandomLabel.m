function imgLabels = assignRandomLabel(img)
bdry = seg2bdry(img);
% Label the image
STATS = regionprops(~bdry, 'PixelIdxList');

imgLabels = zeros(size(bdry)); % label image
cind = randperm(length(STATS));
for i=1:length(STATS)
    imgLabels(STATS(i).PixelIdxList) = cind(i)+1;
end