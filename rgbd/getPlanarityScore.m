function outlierRatio = getPlanarityScore(neighbors, STATS, featVect)
% neighbors
pointPlaneTh = 0.05;
maxNumPoints = 10000;

indices1 = STATS(neighbors(1)).PixelIdxList;
indices3 = indices1;
if(length(indices1)>maxNumPoints)
    indices3 = randperm(length(indices1),maxNumPoints)';
end

indices2 = STATS(neighbors(2)).PixelIdxList;
indices4 = indices2;
if(length(indices2)>maxNumPoints)
    indices4 = randperm(length(indices2),maxNumPoints)';
end

indices = [indices1; indices2]; % index to measure the planarity
indicesT = [indices3; indices4]; % index to obtain the plane

% length(indices)

tempXYZ = featVect(indicesT, :);
[~, ~, inliers] = ransacfitplane(tempXYZ', pointPlaneTh);
outlierRatio = length(inliers) / length(indices);