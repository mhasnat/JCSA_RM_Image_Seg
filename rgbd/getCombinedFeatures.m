function combFeat = getCombinedFeatures(allInfoSI, scD)

if(nargin<2)
    scD = 1;
end

% Feature: Normalized L*a*b* values
labImg = rgb2lab(allInfoSI.rgbImg);
colinfosc = double(labImg(1:scD:end, 1:scD:end, :));
[r,c,d] = size(colinfosc);
featVecLab = normalizeandscale(double(reshape(colinfosc, r*c, d)));

% Feature: 3D points
% Convert to 3D points
[info3D, ~] = DepthtoCloud(allInfoSI.depImg);
info3Dsc = info3D(1:scD:end, 1:scD:end, :);
[r,c,d] = size(info3Dsc);
featVec3D = normalizeandscale(reshape(info3Dsc, r*c, d));

% Feature: image normal
normalssc = allInfoSI.imgNormals(1:scD:end, 1:scD:end, :);
[r,c,d] = size(normalssc);
featVecNormal = reshape(normalssc, r*c, d);

combFeat = [featVecLab featVec3D featVecNormal];
end