clear all; close all; clc;

MethodType = 'rb_jcsa_rm'; % 'rb_jcsa'

% Load image and other information
fileName = 'rgbd_info_1.mat'; % this file contains the color image, depth image and image normals
load(fileName);

rgbImg = rgbd_data.rgbImg;
depImg = rgbd_data.depImg;
imgNormals = rgbd_data.imgNormals;
allInfo = rgbd_data;

% Display images
subplot(2, 3, 1); imshow(rgbImg); title('Color Image');
subplot(2, 3, 2); imshow(depImg, []); title('Depth Image');
subplot(2, 3, 3); imshow(imgNormals); title('Image Normals');

% Load options and the threshold values
sc = 2;
kMax = 20;

thOptions.thDivNormalMax = 2;
thOptions.thDivNormalMin = 1;
thOptions.planarityTh = 0.9;
thOptions.thKappa = 5;
thOptions.edgeStrengthTh = 0.2;

opt.showLLH = 0;
opt.showIt = 0;
opt.numiter = 20;

szh(1) = uint16(size(rgbImg(1:sc:end,1:sc:end,:),1));
szh(2) = uint16(size(rgbImg(1:sc:end,1:sc:end,:),2));

% accumulate features
combFeat = getCombinedFeatures(allInfo, sc);

% Generate oversegmentation
display('Applying JCSA Clustering ...');
label = uint8(fusionBD_Color_3D_Axis(combFeat, kMax, [1 1 1], opt));
labelImg = reshape(label, szh(1), szh(2));

% post-processing
display('Applying Post-processing ...');
bdry = seg2bdry_2(labelImg);
regImg = bwlabeln(~bdry);
rgbImg = rgbImg(1:2:end, 1:2:end, :);

[edges, tmp, neighbors, tmpLabImg] = seg2fragments(double(regImg), rgbImg, 10);

display('Applying Region-Megring method ...');
if(strcmp(MethodType, 'rb_jcsa_rm'))
    % Compute parameters of the distributions
    parNormal = getKappaWMM(combFeat(:,7:9), tmpLabImg(:));
    parColor = getParamsDivergenceGMM(combFeat(:,1:3), tmpLabImg(:));
    
    % Set initial parameters for the neighborhood regions
    rgnb.nbsegs = getRegionNeighbors(tmpLabImg);
    rgnb.Div_N = ones(length(rgnb.nbsegs))*5000;
    rgnb.Div_C = ones(length(rgnb.nbsegs))*5000;
    rgnb.KappaMerged = ones(length(rgnb.nbsegs)) * -20;
    
    % Update parameters for the neighborhood regions
    for i=1:length(rgnb.nbsegs)
        adjs = rgnb.nbsegs{i};
        adjEdges = cell(length(adjs),1);
        
        for j = 1:length(adjs)
            commonEdge = intersect(neighbors.segment_fragments{i}, neighbors.segment_fragments{adjs(j)});
            commIndx = [];
            for k=1:length(commonEdge)
                commIndx = [commIndx; fliplr(floor(edges{commonEdge(k)}))];
            end
            rgnb.adjEdges{i,adjs(j)} = commIndx;
            
            % Compute kappa after it is merged with the neighbor cluster
            rgnb.KappaMerged(i,adjs(j)) = mergeKappaWMM(parNormal, [i,adjs(j)]);
            
            % Compute BD_normal among the neighbors
            rgnb.Div_N(i,adjs(j)) = parNormal.DLNF(i) - parNormal.DLNF(adjs(j)) - ( (parNormal.eta(i,:) - parNormal.eta(adjs(j),:)) * parNormal.theta_cl(adjs(j),:)');
            
            % Compute BD_color among the neighbors
            rgnb.Div_C(i,adjs(j)) = compute_Div_C(parColor, i, adjs(j));
        end
    end
    
    % Apply RM method
    img = getRAGMergeSegments_Edge_Plane_gui(combFeat, tmpLabImg, thOptions, parNormal, parColor, rgnb, allInfo, sc);
else
    img = tmpLabImg;
end

display('Displaying segmentation results ...');
segres = label2rgb(assignRandomLabel(img));
subplot(2, 3, [4:6]); imshow(segres); title('Segmented Image');