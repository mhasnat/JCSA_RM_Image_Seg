function [mu, clust] = kmeans_Color_3D_Normal_plus_plus(vectors, k, wt, SP)

showRes=0;
if(nargin>3)
    if(~isempty(SP))
        showRes=1;
    end
end

[D,V] = size(vectors);

% Set the indices
colIndx = 1:3;
spatialIndx = 4:6;
normalIndx = 7:9;

% ----------------------------------------------------------------
% Getting the initial means using different features and different
% clustering algorithms
% ----------------------------------------------------------------

% For color features
% K-means++
[label,~] = kmeanspp(vectors(:, colIndx)', k);
% [label,~] = getPCAbasedCenters(vectors(:, colIndx), k);
for i=1:k
    muT = mean(vectors(find(label==i),colIndx));
    mu(i,colIndx) = muT;
end

% For spatial features (3D points)
% K-means++
[label,~] = kmeanspp(vectors(:, spatialIndx)', k);
% [label,~] = getPCAbasedCenters(vectors(:, spatialIndx), k);
for i=1:k
    muT = mean(vectors(find(label==i),spatialIndx));
    mu(i,spatialIndx) = muT;
end

% For structural features (Image Normal)
% K-means++
[label,~] = kmeansPP_CD(vectors(:, normalIndx)', k);
% [label,~] = getPCAbasedCenters(vectors(:, normalIndx), k);
for i=1:k
    muT = mean(vectors(find(label==i), normalIndx));
    mu(i,normalIndx) = muT ./ norm(muT);
end

% mu

% % --------------------------------------
% % Getting the means from combined kmeans
% % --------------------------------------
diff    = 10000;% 1;
epsilon = 0.1;
value   = 100;
iteration = 1;
maxnumIter = 50;
while (diff > epsilon & maxnumIter>iteration)
    iteration = iteration+1;
    oldvalue      = value;
    
    % Color features distances
    if(wt(1)~=0)
        temp = zeros(D,k);
        for ti = 1:k
            temp(:,ti) = sqrt(sum(bsxfun(@minus, vectors(:, colIndx), mu(ti,colIndx)).^2, 2));
        end
        
        distMat_color =  temp;
    else
        distMat_color = 0;
    end
    
    % Spatial (3D) features distances
    if(wt(2)~=0)
        temp = zeros(D,k);
        for ti = 1:k
            temp(:,ti) = sqrt(sum(bsxfun(@minus, vectors(:, spatialIndx), mu(ti,spatialIndx)).^2, 2));
        end
        
        distMat_spatial =  temp;
    else
        distMat_spatial = 0;
    end
    
    % Structural features distances
    if(wt(3)~=0)
        distMat_normal        =  1 - (vectors(:, normalIndx) * mu(:,normalIndx)').^2; % [0 2] -> range of distance for normal
    else
        distMat_normal = 0;
    end
    
    % assign points to nearest cluster (based on minimum distance)
    distMat = distMat_color.*wt(1) + distMat_spatial.*wt(2) + distMat_normal.*wt(3); % the distance function
    [distMin,clust] =  min(distMat,[],2); % the minimum as assignment criterion
    
    % Display results
    if(showRes==1)
        clstImg = assignClusterLabels2Superpixels(SP, clust);
        imshow(label2rgb(clstImg));
    end
    
    % compute objective function value
    value         = sum(distMin);
    
    % compute new cluster centroids
    for h=1:k
        try
            % The mean color
            mu(h,colIndx) = mean(vectors(find(clust==h),colIndx));
            
            % The mean 3D space
            mu(h,spatialIndx) = mean(vectors(find(clust==h),spatialIndx));
            
            % The mean directions
            sumVec(h,:)  = sum(vectors(find(clust==h),normalIndx),1);
            mu(h,normalIndx)      = sumVec(h,:)./ norm(sumVec(h,:));
        catch
            % display('I am trying to catch the error ...');
            succ = 0;
            return;
        end
    end
    
    % Remove empty clusters if any
    if(~isempty(find(isnan(sum(mu, 2)))))
        tindices = find(isnan(sum(mu, 2)));
        mu(tindices, :) = [];
        k = k-length(tindices);
    end
    
    diff = abs(value - oldvalue);
end