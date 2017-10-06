function params = getParamsDivergenceGMM(vectors, clust, k)

if(nargin<3)
    k = length(unique(clust));
end

dim = size(vectors, 2);
numOfDataSample = size(vectors,1);

% Sufficient statistics for every sample
sst1 = vectors;
% if(size(vectors,2)<=3)
%     sst2 = computeSST(vectors); % mex function
% else
    for i=1:size(vectors,1)
        sst2(:,:,i) = -vectors(i,:)' * vectors(i,:);
    end
% end

% 
for j=1:k
    indx = find(clust==j);
    alpha(j) = length(indx) / numOfDataSample; % weight
    
    % data
    dataClust = vectors(indx, :);
    
    % source parameters
    mu(j,:) = mean(dataClust);
    sigma(:,:,j) = cov(dataClust);
    invSig(:,:,j) = inv(sigma(:,:,j));
    
    % Natural parameters -> Gradiant of dual log normalizing function
    theta1(:,j) = invSig(:,:,j)*mu(j,:)';
    theta2(:,:,j) = 0.5 * invSig(:,:,j);
    
    % Expectation parameters
    eta(j,:) = mu(j,:);
    H(:,:,j) = -(sigma(:,:,j) + (mu(j,:)' * mu(j,:)));
    
    % For Computing Bragman divergence
    tmp1 = inv(theta2(:,:,j));
    tmp2 = tmp1 * theta1(:,j);
    
    % Log normalizing function
    LNF(j) = (0.25 * trace(tmp2 * theta1(:,j)')) - (0.5*log(det(theta2(:,:,j)))) + (dim*0.5*log(pi));
    
    % Gradient of Dual log normalizing function
    GLNF{j,1} = 0.5 * tmp2;
    GLNF{j,2} = -(0.5*tmp1) - (0.25*tmp2*tmp2');
end

params.alpha = alpha;
params.mu = mu;
params.sigma = sigma;
params.theta1 = theta1;
params.theta2 = theta2;
params.eta = eta;
params.H = H;
params.sst1 = sst1;
params.sst2 = sst2;

params.LNF = LNF;
params.GLNF = GLNF;