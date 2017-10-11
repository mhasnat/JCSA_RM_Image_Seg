function params = getParamsVMFBD(vectors, clust, k)

if(nargin<3)
    k = length(unique(clust));
end

numOfDataSample = size(vectors,1);

% 
for i=1:k
    numOfDataPointToThisCluster = length(find(clust==i));
    
    alpha(i) = numOfDataPointToThisCluster / numOfDataSample; % weight
    
    % Expectation parameter computation
    indx = find(clust==i);
    dataClust = vectors(indx, :);
    
    sufStat(i, :) = sum(dataClust);
    
    eta(i, :) = sufStat(i, :) / length(dataClust);
    normEta(i) = sqrt(eta(i, :) * eta(i, :)');
    normTheta(i) = getThetaFromEta(normEta(i));
    
    % Compute R(normTheta)
    R_norm_theta(i) = ((1/tanh(normTheta(i))) - (1/normTheta(i))) / normTheta(i);
    theta_cl(i, :) = eta(i, :) ./ R_norm_theta(i); % natural parameter
    
    % Source parameters
    kappa(i) = normTheta(i);
    mu(i, :) = theta_cl(i, :) ./ normTheta(i);
    
    %For Computing Bragman divergence
    LNF(i) = log((sinh(normTheta(i))) / normTheta(i)); % The log normalizing function   
    DLNF(i) = (eta(i, :) * theta_cl(i, :)') - LNF(i); % Dual (expectation parameter) of the log normalizing function
end

params.alpha = alpha;
params.eta = eta;
params.theta_cl = theta_cl;
params.normTheta = normTheta;
params.kappa = kappa;
params.mu = mu;

params.LNF = LNF;
params.DLNF = DLNF;