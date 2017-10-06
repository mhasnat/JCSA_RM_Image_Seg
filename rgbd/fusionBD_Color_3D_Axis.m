function clust = fusionBD_Color_3D_Axis(vectors, k, wt, opt)

[D,~] = size(vectors);

% Basic information
ColorTerm = vectors(:,1:3); % Color data
DepthTerm = vectors(:,4:6); % 3D data
NormalTerm = vectors(:,7:9); % Normal data

%% Step 1: Initialization
% Getting the initial means using k-means
[~, clust] = kmeans_Color_3D_Normal_plus_plus(vectors, k, wt);
% [~, clust] = kmeans_Color_3D_Normal(vectors, k, wt);
% clust = kmeans(vectors, k);

% % % % For color features - convert the parameters in natural parameter space
ParamsColor = getParamsGMMBD(ColorTerm, clust, k);
invalidIndx1 = getInvalidIndices(ParamsColor, 'gmm');

% % % % For Depth (3D spatial) features - convert the parameters in natural parameter space
ParamsDepth = getParamsGMMBD(DepthTerm, clust, k);
invalidIndx2 = getInvalidIndices(ParamsDepth, 'gmm');

% % % % For structural features - convert the parameters in expectation parameter space
[ParamsAxis, AxisTerm] = getParamsWMMBD(NormalTerm, clust, k);
invalidIndx3 = getInvalidIndices(ParamsAxis, 'wmm');

invalidIndx = unique([invalidIndx1 invalidIndx2 invalidIndx3]);

if(~isempty(invalidIndx))
    [ParamsColor, ParamsDepth, ParamsAxis] = getAnnihilatedParams_JBD(ParamsColor, ParamsDepth, ParamsAxis, invalidIndx);
end

k = k - length(invalidIndx);

% Prior class probability
alpha = ParamsAxis.alpha;

% % --------------------------------------
% % Start EM alorithm
% % --------------------------------------
MAX_ITERATIONS = opt.numiter;
logLikelihoodNew = 100;
logLikelihoodOld = 0;
logLikelihoodThreshold = 0.001;

iterations = 1;
diffLLH = abs(logLikelihoodOld - logLikelihoodNew);
while(iterations<MAX_ITERATIONS && diffLLH>logLikelihoodThreshold)
    if(opt.showIt)
        display(num2str(iterations));
    end
    
    logLikelihoodOld = logLikelihoodNew;
    
    %% Expectation Step
    DG_Gl = zeros(D, k);
    
    expTerm = zeros(D, k);
    probTerm = zeros(D, k);
        
    % Compute Bragman Divergence for Colors
    DG_C = getBDGMM(ParamsColor);
    
    % Compute Bragman Divergence for Depth
    DG_D = getBDGMM(ParamsDepth);
        
    % Compute Bragman Divergence for Normals
    DG_A = getBDWMM(ParamsAxis, AxisTerm);
    
    % Compute Combined Bregman Divergence and the Posterior Probability  
    for j=1:k
        % Compute Bragman divergence
        DG_Gl(:,j) = DG_C(:,j) + DG_D(:,j) + DG_A(:,j);
        
        % Compute Posterior Probability
        expTerm(:,j) = exp(-DG_Gl(:,j));
        probTerm(:, j) = alpha(j) * expTerm(:,j);
    end
    
    probTerm(isinf(probTerm)) = 1e7; % Just a fix
    
    probTerm = bsxfun(@rdivide, probTerm, sum(probTerm, 2));
    
    % Another fix
    nanIndices = isnan(probTerm);
    probTerm(nanIndices) = 1/k;
    
    [~, clust] = max(probTerm,[], 2);
    
    %% Maximization Step
    
    sumProbTerm = sum(probTerm);
    
    % Update weight
    alpha = sumProbTerm ./ size(probTerm,1);
    ParamsColor.alpha = alpha;
    ParamsDepth.alpha = alpha;
    ParamsAxis.alpha = alpha;
    
    % Annihiliate Components
    indxAnn = find(alpha < 0.001);
    
    if(~isempty(indxAnn))
        % GMM components
        ParamsColor = annihilateComponent(ParamsColor, indxAnn, 'gmm');
        DG_C(:,indxAnn) = [];
        
        ParamsDepth = annihilateComponent(ParamsDepth, indxAnn, 'gmm');
        DG_D(:,indxAnn) = [];
        
        % vMF components
        ParamsAxis = annihilateComponent(ParamsAxis, indxAnn, 'vmfmm');
        DG_A(:,indxAnn) = [];
              
        % common
        probTerm(:,indxAnn) = [];
        sumProbTerm(indxAnn) = [];
        k=k-length(indxAnn);
        alpha(indxAnn) = [];
        alpha = alpha ./ sum(alpha);
        ParamsColor.alpha = alpha;
        ParamsDepth.alpha = alpha;
        ParamsAxis.alpha = alpha;
    end
    
    % Update parameter --> 
    % Color with GMM
    ParamsColor = updateParametersBDGMM(ParamsColor, probTerm);
    invalidIndx1 = getInvalidIndices(ParamsColor, 'gmm');
    
    % Depth with GMM
    ParamsDepth = updateParametersBDGMM(ParamsDepth, probTerm);
    invalidIndx2 = getInvalidIndices(ParamsDepth, 'gmm');
    
    % Update parameter --> for Normal with WMM
    ParamsAxis = updateParametersBDWMM(AxisTerm, ParamsAxis, probTerm);
    invalidIndx3 = getInvalidIndices(ParamsAxis, 'vmfmm');
    
    invalidIndx = unique([invalidIndx1 invalidIndx2 invalidIndx3]);
    
    if(~isempty(invalidIndx))
        [ParamsColor, ParamsDepth, ParamsAxis] = getAnnihilatedParams_JBD(ParamsColor, ParamsDepth, ParamsAxis, invalidIndx);
        alpha = ParamsColor.alpha;
    end
    
    k = k - length(invalidIndx);

    % compute log likelihood
    ColorInfo.Term = ColorTerm;
    ColorInfo.mu = ParamsColor.mu;
    ColorInfo.sigma = ParamsColor.sigma;
    
    DepthInfo.Term = DepthTerm;
    DepthInfo.mu = ParamsDepth.mu;
    DepthInfo.sigma = ParamsDepth.sigma;
        
    AxisInfo.Term = NormalTerm;
    AxisInfo.mu = ParamsAxis.mu;
    AxisInfo.kappa = ParamsAxis.kappa;

    if(iterations==10)
        pp = 0;
    end
    %     iterations
    logLikelihoodNew = getJointLLH_Color_3D_Axis(alpha, ColorInfo, DepthInfo, AxisInfo);
    LLH(iterations) = -logLikelihoodNew; % negative log likelihood
    
    diffLLH = abs(logLikelihoodOld - logLikelihoodNew);
    iterations = iterations+1;
end

if(opt.showLLH)
    figure(1); plot(LLH, '--rs','LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',5);
end