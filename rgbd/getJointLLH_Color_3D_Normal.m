function LLH = getJointLLH_Color_3D_Normal(alpha, ColorInfo, DepthInfo, NormalInfo)

ColorTerm = ColorInfo.Term;
col_mu_gmm = ColorInfo.mu;
col_sigma = ColorInfo.sigma;

DepthTerm = DepthInfo.Term;
dep_mu_gmm = DepthInfo.mu;
dep_sigma = DepthInfo.sigma;

NormalTerm = NormalInfo.Term;
mu_vmf = NormalInfo.mu;
kappa = NormalInfo.kappa;


D= size(ColorTerm,1);
k = length(alpha);

% Compute new log likelihood --> GMM
Color_logGaussian = zeros(D,k);
for j = 1:k
    Color_logGaussian(:,j) = loggausspdf(ColorTerm',col_mu_gmm(j,:)',col_sigma(:,:,j));
end

Depth_logGaussian = zeros(D,k);
for j = 1:k
    Depth_logGaussian(:,j) = loggausspdf(DepthTerm',dep_mu_gmm(j,:)',dep_sigma(:,:,j));
end

% Compute new log likelihood --> vMFMM
logNormTerm = log(kappa) - log(sinh(kappa));
logExpTerm  = bsxfun(@times, kappa,  (mu_vmf * NormalTerm')');
logVMF = bsxfun(@plus, logNormTerm , logExpTerm);

% add all
logGaussVMF = Color_logGaussian + Depth_logGaussian + logVMF;

logPDF = bsxfun(@plus,logGaussVMF,log(alpha));
logMixtureLikelihoodPerComponent = logsumexp(logPDF,2); % log(sum(exp(logPDF),2));

LLH = sum(logMixtureLikelihoodPerComponent);

end