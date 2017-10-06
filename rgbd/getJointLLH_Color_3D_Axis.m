function LLH = getJointLLH_Color_3D_Axis(alpha, ColorInfo, DepthInfo, AxisInfo)

ColorTerm = ColorInfo.Term;
col_mu_gmm = ColorInfo.mu;
col_sigma = ColorInfo.sigma;

DepthTerm = DepthInfo.Term;
dep_mu_gmm = DepthInfo.mu;
dep_sigma = DepthInfo.sigma;

AxisTerm = AxisInfo.Term;
mu_wmm = AxisInfo.mu;
kappa = AxisInfo.kappa;


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

% Compute new log likelihood --> WMM
logNormTerm = zeros(1,k);
d = size(AxisInfo.mu,2);

for j=1:k
    logNormTerm(j)=log(1/chgm(0.5,d/2,kappa(j)));
end
logExpTerm  = bsxfun(@times, kappa,  (mu_wmm * AxisTerm')'.^2);
logWMM = bsxfun(@plus, logNormTerm , logExpTerm);
ind = logWMM>700;
logWMM(ind) = logWMM(ind) - 1000;

% add all
logGaussWMM = Color_logGaussian + Depth_logGaussian + logWMM;

logPDF = bsxfun(@plus,logGaussWMM,log(alpha));
logMixtureLikelihoodPerComponent = logsumexp(logPDF,2); % log(sum(exp(logPDF),2));

LLH = sum(logMixtureLikelihoodPerComponent);

end