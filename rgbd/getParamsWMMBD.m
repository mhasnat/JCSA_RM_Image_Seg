function [params, X_ex] = getParamsWMMBD(X, clust, k)

% build extension matrices for X and mu : extension matrix
[n, d] = size(X);
a = 0.5;
c = d/2;

for j=1:k
     muT = mean(X(clust==j,:));
     mu(j,:) = muT ./ norm(muT);
end

[X_ex, ~] = getExtensionMatrices(X, mu, d);

alpha = zeros(1,k);

for j = 1:k
  indices = clust==j;  
  alpha(j) = sum(indices)/n;
  
  % parameters W (from mu) and D (from x)
  eta(j,:) = mean(X_ex(indices, :)); 
  normTheta(j) = getNormThetaNR(norm(eta(j,:)), d);
  
  % Compute R(normTheta)
  g_norm_theta(j) = (a/c) * (chgm(a+1, c+1, normTheta(j)) / chgm(a, c, normTheta(j)));
  R_norm_theta(j) = g_norm_theta(j) / normTheta(j);
  theta_cl(j, :) = eta(j, :) ./ R_norm_theta(j); % natural parameter
end

params.alpha = alpha;
params.eta = eta;
params.kappa = normTheta;
params.mu = mu;
params.theta_cl = theta_cl;
params.normTheta = normTheta;