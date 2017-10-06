function Params = updateParametersBDWMM(AxisTerm, Params, probTerm)

d = size(Params.mu,2);
a = 0.5;
c = d/2;

k=length(Params.alpha);

sumProbTerm = sum(probTerm);

R_norm_theta = zeros(1,k);
g_norm_theta = zeros(1,k);

for j=1:k
    Params.eta(j, :) = sum(bsxfun(@times, probTerm(:,j) , AxisTerm)) ./ sumProbTerm(j);
    
    % Convert to source parameters (mu, kappa) from expectation parameter (eta)
    Params.normTheta(j) = getNormThetaNR(norm(Params.eta(j,:)), d);
    
    % Compute R(normTheta)
    g_norm_theta(j) = (a/c) * (chgm(a+1, c+1, Params.normTheta(j)) / chgm(a, c, Params.normTheta(j)));
    R_norm_theta(j) = g_norm_theta(j) / Params.normTheta(j);
    Params.theta_cl(j, :) = Params.eta(j, :) ./ R_norm_theta(j); % natural parameter
    
    % update sorce parameters
    Params.kappa(j) = Params.normTheta(j);
    W = Params.theta_cl(j, :) ./ Params.normTheta(j);
    Params.mu(j, :) = getMufromW(W, size(Params.mu,2));
end