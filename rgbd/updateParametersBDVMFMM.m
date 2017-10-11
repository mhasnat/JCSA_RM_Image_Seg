function Params = updateParametersBDVMFMM(NormalTerm, Params, probTerm)

k=length(Params.alpha);
sumProbTerm = sum(probTerm);
% D = size(probTerm,1);
% dim = size(Params.eta,2);

normEta = zeros(1,k);
R_norm_theta = zeros(1,k);

for j=1:k
    Params.eta(j, :) = sum(bsxfun(@times, probTerm(:,j) , NormalTerm)) ./ sumProbTerm(j);
    
    % Convert to source parameters (mu, kappa) from expectation parameter
    % (eta_vmf)
    normEta(j) = sqrt(Params.eta(j, :) * Params.eta(j, :)');
    Params.normTheta(j) = getThetaFromEta(normEta(j));
    
    % Compute R(normTheta)
    R_norm_theta(j) = ((1/tanh(Params.normTheta(j))) - (1/Params.normTheta(j))) / Params.normTheta(j);
    Params.theta_cl(j, :) = Params.eta(j, :) ./ R_norm_theta(j);
    
    Params.kappa(j) = Params.normTheta(j);
    Params.mu(j, :) = Params.theta_cl(j, :) ./ Params.normTheta(j);
end