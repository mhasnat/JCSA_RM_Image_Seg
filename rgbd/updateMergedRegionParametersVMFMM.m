function Up_Params = updateMergedRegionParametersVMFMM(Up_Params, Params, subset_cl, indx)

% Right sided Bregman centroid computation
Up_Params.alpha(indx) = sum(Params.alpha(subset_cl));

% Expectation (from natural) update
Up_Params.eta(indx,:) = sum(bsxfun(@times, Params.alpha(subset_cl)', Params.eta(subset_cl,:))) ./ Up_Params.alpha(indx);

% Convert to source parameters (mu, kappa) from expectation parameter
% (eta)
Up_Params.normTheta(indx) = getThetaFromEta(norm(Up_Params.eta(indx, :)));

% Compute R(normTheta)
R_norm_Up_theta = ((1/tanh(Up_Params.normTheta(indx))) - (1/Up_Params.normTheta(indx))) / Up_Params.normTheta(indx);
Up_Params.theta_cl(indx, :) = Up_Params.eta(indx, :) ./ R_norm_Up_theta; % natural parameter


% Update Source Parameters
Up_Params.kappa(indx) = Up_Params.normTheta(indx);
Up_Params.mu(indx, :) = Up_Params.theta_cl(indx, :) ./ Up_Params.normTheta(indx);
  
%For Computing Bragman divergence
Up_Params.LNF(indx) = log((sinh(Up_Params.normTheta(indx))) / Up_Params.normTheta(indx)); % The log normalizing function   
Up_Params.DLNF(indx) = (Up_Params.theta_cl(indx, :) * Up_Params.eta(indx,:)') - Up_Params.LNF(indx); % Dual (expectation parameter) of the log normalizing function