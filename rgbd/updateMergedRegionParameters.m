function Up_Params = updateMergedRegionParameters(Up_Params, Params, subset_cl, indx)


d = size(Params.mu,2);
a = 0.5;
c = d/2;

% Right sided Bregman centroid computation
Up_Params.alpha(indx) = sum(Params.alpha(subset_cl));

% Expectation (from natural) update
Up_Params.eta(indx,:) = sum(bsxfun(@times, Params.alpha(subset_cl)', Params.eta(subset_cl,:))) ./ Up_Params.alpha(indx);

% Convert to source parameters (mu, kappa) from expectation parameter
% (eta)
Up_Params.normTheta(indx) = getNormThetaNR(norm(Up_Params.eta(indx, :)), d);

% Compute R(normTheta)
g_norm_theta_Up = (a/c) * (chgm(a+1, c+1, Up_Params.normTheta(indx)) / chgm(a, c, Up_Params.normTheta(indx)));
R_norm_Up_theta = g_norm_theta_Up / Up_Params.normTheta(indx);
Up_Params.theta_cl(indx, :) = Up_Params.eta(indx, :) ./ R_norm_Up_theta; % natural parameter

% update source parameters
Up_Params.kappa(indx) = Up_Params.normTheta(indx);
W = Up_Params.theta_cl(indx, :) ./ Up_Params.normTheta(indx);
Up_Params.mu(indx, :) = getMufromW(W, d);

%For Computing Bragman divergence
Up_Params.LNF(indx) = log(chgm(0.5, d/2, Up_Params.normTheta(indx))); % The log normalizing function
Up_Params.DLNF(indx) = (Up_Params.theta_cl(indx, :) * Up_Params.eta(indx,:)') - Up_Params.LNF(indx); % Dual (expectation parameter) of the log normalizing function