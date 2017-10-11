function Up_kappa = mergeKappaVMFMM(Params, subset_cl)

% Right sided Bregman centroid computation
Up_Params.alpha = sum(Params.alpha(subset_cl));

% Expectation (from natural) update
Up_Params.eta = sum(bsxfun(@times, Params.alpha(subset_cl)', Params.eta(subset_cl,:))) ./ Up_Params.alpha;

% Convert to source parameters (mu, kappa) from expectation parameter (eta)
Up_kappa = getThetaFromEta(norm(Up_Params.eta));