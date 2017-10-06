function Up_kappa = mergeKappaWMM(Params, subset_cl)

d = size(Params.mu,2);
a = 0.5;
c = d/2;

% Right sided Bregman centroid computation
Up_Params.alpha = sum(Params.alpha(subset_cl));

% Expectation (from natural) update
Up_Params.eta = sum(bsxfun(@times, Params.alpha(subset_cl)', Params.eta(subset_cl,:))) ./ Up_Params.alpha;

% Convert to source parameters (mu, kappa) from expectation parameter (eta)
Up_kappa = getNormThetaNR(norm(Up_Params.eta), d);