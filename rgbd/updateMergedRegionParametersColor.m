function Up_Params = updateMergedRegionParametersColor(Up_Params, Params, subset_cl, indx)
dim = size(Params.mu,2);

% Right sided Bregman centroid computation
Up_Params.alpha(indx) = sum(Params.alpha(subset_cl));

% Natural params centroid merge
tmp1 = zeros(size(Params.theta1,1),1);
tmp2 = zeros(size(Params.theta2,1),size(Params.theta2,2));
for cnt=1:length(subset_cl)
    tmp1 = tmp1 + (Params.alpha(subset_cl(cnt)) .* Params.theta1(:, subset_cl(cnt)));
    tmp2 = tmp2 + (Params.alpha(subset_cl(cnt)) .* Params.theta2(:,:, subset_cl(cnt)));
end
Up_Params.theta1(:,indx) = tmp1 ./ Up_Params.alpha(indx);
Up_Params.theta2(:,:,indx) = tmp2 ./ Up_Params.alpha(indx);

% Expectation (from natural) update
th2Invth1 = inv(Up_Params.theta2(:,:,indx)) * Up_Params.theta1(:,indx);
Up_Params.eta(indx,:) = 0.5 * th2Invth1;
Up_Params.H(:,:,indx) = -0.5*inv(Up_Params.theta2(:,:,indx)) - 0.25*th2Invth1*th2Invth1';

% Source (from expectation) update
Up_Params.mu(indx,:) = Up_Params.eta(indx,:);
Up_Params.sigma(:,:,indx) = -(Up_Params.H(:,:,indx) + Up_Params.eta(indx,:)'*Up_Params.eta(indx,:));

% For Computing Bragman divergence
tmp1 = inv(Up_Params.theta2(:,:,indx));
tmp2 = tmp1 * Up_Params.theta1(:,indx);

% Log normalizing function
Up_Params.LNF(indx) = (0.25 * trace(tmp2 * Up_Params.theta1(:,indx)')) - (0.5*log(det(Up_Params.theta2(:,:,indx)))) + (dim*0.5*log(pi));

% Gradient of Dual log normalizing function
Up_Params.GLNF{indx,1} = 0.5 * tmp2;
Up_Params.GLNF{indx,2} = -(0.5*tmp1) - (0.25*tmp2*tmp2');

% % % % Right sided Bregman centroid computation
% % % Up_Params.alpha(indx) = sum(Params.alpha(subset_cl));
% % %
% % % % Expectation (from natural) update
% % % Up_Params.eta(indx,:) = sum(bsxfun(@times, Params.alpha(subset_cl)', Params.eta(subset_cl,:))) ./ Up_Params.alpha(indx);
% % %
% % % % Convert to source parameters (mu, kappa) from expectation parameter
% % % % (eta)
% % % Up_Params.normTheta(indx) = getNormThetaNR(norm(Up_Params.eta(indx, :)), d);
% % %
% % % % Compute R(normTheta)
% % % g_norm_theta_Up = (a/c) * (chgm(a+1, c+1, Up_Params.normTheta(indx)) / chgm(a, c, Up_Params.normTheta(indx)));
% % % R_norm_Up_theta = g_norm_theta_Up / Up_Params.normTheta(indx);
% % % Up_Params.theta_cl(indx, :) = Up_Params.eta(indx, :) ./ R_norm_Up_theta; % natural parameter
% % %
% % % % update source parameters
% % % Up_Params.kappa(indx) = Up_Params.normTheta(indx);
% % % W = Up_Params.theta_cl(indx, :) ./ Up_Params.normTheta(indx);
% % % Up_Params.mu(indx, :) = getMufromW(W, d);