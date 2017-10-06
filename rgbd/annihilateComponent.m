function Params = annihilateComponent(Params, indxAnn, type)

if(strcmp(type, 'gmm'))
    % Gaussian components
    Params.eta(indxAnn, :) = [];
    Params.theta1(:, indxAnn) = [];
    Params.theta2(:, :, indxAnn) = [];
    Params.H(:, :, indxAnn) = [];
    
    Params.mu(indxAnn,:) = [];
    Params.sigma(:,:,indxAnn) = [];
elseif(strcmp(type, 'vmfmm'))
    % vMF components
    Params.eta(indxAnn, :) = [];
    Params.theta_cl(indxAnn, :) = [];
    Params.mu(indxAnn, :) = [];
    Params.kappa(indxAnn) = [];
    Params.normTheta(indxAnn) = [];
end

Params.alpha(indxAnn) = [];