function [ParamsColor ParamsDepth ParamsNormal] = getAnnihilatedParams_JBD(ParamsColor, ParamsDepth, ParamsNormal, invalidIndx)

ParamsColor.alpha(invalidIndx) = [];
ParamsColor.alpha = ParamsColor.alpha / sum(ParamsColor.alpha);

ParamsColor.mu(invalidIndx,: ) = [];
ParamsColor.sigma(:,:, invalidIndx) = [];
ParamsColor.theta1(:,invalidIndx) = [];
ParamsColor.theta2(:,:, invalidIndx) = [];
ParamsColor.eta(invalidIndx,:) = [];
ParamsColor.H(:,:, invalidIndx) = [];

ParamsDepth.alpha(invalidIndx) = [];
ParamsDepth.alpha = ParamsDepth.alpha / sum(ParamsDepth.alpha);

ParamsDepth.mu(invalidIndx,: ) = [];
ParamsDepth.sigma(:,:, invalidIndx) = [];
ParamsDepth.theta1(:,invalidIndx) = [];
ParamsDepth.theta2(:,:, invalidIndx) = [];
ParamsDepth.eta(invalidIndx,:) = [];
ParamsDepth.H(:,:, invalidIndx) = [];

ParamsNormal.alpha(invalidIndx) = [];
ParamsNormal.alpha = ParamsNormal.alpha / sum(ParamsNormal.alpha);

ParamsNormal.mu(invalidIndx,: ) = [];
ParamsNormal.kappa(invalidIndx) = [];
ParamsNormal.eta(invalidIndx,:) = [];
ParamsNormal.theta_cl(invalidIndx,:) = [];
ParamsNormal.normTheta(invalidIndx) = [];

end