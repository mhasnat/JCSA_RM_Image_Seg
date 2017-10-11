function DG_N = getBDVMFMM(Params, NormalTerm)

k = length(Params.alpha);
D = size(NormalTerm, 1);

LNF_vmf = zeros(1,k);
DLNF_vmf = zeros(1,k);
GDLNF_vmf = zeros(k, 3);
SST_MIN_EP_vmf = zeros(D, 3, k);
innerProdTerm = zeros(D, k);

DG_N = zeros(D, k);

% Compute Bragman Divergence for Normals
for j=1:k
    % Compute Bragman divergence
    LNF_vmf(j) = log((sinh(Params.normTheta(j))) / Params.normTheta(j));
    DLNF_vmf(j) = (Params.eta(j, :) * Params.theta_cl(j, :)') - LNF_vmf(j);
    GDLNF_vmf(j,:) = Params.theta_cl(j, :);
    
    SST_MIN_EP_vmf(:, :, j) = bsxfun(@minus, NormalTerm , Params.eta(j, :));
    innerProdTerm(:,j) = SST_MIN_EP_vmf(:, :, j) * GDLNF_vmf(j,:)';
    
    % Compute divergence
    DG_N(:,j) = -(DLNF_vmf(j) + innerProdTerm(:,j));
end