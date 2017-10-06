function DG_C = getBDUGMM(ParamsColor)

dim_gmm = size(ParamsColor.eta, 2);
k = length(ParamsColor.alpha);
D = size(ParamsColor.sst1,1);

DLNF_gmm = zeros(1,k);
GDLNF_gmm = cell(k,2);
SST_MIN_EP_gmm = cell(k,2);
traceTerm = zeros(D, k);
innerProdTerm = zeros(D, k);
DG_C = zeros(D, k);


% Compute Bragman Divergence for Colors
for j=1:k
    % Dual log normalizing function
    DLNF_gmm(j) = -0.5 * (log(1.00001 +(ParamsColor.eta(j,:)*inv(ParamsColor.H(:,:,j))*ParamsColor.eta(j,:)')) + log(det(-ParamsColor.H(:,:,j))) + (dim_gmm*log(2*pi*2.71)));
    
    % Gradient of Dual log normalizing function
    GDLNF_gmm{j,1} = ParamsColor.theta1(:,j);
    GDLNF_gmm{j,2} = ParamsColor.theta2(:,:,j);
    
    % Sufficien Statistics minus Expectation Parameters
    SST_MIN_EP_gmm{j,1} = bsxfun(@minus, ParamsColor.sst1, ParamsColor.eta(j, :));
    SST_MIN_EP_gmm{j,2} = bsxfun(@minus, ParamsColor.sst2, ParamsColor.H(:,:,j));
    
    % need to efficiently handle
    %     for i=1:D
    %         traceTerm(i,j) = trace(SST_MIN_EP_gmm{j,2}(:,:,i)*GDLNF_gmm{j,2}');
    %     end
    traceTerm(:,j) = computeTraceTerm(SST_MIN_EP_gmm{j,2}, GDLNF_gmm{j,2}');
    innerProdTerm(:,j) = traceTerm(:,j) + SST_MIN_EP_gmm{j,1}*GDLNF_gmm{j,1};
    
    % Compute divergence
    DG_C(:,j) = -(DLNF_gmm(j) + innerProdTerm(:,j));
end