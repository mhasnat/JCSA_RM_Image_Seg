function DG_A = getBDWMM(Params, AxisTerm)

k = length(Params.alpha);
[N, p] = size(AxisTerm);
d = size(Params.mu,2);

LNF = zeros(1,k);
DLNF = zeros(1,k);
GDLNF = zeros(k, p);
sufStat_minus_expectParam = zeros(N, p, k);
innerProdTerm = zeros(N, k);

DG_A = zeros(N, k);

% Compute Bragman Divergence for Axes
for j=1:k
    % Compute Bragman divergence
    LNF(j) = log(chgm(0.5, d/2, Params.normTheta(j))); % The log normalizing function
    DLNF(j) = (Params.theta_cl(j,:) * Params.eta(j,:)') - LNF(j); % Dual (expectation parameter) of the log normalizing function
    
    GDLNF(j,:) = Params.theta_cl(j, :);
    
    sufStat_minus_expectParam(:, :, j) = bsxfun(@minus, AxisTerm , Params.eta(j, :));
    innerProdTerm(:,j) = sufStat_minus_expectParam(:, :, j) * GDLNF(j,:)';
    
    % Compute divergence
    DG_A(:,j) = -(DLNF(j) + innerProdTerm(:,j));
end