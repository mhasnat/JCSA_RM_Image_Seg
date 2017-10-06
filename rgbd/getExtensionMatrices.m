function [X_ex M_ex] = getExtensionMatrices(X, mu, d)
% build extension matrices for X and mu : extension matrix
B = sortrows(combnk(1:d,2));

X_ex = X.^2;
M_ex = mu.^2;

for ii=1:size(B,1)
    X_ex(:,d+ii) = sqrt(2) * (X(:,B(ii,1)) .* X(:,B(ii,2)));
    M_ex(:,d+ii) = sqrt(2) * (mu(:,B(ii,1)) .* mu(:,B(ii,2)));
end