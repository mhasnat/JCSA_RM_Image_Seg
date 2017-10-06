function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);

[U,p]= chol(Sigma);

if p ~= 0
    %     error('ERROR: Sigma is not PD.');
    [V,D] = eig(Sigma);
    Sigma = V * diag(max(diag(D),eps)) / V;
    [U,p]= chol(Sigma);
end

Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;