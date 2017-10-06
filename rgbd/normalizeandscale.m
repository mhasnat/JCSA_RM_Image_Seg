function opvec = normalizeandscale(vector)

numCol = size(vector, 2);

% Normalize
sigma_c = numCol/ sum(std(vector));
vector = sigma_c .* vector;

% Scale
vector = bsxfun(@minus, vector, min(vector));
opvec = bsxfun(@rdivide, vector, max(vector));