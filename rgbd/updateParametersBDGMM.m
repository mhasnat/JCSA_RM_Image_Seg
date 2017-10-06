function Params = updateParametersBDGMM(Params, probTerm)

k=length(Params.alpha);
D = size(probTerm,1);
sumProbTerm = sum(probTerm);
dim = size(Params.eta,2);

for j=1:k
    % Expectation parameters
    Params.eta(j,:) = sum(bsxfun(@times, probTerm(:,j) , Params.sst1)) ./ sumProbTerm(j);
    
    tmp1 = reshape(Params.sst2, dim^2, D);
    tmp2 = sum(bsxfun(@times, tmp1', probTerm(:,j)))./ sumProbTerm(j);
    Params.H(:,:,j) = reshape(tmp2, dim, dim);
    clear tmp1 tmp3
    
    % Natural parameters
    tmp3 = inv(Params.H(:,:,j) + Params.eta(j,:)'*Params.eta(j,:));
    Params.theta1(:,j) = - tmp3 * Params.eta(j,:)';
    Params.theta2(:,:,j) = - 0.5 * tmp3;
    
    % Source parameters
    Params.mu(j,:) = Params.eta(j,:);
    Params.sigma(:,:,j) = -(Params.H(:,:,j) + Params.eta(j,:)'*Params.eta(j,:));
end