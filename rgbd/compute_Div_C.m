function Div_C = compute_Div_C(params, i, j)

diff_ep_vec = params.theta1(:,i) - params.theta1(:,j);
diff_ep_Mat = params.theta2(:,:,i) - params.theta2(:,:,j);
innerProdTerm = trace(diff_ep_Mat*params.GLNF{j,2}') + diff_ep_vec'*params.GLNF{j,1};

Div_C = params.LNF(i) - params.LNF(j) - innerProdTerm;