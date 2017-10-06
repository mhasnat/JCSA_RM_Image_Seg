function Div_N = compute_Div_N(parNormal, i, j)

Div_N = parNormal.DLNF(i) - parNormal.DLNF(j) - ( (parNormal.eta(i,:) - parNormal.eta(j,:)) * parNormal.theta_cl(j,:)');