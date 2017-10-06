function mu = getMufromW(W, d)

% get mu initial
mu_init = sqrt(W(1:d));
W = W>0;
% set the sign of mu
mu = mu_init;
if(W(4)==1 & W(5)==1 & W(6)==1)
elseif(W(4)==1 & W(5)==0 & W(6)==0)
    mu(3) = -mu_init(3);
elseif(W(4)==0 & W(5)==0 & W(6)==1)
    mu(1) = -mu_init(1);
elseif(W(4)==0 & W(5)==1 & W(6)==0)
    mu(2) = -mu_init(2);
end

mu = mu / norm(mu);

end