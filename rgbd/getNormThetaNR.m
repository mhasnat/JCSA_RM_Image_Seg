function x = getNormThetaNR(v, d)
a = 0.5;
c = d/2;

kappa_upper_bound = 500;

x = 1; % initial value

diff = 100;
x_v(1) = 500000;
ii=2;

while(diff>0.001 && x<kappa_upper_bound)
    % calculate g(x)
    g_x = (a/c) * (chgm(a+1, c+1, x) / chgm(a, c, x));
    grad_g_x = ((1-(c/x))*g_x) + (a/x) - (g_x^2);
    
    x = x - ((g_x - v)/grad_g_x);
    
    x_v(ii) = x;
    
    diff = abs(x_v(ii-1) - x);
    
    ii = ii + 1;
end

end