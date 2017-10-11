function x = getThetaFromEta(D)
x = 1; % initial value

diff = 100;
x_v(1) = 500000;
ii=2;

while(diff>0.001)
    a = tanh(x).^(-1);
    b = x.^(-1);
    
    fx = a - b - D;
    f1x = 1 - a^2 + b^2;
    
    x = x - fx/f1x;
    
    x_v(ii) = x;
    
    diff = abs(x_v(ii-1) - x);
    
    ii = ii + 1;
end

end