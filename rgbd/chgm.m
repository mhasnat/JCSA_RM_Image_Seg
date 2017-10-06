function hg=chgm(a,b,x)

%     M(a,b,x)
%     Input  : a  --- Parameter
%     b  --- Parameter ( b <> 0,-1,-2,... )
%     x  --- Argument
%     Output:  HG --- M(a,b,x)
%     Routine called: GAMMA for computing ג(x)
%     ===================================================



ta=[];tb=[];xg=[];tba=[];
pi=3.141592653589793d0;
a0=a;
a1=a;
x0=x;
hg=0.0d0;
if (b == 0.0d0|b == -abs(fix(b))) ;
    hg=1.0d+300;
elseif (a == 0.0d0|x == 0.0d0);
    hg=1.0d0;
elseif (a == -1.0d0);
    hg=1.0d0-x./b;
elseif (a == b);
    hg=exp(x);
elseif (a-b == 1.0d0);
    hg=(1.0d0+x./b).*exp(x);
elseif (a == 1.0d0&b == 2.0d0);
    hg=(exp(x)-1.0d0)./x;
elseif (a == fix(a)&a < 0.0d0);
    m=fix(-a);
    r=1.0d0;
    hg=1.0d0;
    for  k=1:m;
        r=r.*(a+k-1.0d0)./k./(b+k-1.0d0).*x;
        hg=hg+r;
    end;
end;
if (hg ~= 0.0d0) return; end;
if (x < 0.0d0) ;
    a=b-a;
    a0=a;
    x=abs(x);
end;
if (a < 2.0d0) nl=0; end;
if (a >= 2.0d0) ;
    nl=1;
    la=fix(a);
    a=a-la-1.0d0;
end;
for  n=0:nl;
    if (a0 >= 2.0d0) a=a+1.0d0; end;
    if (x <= 30.0d0+abs(b)|a < 0.0d0) ;
        hg=1.0d0;
        rg=1.0d0;
        for  j=1:500;
            rg=rg.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*x;
            hg=hg+rg;
            if (abs(rg./hg) < 1.0d-15) break; end;
        end;
    else;
        [a,ta]=gamma(a,ta);
        [b,tb]=gamma(b,tb);
        xg=b-a;
        [xg,tba]=gamma(xg,tba);
        sum1=1.0d0;
        sum2=1.0d0;
        r1=1.0d0;
        r2=1.0d0;
        for  i=1:8;
            r1=-r1.*(a+i-1.0d0).*(a-b+i)./(x.*i);
            r2=-r2.*(b-a+i-1.0d0).*(a-i)./(x.*i);
            sum1=sum1+r1;
            sum2=sum2+r2;
        end;
        hg1=tb./tba.*x.^(-a).*cos(pi.*a).*sum1;
        hg2=tb./ta.*exp(x).*x.^(a-b).*sum2;
        hg=hg1+hg2;
    end;
    if (n == 0) y0=hg; end;
    if (n == 1) y1=hg; end;
end;
if (a0 >= 2.0d0) ;
    for  i=1:la-1;
        hg=((2.0d0.*a-b+x).*y1+(b-a).*y0)./a;
        y0=y1;
        y1=hg;
        a=a+1.0d0;
    end;
end;
if (x0 < 0.0d0) hg=hg.*exp(x0); end;
a=a1;
x=x0;
return;



function [x,ga]=gamma(x,ga);

%     Input :  x  --- Argument of ג(x)
%     ( x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================


%
%
%
g=zeros(26,1);

pi=3.141592653589793d0;
if (x == fix(x)) ;
    if (x > 0.0d0) ;
        ga=1.0d0;
        m1=x-1;
        for  k=2:m1;
            ga=ga.*k;
        end;
    else;
        ga=1.0d+300;
    end;
else;
    if (abs(x) > 1.0d0) ;
        z=abs(x);
        m=fix(z);
        r=1.0d0;
        for  k=1:m;
            r=r.*(z-k);
        end;
        z=z-m;
    else;
        z=x;
    end;
    g=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
    gr=g(26);
    for  k=25:-1:1;
        gr=gr.*z+g(k);
    end;
    ga=1.0d0./(gr.*z);
    if (abs(x) > 1.0d0) ;
        ga=ga.*r;
        if (x < 0.0d0) ga=-pi./(x.*ga.*sin(pi.*x)); end;
    end;
end;
return;