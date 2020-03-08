function [f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit)
% code for netwon method, including backtracking line search

%iteration counter
itc = 1;
xc = x0;
[fc, gc, hc] = fun(xc);
alpha = 1./4;
beta = 0.5;

ithist=zeros(maxit,5);
ithist(1,1) = fc;
ithist(1,2) = norm(gc);
ithist(1,3) = 0;
ithist(1,4) = 1;
ithist(1,5) = itc-1;

while (norm(gc) > tol) && (itc <maxit)
    %
    %  the original step t=1
    t = 1;
    
    dx = -hc\gc;
    xt = xc + t * dx;
    ft = fun(xt);
    fgoal = fc + alpha * t * (gc'* dx);
    j = 0;
    %  backtrack line search
    while ft > fgoal
        t = beta * t;
        xt = xc + t * dx;
        ft = fun(xt);
        fgoal = fc + alpha * t * (gc'* dx);
        j = j + 1;
    end
    xc=xt;
    [fc, gc]=fun(xc);
    itc = itc+1;
    ithist(itc,1)=fc;
    ithist(itc,2) = norm(gc);
    ithist(itc,3)=j;
    ithist(itc,4)=t;
    ithist(itc,5)=itc-1;
end
f_all = ithist(:,1);
gnorm_all = ithist(:,2);