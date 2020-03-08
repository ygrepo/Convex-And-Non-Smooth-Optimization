function [f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit)
% code for gradient method, including backtracking line search, goes here,
% instead of the following dummy code, which just takes one gradient step
% and does not even check whether fun(x) < fun(x0).
% Note that the input parameter 'fun' is an anonymous function.
%
%[f0,g0] = fun(x0);  % call the anonymous function to get f0 and g0
%x = x0 - g0;        % a gradient step
%f_all(1) = fun(x);         % evaluate fun at new x
%gnorm_all(1) = norm(g0);   % we took one iteration with no line search


%iteration counter
itc = 1;
xc = x0;
[fc,gc] = fun(xc);
alpha = 1./4;
beta = 0.5;

ithist=zeros(maxit,5);
ithist(1,1)=fc;
ithist(1,2) = norm(gc);
ithist(1,3)=0;
ithist(1,4)=1;
ithist(1,5)=itc-1;

while (norm(gc) > tol) && (itc <maxit)
    %
    %  the original step t=1
    t = 1;
    
    xt = xc - t * gc;
    ft = fun(xt);
    fgoal = fc - alpha * t * (gc'* gc);
    j = 0;
    %  backtrack line search
    while ft > fgoal
        t = beta * t;
        xt = xc - t * gc;
        ft = fun(xt);
        fgoal = fc - alpha * t * (gc'* gc);
        j = j + 1;
    end
    xc=xt;
    [fc,gc]=fun(xc);
    itc = itc+1;
    ithist(itc,1)=fc;
    ithist(itc,2) = norm(gc);
    ithist(itc,3)=j;
    ithist(itc,4)=t;
    ithist(itc,5)=itc-1;
end
f_all = ithist(:,1);
gnorm_all = ithist(:,2);