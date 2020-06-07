function [f, g, x] = nestex_fun(x, m, M)
% Objective function of Nesterov's worst case example (p.67)
% This function returns the evaluation, gradient and hessian and x

% Sizes
n = size(x,1);

% Value of objective
f = ( (M-m)/8 )*( x(1)^2 - 2*x(1) ) + (m/2)*norm(x,2)^2;

for i=1:n-1
    f = f + ( (M-m)/8 ) * ( x(i) - x(i+1) )^2;
end
T= zeros(n,n);
T(1:1+n:n*n) = 2;
T(n+1:1+n:n*n) = -1;
T(2:1+n:n*n) = -1;
T = sparse(T);

e = zeros(n,1);
e(1) = 1;
g = sparse((((M-m)/4)*T + m*eye(n))*x - ((M-m)/4)*e);
end

