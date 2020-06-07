function [f_all, gnorm_all, x_all] = nesterov(fun, x0, M, q, tol, maxit)
% Code for Nesterov algorithm
%Input
%   fun : anonymous function x0: starting point
%   x0: starting point
%   M: maximum eigenvalue
%   q: Nesterov parameter
%   tol: tolerance on the norm of the gradient 
%   maxit: max number of iterations
%Output
%   f all : vector of function values f(x(k))
%   gnorm all : corresponding vector of gradient norms
%   x_all: position points x

% Initialization
k = 1; % iteration index
f_all = [];
gnorm_all = [];
x_all = {};

xk = x0;
yk = x0;


% Main loop of the gradienst descent
while k  <= maxit % controls the number of iterations max.
    [f, g, ~] = fun(yk);
    
    % Record yk and its value
    f_all = [f_all; f];
    gnorm_all = [gnorm_all; norm(g)];
    if gnorm_all(k) < tol % check tolerance on the norm of the gradient
        break
    end
    x_all{end + 1} = yk;
    
    % Take the step for x and y
    xkp1 = yk - 1/M * g;
    ykp1 = xkp1 + q * (xkp1 - xk);
    xk = xkp1;
    yk = ykp1;
    
    k = k + 1;    
end

end