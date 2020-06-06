function [f_all,gnorm_all, x_end] = netwmeth(fun, x0, tol, maxit)
%code for gradient method, including backtracking line search
%Input
%   fun : anonymous function x0: starting point
%   tol: tolerance on the norm of the gradient maxit: max number of iterations
%Output
%   f all : vector of function values f(x(k))
%   gnorm all : corresponding vector of gradient norms

% Initialization
i = 1; % iteration index
f_all = [];
gnorm_all = [];

x = x0;
x_end = x0;

while i <= maxit % controls the number of iterations max.
    [f, g, h] = fun(x);
    dx = -h\g; % Use negative newtown update as descent direction
    f_all(i) = f;
    gnorm_all(i) = norm(g);
    
    % Backtracking line search
    t = BTLS(fun, x,f,g,dx); % Determine proper step size
    x = x + t.*dx; % update x
    x_end = x;
    
    if gnorm_all(i) < tol % check tolerance on the norm of the gradient
        break
    end
    
    
    i = i + 1;
end
end



function t = BTLS(fun, x,f,g,dx)
% Backtracking line search
% Input:
%   fun: anonymous function
%   x: current position
%   f: function value
%   g: gradient at x
%   dx: descent direection
%
% Output
% t: step size in gradient descent
%
% Parameters

alpha = 0.25;
beta = 0.5;

t = 1; % Initial step size
while 1
    f_new = fun(x + t.*dx);
    if f_new > f + alpha.* t.*g'*dx
        t = beta*t;
    else
        break;
    end
end
end

