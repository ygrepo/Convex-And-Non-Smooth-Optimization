function [f_all, gnorm_all, x_all] = gradmeth(fun, x0, t, tol, maxit)

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
x_all = {};

x = x0;

% Main loop of the gradient descent
while i <= maxit % controls the number of iterations max.
    [f, g, ~] = fun(x);
    dx = -g; % Use negative gradient descent as descent direction
    f_all(i) = f;
    gnorm_all(i) = norm(g);
    if gnorm_all(i) < tol % check tolerance on the norm of the gradient
        break
    end
    x_all{end + 1} = x;

    % Backtracking line search
    %t = BTLS(fun, x,f,g,dx); % Determine proper step size
    x = x + t*dx; % update x
    
    i = i + 1;
end
end

