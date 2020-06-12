function [x_star, p_star, f_all, x_all, z_all, u_all, r_norm, s_norm, ...
    eps_pri, eps_dual, max_iter] = lasso_admm(A, b, lambda, rho)

% Description
% Solve the problem:
% Minimize (1/2) || A x -b ||_2^2 + \lambda ||x||_1
% using the ADMM method.
% Input:
% A: M x N matrix
% b: N x 1 matrix
% lambda: real number > 0.
% rho: real > 0, penalty parameter.
% Output:
% f_all: objective values, minimum to Lasso problem at iteration k.
% x_all, z_all, u_all: variables of ADMM algorithm.
% r_norm, s_norm: primary and dual residual norms.
% x_star: computed argmin to the Lasso problem.
% p_star: computed argmin to the lasso problem of objective function.

% Global constants and defaults

% MAX_ITER = 100;
% ABSTOL   = 1e-4;
% RELTOL   = 1e-2;

MAX_ITER = 1000;
ABSTOL   = 1e-7;
RELTOL   = 1e-4;

% cached computations
AtA = A'*A;
Atb = A'*b;


[m, n] = size(A);


x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

[L, U] = factor(A, rho);

f_all = [];
x_all = {};
z_all = {};
u_all = {};
r_norm = [];
s_norm = [];
eps_pri = [];
eps_dual = [];
max_iter = 0;

for k = 1:MAX_ITER
    
    % x-update
    q = Atb + rho*(z - u);
    if m >= n
        x = U \ (L \ q);
    else
        x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
    end
    x_all{end + 1} = x;
    
    % z-update
    zold = z;
    z = prox_l1(x + u, lambda/rho);
    z_all{end + 1} = z;
    
    % u-update
    u = u + x - z;
    u_all{end + 1} = u;
    
    % diagnostics, reporting, termination checks
    p_x = 0.5*norm(A*x - b, 2)^2 + lambda*norm(x,1);
    p_z = 0.5*norm(A*z - b, 2)^2 + lambda*norm(z,1);
    [p_star, idx ] = min([p_x p_z]);
    x_star = (idx == 1) * x + (idx==2) * z;
    f_all(k) = p_star;
    
    % Stopping criteria p19
    r_norm(k)   = norm(x - z);
    s_norm(k)   = norm(-rho*(z - zold));
    eps_pri(k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    eps_dual(k) = sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
    max_iter = k;
    
    if r_norm(k) < eps_pri(k) && s_norm(k) < eps_dual(k)
        break;
    end
    
end

p_star = f_all(end);

end
