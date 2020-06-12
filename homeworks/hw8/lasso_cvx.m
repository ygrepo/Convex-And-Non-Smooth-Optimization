function [x_star, p_star] = lasso_cvx(A,b, lambda)
% Description
% Solve the problem:
% Minimize (1/2) || A x -b ||_2^2 + \lambda ||x||_1
% using  CVX.
% Input:
% A: M x N matrix
% b: N x 1 matrix
% lambda: real number > 0.
% Output:
% x_star: computed argmin to the Lasso problem.
% p_star: computed argmin to the lasso problem of objective function.

[~,n] = size(A);
cvx_begin quiet
    variable x(n)
    minimize(0.5*square_pos( norm(A*x - b, 2) ) + lambda *norm(x,1))
cvx_end
x_star = x; % set optimal value x
p_star = cvx_optval; % set value of objective function
