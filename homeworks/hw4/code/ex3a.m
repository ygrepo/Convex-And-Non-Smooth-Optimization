% Solve the SDP given in BV (5.114) using data set 1 and data set 2
% In addition solve BV (5.115) using data set 1:
% 1- Compare computed optimal and dual variables from (5.115)
%    with computed optimal primal variables from (5.114) and vice-versa
% 2- Compute approximate ranks of the computed optimal primal and dual
%    matrices
% 3- Test for complimentarity of these matrices.


clear
cvx_quiet true

W = load("hw4data1.mat").W;
[opt_val, X1, nu1] = solve_dual(W);
fprintf("Found for data set 1, Lagrange dual problem, optimal value:%4.2f\n",  opt_val)

[opt_val, X2, nu2, lambda] = solve_sdp_relaxation(W);
fprintf("Found for data set 1, SDP relaxation of two-way partitioning problem,optimal value:%4.2f\n",  opt_val)

if is_equal(abs(nu1), abs(nu2))
    disp("Computed dual variables from 5.115 same as computed primal variables from 5.114")
else
    disp("Computed dual variables from 5.115 different from the computed primal variables from 5.114")
end

if is_equal(X1, X2)
    disp("Computed primal variables from 5.115 same as computed dual variables from 5.114")
else
    disp("Computed primal variables from 5.115 different than computed dual variables from 5.114")
end

fprintf("Rank(optimal primal matrix): %d, Rank(optimal dual matrix): %d\n",...
    compute_rank(X2) , compute_rank(lambda))

tol = 1e-4;
if are_complementary(X2, lambda, tol)
    fprintf("Computed optimal and dual matrices are complementary, tol:%2.5f\n", tol)
else    
    fprintf("Computed optimal and dual matrices are not complementary, tol:%5.2f\n", tol)
end

clear
W = load("hw4data2.mat").W;
[opt_val, ~, ~] = solve_dual(W);
fprintf("Found for data set 2, Lagrange dual problem, optimal value:%4.2f\n",  opt_val)
[opt_val, ~, ~, ~] = solve_sdp_relaxation(W);
fprintf("Found for data set 2, SDP relaxation of two-way partitioning problem,optimal value:%4.2f\n",  opt_val)

function [opt_val, X, nu] = solve_dual(W)

[n,~]=size(W);

cvx_begin sdp
    variable nu(n)
    dual variable X
    maximize (-sum(nu, 1))
    X: W + diag(nu) == semidefinite(n);
cvx_end
opt_val = cvx_optval;
end

function [opt_val, X, nu, lambda]  = solve_sdp_relaxation(W)

[n,~]=size(W);

cvx_begin sdp
    variable X(n, n) symmetric
    dual variables nu lambda 
    minimize trace(W * X)
    nu: diag(X) == 1;
    lambda: X == semidefinite(n);
cvx_end
opt_val = cvx_optval;

end

function t = is_equal(A, B)
% We are expecting that A and B have same dimensions.

    [m, n] = size(A);
    if (sum(A == B, "all") == m * n)
        t = true;
    else
        t = false;
    end
end

function rk = compute_rank(A)
% Assume A is symmetric
D = eig(A);
[n, ~] = size(A);
k = 0;
for i=1:n
    if D(i,1) == 0
        k = k + 1;
    end
end   
rk = n - k;
end

function t = are_complementary(A, B, tol)
% Testing if A * B ~ 0
    t = isempty(find(A* B > tol));
end