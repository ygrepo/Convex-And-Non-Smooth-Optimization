%% experiment

clear all
clc



M = 100000;
density = 2/M;
N = 10000;
A = sprandn(M, N, density);
b=randn(M,1);

lambda = 0.2;
rho = 2;

h = lasso(A, b, lambda, rho);

h.admm_iter = length(h.admm_optval);
K = max([h.admm_iter]);
h.cvx_optval  = h.p_cvx*ones(K,1);
h.admm_optval = padarray(h.admm_optval', K-h.admm_iter, h.p_admm, 'post');
fig = figure;

plot(1:K, h.cvx_optval,  'k--', ...
     1:K, h.admm_optval, 'r-');
 %    1:K, h.prox_optval, 'r-', ...
  %   1:K, h.fast_optval, 'g-', ...

%xlim([0 45]);
legend('True', 'ADMM');
print -depsc lasso_comp.eps;
