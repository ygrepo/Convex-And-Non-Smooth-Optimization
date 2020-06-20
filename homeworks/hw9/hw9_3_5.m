%% Problem 5
clear all;
clc;


% This time we perform nuclear norm minimization using ADMM

% Generate random matrix of size q and rank 3 (X_low)
% Matrix size 
s = 10;
% Matrix desired rank
rnk = 1;

% Random matrix of rank 1
% Note that rank can be checked by counting the non−zero entries of svd(X)
tmp = randn(s, rnk);

% Create square matrix of rank rnk
XDATA = tmp * tmp';

% Tolerance on recovery error (Froeb norm)
tol = 1e-4;

% norm error across SDPs in a log plot
[Xs, ids, fnorms] = run_admm_sdps(XDATA, tol);

figure
semilogy(fnorms,"b-o","LineWidth", 2);
grid on
title("Froebenius norm error of the SDP solution vs number of constraints","FontSize", 14);
xlabel('Number of specified entries');
ylabel("Froebenius norm");
legend('$X-X$','FontSize', 14, "Location", "NorthEast");
