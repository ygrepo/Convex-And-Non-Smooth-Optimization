%% hw9 - 3

clear;
clc;

% For each of X1, X2, X3 we will perform N trials to average the results
% Each trial = run successive nuclear norm minimization SDPs, increasing
% the number of constraints each time. We stop the SDPs when Xi has been
% approximately recovered, i.e. when solution X of SDP is close enough
% in a sense of the Froebenius norm (parameter 'tol')

% Froebenius norm tolerance: did we recover X1 accurately enough?
tol = 1e-6;

%% Load data
X1 = load('Xdata.mat').X1;

%% Perform a single
% norm error across SDPs in a log plot
[~, ~, fnorms] = low_rank_completion(X1, tol);
figure
semilogy(fnorms,"b-o","LineWidth", 2);
grid on
title("Froebenius norm error of the SDP solution vs number of constraints","FontSize", 14);
xlabel('Number of specified entries');
ylabel("Froebenius norm");
legend('$X-X_1$','FontSize', 14, "Location", "NorthEast");

%% Analysis of repeated trials on X1 ,X2 and X3.
% Number of trials for averaging
N = 10;

% Problem 3 - X1
% To store the number of constraints required to converge for each trial
n_ctrs_1 = zeros(N,1);
% Store the 'final' recovered X (i.e. within 'tol' from X1)
X_final_1 = {};
% Nuclear norm of the final X (to be compared with nuclear norm of X1)
norm_final_1 = zeros(N,1);

for i=1:N
    [Xs, ids, fnorms] = low_rank_completion(X1, tol);
    X_final_1{i} = Xs{end};
    norm_final_1(i) = fnorms(end);
    n_ctrs_1(i) = length(Xs);
end    

% Average number of constraints required to recover X1 within tolerance
avg_1 = sum(n_ctrs_1)/N;

% Plot number of iterations
figure
grid on
plot(n_ctrs_1,"b-o","LineWidth", 2);
hold on
plot(avg_1 * ones(N,1), "r--","LineWidth", 2);  
title("Average number of constraints required to recover X1 through nuclear norm minimization",...
    "FontSize", 14)
xlabel("Trials")
ylabel("Number of constraints")
legend("Measured",....
       "Average",...
       "FontSize", 14, "Location", "NorthEast");
saveas(gcf,"x1_average_nb_constraints",'pdf')

%% Problem 4

% Generate random matrix of size q and rank 3 (X_low)
% Matrix size 
s = 15;
% Matrix desired rank
rnk = 3;

% Random matrix of rank 3
% Note that rank can be checked by counting the non−zero entries of svd(X)
tmp = randn(s, rnk);

% Create square matrix of rank rnk
X_high = tmp * tmp';

% norm error across SDPs in a log plot
[~, ~, fnorms] = low_rank_completion(X_high, tol);
figure
semilogy(fnorms,"b-o","LineWidth", 2);
grid on
title("Froebenius norm error of the SDP solution vs number of constraints","FontSize", 14);
xlabel('Number of specified entries');
ylabel("Froebenius norm");
legend('$X-X_low$','FontSize', 14, "Location", "NorthEast");


% Random matrix of rank higher 10 (X_high)

% Matrix size 
s = 15;
% Matrix desired rank
rnk = 10;

tmp = randn(s, rnk);

% Create square matrix of rank rnk
X_high = tmp * tmp';

% norm error across SDPs in a log plot
[~, ~, fnorms] = low_rank_completion(X_high, tol);
figure
semilogy(fnorms,"b-o","LineWidth", 2);
grid on
title("Froebenius norm error of the SDP solution vs number of constraints","FontSize", 14);
xlabel('Number of specified entries');
ylabel("Froebenius norm");
legend('$X-X_high$','FontSize', 14, "Location", "NorthEast");
