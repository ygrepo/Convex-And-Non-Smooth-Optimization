%% TestLASSO
% In this script you will be implementing LASSO regression using CVX

clear all
clc
load LASSO_data.mat

[n,p]=size(X);

n_lambda = length(lambdas) - 1;
B=zeros(p,n_lambda);


%% YOUR CVX CODE HERE!
% Write a loop, and save the solution for each lambda to one column of B.

for i=1:20
%INPUT: X,y, lambdas
cvx_begin quiet
    cvx_precision low
    variable x(p)
    minimize(0.5*sum_square(y - X * x) + lambdas(i) * norm(x,1))
cvx_end
B(:,i) = x;
end


%OUTPUT: B


%% Script for them to plot results

% Plot the solution path against tuning parameter lambda
figure(1)
plot(lambdas(1:20),B(:,1:20)');

% % Another way to visualize the solution path (for anyone interested)
% xx=sum(abs(B));
% xx=xx/xx(end);
% plot(xx,B','-');
xlabel('lambda')
ylabel('coefficients')

figure(2)
mse=bsxfun(@minus,B,b_gnd);

% Plot the MSE of recovery
plot(sum(mse.*mse)/p,'*-');
grid on
ylim([0.0008,0.0015])
xlabel('index of lambda')
ylabel('MSE')
