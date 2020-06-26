%% hw10 - 1
clear all;
clc;

x = randi(10,1,2)
b = randi(10,1,2)
t = vec_solver(x,b);
fprintf('t = %.2f\n', t);

%% hw10 - 2
n = 3;
X = rand(n,n);
X = X * X';
X = X + eye(n);
B = [1 2 3; 2 2 1; 3 1 3]
% B = rand(n,n);
% B = B * B';
t = mat_solver(X,B);
t