%% 1.(b)
clear all;
close all;
clc;
load Adata.mat

% Initial size
% n = 100
% A = load("Adata.mat").A;

[n m] = size(A);
fun = @(x)logfunct(x, A);

% now we call the gradient method
x0 = zeros(n,1);
tol = 1e-6;
maxit = 100;

[f_all,gnorm_all, x_end] = gradmeth(fun, x0, tol, maxit);

pstar = f_all(end);

estimate(fun, x_end, A)

% plot the function values
tiledlayout(1,2);
nexttile
semilogy(f_all - pstar);
xlim([0 80]);
xlabel("k");
ylabel("Log value")
title('Log value of $f(x^{k}) -p^*$','Interpreter','latex');
legend('$f(x^{k}) -p^*$','Interpreter','latex');

% plot the gradient
nexttile
semilogy(gnorm_all);
xlim([0 80]);
xlabel("k");
ylabel("Log value")
title("Log value of $\Vert \nabla f \Vert k$","Interpreter","latex"); 
legend("$\Vert \nabla f \Vert k$","Interpreter","latex"); 
