%% 2.(b)
clear all;
close all;
clc;
load Adata.mat

% Initial size
[n m] = size(A);
fun = @(x)logfunct(x, A);

% now we call the gradient method
x0 = zeros(n,1); % starting point
tol = 1e-8;
maxit = 100;

[f_all,gnorm_all, ~] = gradmeth(fun, x0, tol, maxit);
[f2_all,gnorm2_all, x_end] = netwmeth(fun, x0, tol, maxit);

pstar = f2_all(end); % optimal value

estimate(fun, x_end, A)
 
% plot the function values
tiledlayout(1,2);
nexttile
semilogy(f_all - pstar);
hold on
semilogy(f2_all - pstar);
xlim([0 72]);
xlabel("k");
ylabel("Log value")
title('Log value of $f(x^{k}) -p^*$','Interpreter','latex');
legend('Gradient','Newtown', 'Interpreter', 'latex');

% plot the gradient
nexttile
semilogy(gnorm_all);
hold on
semilogy(gnorm2_all);
xlim([0 72]);
xlabel("k");
ylabel("Log value")
title("Log value of $\Vert \nabla f \Vert k$","Interpreter","latex"); 
legend('Gradient','Newtown', 'Interpreter', 'latex');
