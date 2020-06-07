%% HW 6 − BV 9.30 (p.519)

A = load('Adata.mat').A;
n = size(A,1);

% The objective function of BV 9.30 is implemented in bv930 fun.m
fun = @(x)logfunct(x, A);


% starting point and parameters
x0 = zeros(n,1);
tol = 1e-6;
maxit = 100;

% Condition number of A
[f, g, h] = fun(x0);
M = max(eig(h));
m = min(eig(h));
% Launch grad descent with t=1/M to re−estimate bounds from sequence
[f_all, gnorm_all, x_all] = gradmeth(fun, x0, 1/M, tol, maxit);
[M_max, m_min, L_max] = estimate(fun, x_all{end}, A);

% min_eigs = [];
% max_eigs = [];
% for i=1:length(x_all)
%     h = hessian(A, x_all{i});
%     min_eigs = [min_eigs, min(eig(h))];
%     max_eigs = [max_eigs, min(eig(h))];
% end    
M = M_max;
m = m_min;
kappa = M / m;

% Define optimal value to be last element of Newton sequence
p_star = f_all(end);
x_star = x_all(end);

% Nesterov's q parameter
q = (1 - sqrt(1/kappa)) /  (1 + sqrt(1/kappa));

% Call Nesterov algorithm
[f_all_nest, gnorm_all_nest,~] = nesterov(fun, x0, M, q, tol, maxit);

% Estimate optimal value by last element of Nesterov sequence
p_star = f_all_nest(end);

% x-axis values
t = linspace(1,maxit, maxit);

% Verify the lower bound of Nesterov
figure
semilogy(f_all_nest - pstar, "r-*", 'LineWidth',1,'MarkerSize',3);
hold on
% plot lower bound
lower_bd = (m/2) * q.^(2*t) * (norm(x0 - xstar)^2);
semilogy(lower_bd, "k-", 'LineWidth',1);

grid on
title('Objective error and lower bound for Nesterov method',...
'FontSize', 14);
xlabel('Iterations');
ylabel('Objective value');
legend('$f(x^{k}) -p^*$', ...
'Lower bound', ...
'FontSize', 14, 'Location', 'NorthEast','Interpreter','latex');

%% Gradient descent with fixed step size t = 1/M
% Step size
t1 = 1/M;

% Launch the gradient descent implemented in gradmeth.m
[f_all_grad1, gnorm_all_grad1, ~] = gradmeth(fun, x0, t1, tol, maxit);

% Verify convergence
% Plot objective error ratio and upper bound in a LOG plot
figure
% Plot error ratio
err_ratio_grad1 = ( f_all_grad1 - pstar ) / (f_all_grad1(1) - pstar);
semilogy(err_ratio_grad1, "b-*", 'LineWidth',1,'MarkerSize',3);
hold on

% Plot upper bound
upper_bd_grad1 = (1-1/kappa).^t;
semilogy(upper_bd_grad1, "k-", 'LineWidth',1);

% Custom
grid on
title('Objective error ratio and upper bound for gradient method 1 ( t=1/M )',...
'FontSize', 14);
xlabel('Iterations');
ylabel('Objective value');
legend('$\frac{f(x^{k}) -p^*} {f(x^{0}) -p^*}$', ...
'Upper bound', ...
'FontSize', 14, 'Location', 'NorthEast','Interpreter','latex');

%% Gradient descent with fixed step size t = 2/(m + M)
% Step size
t2 = 2/(m + M);

% Launch the gradient descent implemented in gradmeth.m
[f_all_grad2, gnorm_all_grad2, ~] = gradmeth(fun, x0, t2, tol, maxit);

% Verify convergence
% Plot objective error ratio and upper bound in a LOG plot
figure
% Plot error ratio
err_ratio_grad2 = ( f_all_grad2 - pstar ) / (f_all_grad2(1) - pstar);
semilogy(err_ratio_grad2, "b-*", 'LineWidth',1,'MarkerSize',3);
hold on

% Plot upper bound
upper_bd_grad2 = kappa *( (1  - 1/kappa) / (1+1/kappa) ).^t;
semilogy(upper_bd_grad2, "k-", 'LineWidth',1);

% Custom
grid on
title('Objective error ratio and upper bound for gradient method 2 ( t=2/(m +M) )',...
'FontSize', 14);
xlabel('Iterations');
ylabel('Objective value');
legend('$\frac{f(x^{k}) -p^*} {f(x^{0}) -p^*}$', ...
'Upper bound', ...
'FontSize', 14, 'Location', 'NorthEast','Interpreter','latex');


%% Plot objective for 3 methods in same plot
% Objective
figure, hold on, grid on
% Nesterov
plot(t, f_all_nest, "b-*", 'LineWidth',1,'MarkerSize',3); % Nesterov
plot(t, f_all_grad1, "b-x", 'LineWidth',1,'MarkerSize',3); % Gradient 1
plot(t, f_all_grad2, "c-o", 'LineWidth',1,'MarkerSize',3); %  Gradient 2
% Legend, axis and title
title('Objective decrease for the 3 methods','FontSize', 14);
xlabel('Iterations'); 
ylabel('Objective value'); 
legend('f(xˆ{(k)}) (Nesterov)', ...
'f(xˆ{(k)}) (Gradient, t=1/M)', ... 
'f(xˆ{(k)}) (Gradient, t=2/(M+m)', ... 
'pˆ{*} (Optimal)', ...
'FontSize', 14, 'Location', 'NorthEast');


%% Plot objective error for 3 methods in same plot
% Objective
figure, hold on, grid on
% Nesterov
semilogy(t, f_all_nest - pstar, "b-*", 'LineWidth',1,'MarkerSize',3); % Nesterov
semilogy(t, f_all_grad1  - pstar, "b-x", 'LineWidth',1,'MarkerSize',3); % Gradient 1
semilogy(t, f_all_grad2  - pstar, "c-o", 'LineWidth',1,'MarkerSize',3); %  Gradient 2
% Legend, axis and title
title('Objective error for the 3 methods','FontSize', 14);
xlabel('Iterations'); 
ylabel('Objective value'); 
legend('f(xˆ{(k)}) (Nesterov)', ...
'f(xˆ{(k)}) (Gradient, t=1/M)', ... 
'f(xˆ{(k)}) (Gradient, t=2/(M+m)', ... 
'pˆ{*} (Optimal)', ...
'FontSize', 14, 'Location', 'NorthEast');


%% Plot gradient norm for 3 methods in same LOG plot
% Objective
figure, hold on, grid on
% Nesterov
semilogy(t, gnorm_all_nest, "b-*", 'LineWidth',1,'MarkerSize',3); % Nesterov
semilogy(t, gnorm_all_grad1, "b-x", 'LineWidth',1,'MarkerSize',3); % Gradient 1
semilogy(t, gnorm_all_grad2, "c-o", 'LineWidth',1,'MarkerSize',3); %  Gradient 2
% Legend, axis and title
title('Gradient norm for the 3 methods','FontSize', 14);
xlabel('Iterations'); 
ylabel('Gradient norm'); 
legend('f(xˆ{(k)}) (Nesterov)', ...
'f(xˆ{(k)}) (Gradient, t=1/M)', ... 
'f(xˆ{(k)}) (Gradient, t=2/(M+m)', ... 
'pˆ{*} (Optimal)', ...
'FontSize', 14, 'Location', 'NorthEast');
