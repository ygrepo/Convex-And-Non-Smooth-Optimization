% this is a script showing an example of how to call quad.m via gradmeth.m
% first we need to define A and b which will be used in quad.m
n = 5;
A = hilb(n); % the Hilbert matrix, which is a very ill conditioned matrix
b = ones(n,1);
% now we define the anonymous function.  Note that the first x indicates
% that we are going to be calling the function repeatedly with different x
% but with A and b fixed
fun = @(x)quad(x, A, b);
% now we call the gradient method, which will just take one gradient step 
x0 = ones(n,1);
tol = 1e-6;
maxit = 100;
xopt = -A\b;
popt = fun(xopt);
alpha = 1./4;
beta = 0.5;
eigv = eig(A);
m = eigv(1);
M = eigv(end);
k = m / M;

t = 1 /M;
upper_bound = (1 - k)^maxit;
[f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit, t);
fxl = f_all(end);
initial_diff = fun(x0) - popt;
last_diff = fxl - popt;
fprintf("After %d iterations, difference:%8.6f, initial_diff:%8.6f, bound:%8.6f\n",...
    maxit, last_diff, initial_diff, upper_bound * initial_diff)

t = 2 /(m + M);
kinv = 1/k;
maxit = 900;
upper_bound = (1-kinv)/(1+kinv);
upper_bound = upper_bound^(2 * maxit);
upper_bound = k * upper_bound;
[f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit, t);
fxl = f_all(end);
initial_diff = fun(x0) - popt;
last_diff = fxl - popt;
fprintf("After %d iterations, difference:%8.6f, initial_diff:%8.6f, bound:%8.6f\n",...
    maxit, last_diff, initial_diff, upper_bound * initial_diff)

