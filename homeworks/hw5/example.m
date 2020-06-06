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
[f_all,gnorm_all,~] = gradmeth(fun, x0, tol, maxit);

fxl = f_all(end);
xopt = -A\b;
popt = fun(xopt);
alpha = 1./4;
beta = 0.5;
eigv = eig(A);
m = eigv(1);
M = eigv(end);
tv = (1 - 2 * m * alpha * min(1, beta/M));
fprintf("After %d iterations, ratio:%8.6f, should be smaller than theoretical value:%8.6f\n",...
    maxit, (fxl - popt) / (fun(x0) - popt), tv)
