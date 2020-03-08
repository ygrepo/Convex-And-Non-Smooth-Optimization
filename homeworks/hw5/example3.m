A = load("Adata.mat").A;
[n,m] = size(A);
fun = @(x)logfunct(x, A);
% now we call the gradient method
x0 = zeros(m,1);
tol = 1e-6;
maxit = 100;
[f_all,gnorm_all] = newtmeth(fun, x0, tol, maxit);
optbval = f_all(end);
 
figure(1)
xvalues = 1:maxit;
semilogy(xvalues, f_all - optbval, "-");
xlabel("k");  
ylabel(['$ f(x^({k})) -p^* $'],'interpreter','latex')
cond(A)
