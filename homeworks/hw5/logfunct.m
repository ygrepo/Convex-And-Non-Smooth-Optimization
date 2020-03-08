function [f,g, hess] = logfunct(x, A)
[~,m] = size(A);
f = -sum(log(1 - A * x)) - sum(log(1+x)) - sum(log(1-x));
d = 1./(1 - A * x);
g = A' * d - 1./(1+x) + 1./(1-x);
hess =  A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
if all(abs(A * x) < 1) && all(abs(x) < 1)
    %disp("less than 1")
else
    %disp("some element not in domain of f")
    f = inf;
    g = nan * ones(m,1);
    hess = nan * ones(m,m);
end
end