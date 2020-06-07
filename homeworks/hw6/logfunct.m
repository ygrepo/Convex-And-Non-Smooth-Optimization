function [f, g, h] = logfunct(x, A)
% objective function for 9.30

[n,m] = size(A);

if all(max(A' * x) < 1) && all(norm(x, inf) < 1) % x in dom f
    % function value
    f = -sum(log(1 - A' * x)) - sum(log(1+x)) - sum(log(1-x));
    % gradient
    g = A * (1./(1 - A' * x)) - 1./(1+x) + 1./(1-x);
    % hessian
    h = zeros(n, n); 
    for i=1:m
        h = h + (A(:, i) * A(:, i)') / (1-A(:, i)' *x).^2; 
    end
    h = h + 2* diag( (ones(n,1)+diag(x) * x) ./ (ones(n,1)- diag(x)* x).^2);
else % x not in dom f
    %disp("some element not in domain of f")
    % function value
    f = inf;
     % gradient
    g = nan * ones(m,1);
    % hessian
    h = nan * zeros(n, n);
end

end