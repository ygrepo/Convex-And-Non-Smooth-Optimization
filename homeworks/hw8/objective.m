function p = objective(A, b, lambda, x)
% objective function
    p = 0.5*norm(A*x - b, 2)^2 + lambda*norm(x,1);
end