function [t] = mat_solver(X,B)
[n,~] = size(X);
b = eig(B);
t = inf;
sum(b>=0)
if sum(b>=0) ~= n
    L = chol(X);
    I=eye(n);
    L=I/L;
    B=-B;
    t =max(eig(L'*B*L));
    t = 1 / t;
end
end