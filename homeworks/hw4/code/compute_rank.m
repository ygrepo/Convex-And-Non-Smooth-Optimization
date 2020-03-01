function rk = compute_rank(A)
% Assume A is symmetric
D = eig(A);
[n, ~] = size(A);
k = 0;
for i=1:n
    if D(i,1) == 0
        k = k + 1;
    end
end   
rk = n - k;
end