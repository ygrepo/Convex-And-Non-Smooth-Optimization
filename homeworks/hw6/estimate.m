function [M_max, m_min, L_max] = estimate(fun, xend, A)
% Input
%   x_all: the point x at each iteration

[n m] = size(A);
% Obtain k feasible points (along with the inital / final point)
x = cell(0);
%xfin = x_all(end) ;
xfin = xend;
x(1) = {0};
x(2) = {xfin};

k=50;
for i = 3: k+2
    xrad = rand(n, 1) - 0.5;
    while 1
        if max(A'*xrad) >= 1 || norm(xrad, inf) >= 1 % not feasible
            xrad = 0.1 * xrad;
        else
            break
        end
    end
    x(i) = {xrad};
end


% Compute M and m for all x
M = zeros(k+2, 1);
m= zeros(k+2, 1);
for i = 1: k+2
    [~, ~,h] = fun(x{i});
    m(i) = min(eig(h));
    M(i) = max(eig(h));
end
% Find the largest M and smallest m
m_min=min(m); M_max=max(M);


C = combnk(1:k+2,2);
size_C = size(C,1);
L = zeros(size_C,1);

for i =1:size_C
    p1 = x{C(i,1)};
    p2 = x{C(i,2)};
    [~,~,h1] = fun(p1);
    [~,~,h2] = fun(p2);
    L(i) = norm(h1-h2,2) / norm(p1-p2,2);
end

L_max = max(L);

end