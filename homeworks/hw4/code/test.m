W = load("hw4data1.mat").W;
[n, ~] = size(W)
% Lagrange dual
fprintf(1,'Solving the dual of the two-way partitioning problem...');

cvx_begin sdp
    variable nu(n)
    maximize -sum(nu, 1)
    %W + diag(nu) >= 0;
    W + diag(nu) == semidefinite(n);
cvx_end

fprintf(1,'Done! \n');
opt1 = cvx_optval;