randn('state',0);
n = 10;
W = randn(n); W = 0.5*(W + W');

% Lagrange dual
fprintf(1,'Solving the dual of the two-way partitioning problem...');

cvx_begin sdp
    variable W(2,2) symmetric
    W(1,2) == 1
    W(2,1) == 1
    minimize ( W(1,1) )
    W > 0;
cvx_end

fprintf(1,'Done! \n');
opt1 = cvx_optval;

% SDP relaxation
fprintf(1,'Solving the SDP relaxation of the two-way partitioning problem...');

