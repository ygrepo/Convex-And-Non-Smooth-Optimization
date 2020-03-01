W = load("hw4data1.mat").W;
[n, ~] = size(W)

% SDP relaxation
fprintf(1,'Solving the SDP relaxation of the two-way partitioning problem...');

cvx_begin sdp
    variable X(n,n) symmetric
    minimize ( trace(W*X) )
    diag(X) == 1;
    X >= 0;
cvx_end

fprintf(1,'Done! \n');
opt2 = cvx_optval;
disp(opt2)