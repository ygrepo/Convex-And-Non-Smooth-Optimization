%% hw9 - 1.

clear all;
clc;

M = 3;
N = 2;
X = randn(M,N);
svd(X)

% compute primal SDP
cvx_begin sdp
    variable W(M+N,M+N) 
    minimize(0.5 * trace(eye(M+N)' * W));
    W(1:M,M+1:end)==X;
    W>=0;
cvx_end
pstar = cvx_optval;

%compute dual SDP
cvx_begin sdp
    variable Y(M,N)
    maximize(trace(X' * Y));
   norm(Y,2) <= 1;
cvx_end
dstar = cvx_optval;

% Check primal solution
X
W

% check dual solution
trace(X' * Y)
[U,S,V]= svd(X,0);
trace(S)

% check strong duality
pstar
dstar
norm(svd(X),1)

% check complementarity
S = 0.5 * [eye(M), -Y; -Y' eye(N)];
W * S
