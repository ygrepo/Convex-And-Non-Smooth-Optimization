function [Xs, ids, fnorms] = low_rank_completion(XDATA, tol)

% This function runs nuclear norm minimization SDPs using CVX adding 
% information from XDATA. We stop the iterations when X has become close 
% enough to XDATA (within 'tol' in Froebennius norm).
%   INPUT
%       XDATA: data set
%       tol : tolerance (stopping condition)
%   OUTPUT
%       Xs: sequence of solutions to SDP
%       ids: pairs of indices (i,j)
%       fnorms:froebenius norm errors of X's


% Get size
[m,n] = size(XDATA);
N = m * n;

% to store optimal solutions X's
Xs = {};

% to store Froebenius norm erros of X's
fnorms = [];

% Generate random indices (i,j) from [1,m]*[1,m] without repetition
ids = [];
while size(ids,1) ~= N
    ids = unique([ids; ceil(m*rand), ceil(n*rand)], "row");
end
% Shuffle indices
ids = ids(randperm(size(ids,1)),:);

% Nb of constraints = nb of iterations
n_iters = 1;

% While there are still constraints to add
while n_iters <= N
    
    % solve SDP
    cvx_begin sdp
        variable W1(m,m) symmetric
        variable W2(n,n) symmetric
        variable X(m,n)
        minimize(0.5 * (trace(W1) + trace(W2)) );
        [W1 , X; X',W2]>=0;
        % add constraints
        for k=1:n_iters
            r = ids(k,1);
            c = ids(k,2);
            X(r,c) == XDATA(r,c);
        end
    cvx_end
    
    % Record X
    Xs{end + 1} = full(X);
    % Record Froebenius norm error
    fnorms = [fnorms; norm(full(X) - XDATA, "fro")];
    
    % If X close enough to XDATA , we stop using Froebenius norm
    if norm(full(X) - XDATA, "fro") <= tol
        break
    end
    
    % update counter
    n_iters = n_iters + 1;
end


end