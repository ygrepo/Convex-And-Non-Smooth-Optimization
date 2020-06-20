function [f_all, g_all, x_all, z_all, u_all, r_all, s_all] = ...
admm_sdp(fun, ids, X_data, rho, x0, z0, u0, tol, maxit)

%ADMM for SDP
%   INPUT : 
%          fun: objective function
%          rho: quadratic penalization term
%          x0,z0,u0 : initial primal and dual iterates
%          tol      : absolute tolerance on residuals
%          maxit    : maximum number of iterations
%
%   OUTPUT : 
%           f_all, g_all: objective function in ADMM form
%           x_all, z_all: primal variables
%           u_all: dual variable
%           r_all: primal residuals (2-norm)
%           s_all: dual residuals (2-norm)

% To store values of objective (x,z,u) and residual r
f_all = [];
g_all = [];
x_all = {};
z_all = {};
u_all = {};
r_all = [];

% Initialize current xk and yk at the starting point x0
xk= x0;
zk = z0;
uk = u0;

% Dimension
[m, n]= size(X_data);
N=m+n;

% Useful identity matrices
IN = eye(N);
IN2 = eye(N^2);

% Nb of constraints
p = size(ids,1);

% Big H matrix for big KKT system during x ADMM update
H = zeros(N^2 + p);

% Upper left block
H(1:N^2,1:N^2) = rho*IN2;

for k=1:p
    % Extract the entry location (i,j) of the kˆth constraint
    i = ids(k,1);
    j = ids(k,2);

    % A_k is defined as zero everywhere except at specified entry (i,j)
    A_k = zeros(m,n);
    A_k(i,j) = 1;

    % Construct block matrix A tilde k encoding the kˆth equality constraint on x
    Atilde_k = 0.5*[zeros(m,m), A_k; A_k', zeros(n,n)];

    % Fill upper right and bottom left block of big KKT matrix H
    H(1:N^2,N^2+k) = vec(Atilde_k);
    H(N^2+k,1:N^2) = vec(Atilde_k)';
end    
    

% Relative tolerance
eps_rel= 1e-4;

% Main loop of the gradient
k = 1;
while (k <= maxit)

    % Record objective value
   [fk, gk] = fun(xk, zk);
   f_all(k) = fk;
   g_all(k) = gk;
   
   % Primal iterates
   x_all{end + 1} = xk;
   z_all{end + 1} = zk;
   
   
   % Dual iterates
   u_all{end + 1} = uk;
   
   % Primal residual
   rk = xk - zk;
   r_all(k) = norm(rk, 2);
   
   % RHS vector of KKT system
   h = zeros(N^2+p);
   % Dual part
   h(1:N^2) = vec(zk - uk -(1/rho)*IN);
   % Primal part (vector "b" with specified entries m ij's)
   for kp=1:p
       i = ids(kp,1);
       j = ids(kp,2);
       h(N^2 + p ) = X_data(i,j);
   end
   
   tmp = real(H\h);
   xkp1 = reshape(tmp(1:N^2), [N,N]);
   zkp1 = project_cone_psd(xkp1 + uk);
   ukp1 = uk + xkp1 - zkp1;
   
   % Dual residual
   sk = -rho * eye(N) * (zkp1 - zk);
   s_all(k) = norm(sk,2);
   
   % Stopping criterion
   eps_pri = sqrt(N)*tol + eps_rel*max(norm(xk,2), norm(zk,2));
   eps_dua = sqrt(N)*tol + eps_rel*norm(rho*uk, 2);
   if all((r_all(k) <= eps_pri)) && all((s_all(k) <= eps_dua))
       break
   end
   
   
   %Update iterates
   xk = xkp1;
   zk = zkp1;
   uk = ukp1;
   k = k +1;
    
end    

end

    