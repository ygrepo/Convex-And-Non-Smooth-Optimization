function [f,g] = sdp_admm_fun(x, z)

% f(x) objective : trace of x
f = trace(x);

% g(z) objective : indicator function of Sn+ for z
g = all(eig(z)>= 0) && issymmetric(z);
end