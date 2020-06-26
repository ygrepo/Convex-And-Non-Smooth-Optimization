function [t]= vec_solver(x,b)
[n,~] = size(x)
t = inf;
if sum(b >=0)
   t = max(-b./x);
   t = 1/t;
end    
end