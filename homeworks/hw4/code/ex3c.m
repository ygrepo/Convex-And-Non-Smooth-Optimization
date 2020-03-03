% Implementation of MAXCUT algorithm  
% Using data set 1 and data set 2 to determine
% maxcut and the corresponding partition.

clear
cvx_quiet true

W = load("hw4data1.mat").W;
[n, ~] = size(W);
[X, cut_value] = maxcut(n, W);
V = chol(X);
fprintf("Max. cut value data set 1:%5.2f\n", cut_value)
[partition] = get_partition(V, n);
fprintf("Partition for data set 1\n")
disp(partition')

clear
W = load("hw4data2.mat").W;
[n, ~] = size(W);
[X, cut_value] = maxcut(n, W);
V = chol(X);
fprintf("Max. cut value data set 2:%5.2f\n", cut_value)
[partition] = get_partition(V, n);
fprintf("Partition for data set 2\n")
disp(partition')

function [X, cut_value] = maxcut(n, W)
cvx_begin sdp
    variable X(n,n) symmetric
    maximize 0.25 .* sum(sum(W .* (ones(n,n) - X))) 
    diag(X) == 1;
    X >= 0;
cvx_end
cut_value = cvx_optval;

end

function [partition] = get_partition(V, n)
partition = zeros(n, 1);
r = (1/sqrt(n)) .* ones(n, 1);

for i=1:n
    if dot(V(i,:), r) >= 0
        partition(i) = i;
    else
        partition(i) = 0;
    end
end    
end