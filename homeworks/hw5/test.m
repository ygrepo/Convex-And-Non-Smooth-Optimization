m = 50; n = 100;
A = randn(m,n);
x = ones(n,1);

if all(abs(A * x) < 1)
    disp("less than 1")
else
    disp("some element greater than 1")
    abs(A * x)
end