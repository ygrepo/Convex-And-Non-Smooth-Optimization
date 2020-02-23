
x = 1:1:10;
x = linspace(0, 5);
CM = jet(10); 
for lambda = [0, 2, 4, 6]
%     disp(lambda)
%     disp(L(x, lambda))
    line(x, L(x, lambda), 'color', CM(lambda+1,:),'marker','.')
    hold on
end


function [y] = f(x)
    y = x.^2 + 1;
end

function [y] = L(x, lambda)
    y = f(x) + lambda .* (x-2) .* (x-4);
end