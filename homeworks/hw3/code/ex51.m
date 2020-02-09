
function ex51
x = linspace(0,10);
plot(x, f(x), 'r-')
hold on
plot(x,x, 'b-')
xlabel("x")
ylabel("y")

hold off

function [y] = f(x)
    y = x.^2 + 1;
    