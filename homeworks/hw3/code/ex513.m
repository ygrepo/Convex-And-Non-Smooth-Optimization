x = linspace(-1, 8);
plot(x, f(x), 'r-')
hold all
x = linspace(8, 10);
y = ones(size(x));
plot(x, y, 'r-')

xlabel("u")
ylabel("$p(u^*)$","Interpreter","latex")
xlim([-1 10]) 
ylim auto

legend('Location','southeast')
legend("$p(u^*)$", "Interpreter","latex")
title("Sensitivity analysis","Interpreter","latex")
axis tight;
hold off
saveas(gcf,"sensitivity_analysis",'epsc')

function [y] = f(x)
    y = 11 + x -6 * sqrt(1+x);
end