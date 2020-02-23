function ex512


x = linspace(-0.5, 4);
plot(x, f(x), 'r-')
xlabel("$\lambda$", "Interpreter","latex")
ylabel("$g(\lambda)$", "Interpreter","latex")
xlim auto
%ylim auto
ylim([-10 10])

legend('Location','southeast')
legend("$g(\nu)$", "Interpreter","latex")
title("Lagrange Dual function $g(\nu)$","Interpreter","latex")
axis tight;
hold off

saveas(gcf,"simple_optimization_problem_dual","epsc")
function [y] = f(x)
    y = (-9 * (x.^2)) ./ (x + 1) + 1 + 8 .* x;
end


end