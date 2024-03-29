
function ex51


% Add lines
line([2 2],[0 f(6)]);
line([4 4],[0 f(6)]);

% Add a patch
gray = [0.9 0.9 0.9];
patch([2 4 4 2],[0 0 f(6) f(6)],gray,'HandleVisibility','off')
alpha(0.9);
hold on

x = linspace(0, 5);
plot(x, f(x), 'r-')
text(2, f(2), " optimal point and value", 'FontSize',12)
hold on

x = linspace(0, 4);
plot(x, L(x, 1), 'b-')
plot(x, L(x, 2), 'g-')
plot(x, L(x, 3), 'c-')

xlabel("x")
ylabel("y")
xlim([0 5]) 
ylim([0 25])

legend('Location','northeast')
legend("$f_0(x)$", "$L(x,1)$","$L(x,2)$","$L(x,3)$", "Interpreter","latex")
title("$f_0(x)$ and $L(x,\lambda)$","Interpreter","latex")
axis tight;
hold off

saveas(gcf,"simple_optimization_problem","epsc")
    
function [y] = f(x)
    y = x.^2 + 1;
end

function [y] = L(x, lambda)
    y = f(x) + lambda .* (x-2) .* (x-4);
end
      
end
    