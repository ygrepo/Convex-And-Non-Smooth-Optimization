function [norm_tab, p_values] = plot_norm(A)
    p_values = 0.1:0.3:1;
    p_values =[p_values [2 3]];
    [~, n_p] = size(p_values);
    [~, n_A] = size(A);    
    norm_tab = zeros(n_A, n_p);
    for i=1:n_A
        X = A(:,i);
        for j=1:n_p
            p = p_values(j);
            norm_tab(i,j) = norm(X, p); 
        end
    end

    x_values = 1:1:n_p;
    line(x_values, norm_tab(1,:), 'Color', 'r', ...
        'LineStyle', ':', 'LineWidth', 2);
    hold on 
    line(x_values, norm_tab(2,:), 'Color', 'b', ...
        'LineStyle', '-.', 'LineWidth', 2);
    hold on 
    line(x_values, norm_tab(2,:), 'Color', 'g', ...
        'LineStyle', '--', 'LineWidth', 2);
    hold on 
    xticklabels({'0.1','0.4','0.7','1','2','3'})
    xlabel("p-norm")
    ylabel("norm-values")
    xlim([-inf inf]) 
    ylim([-inf inf])
    axis tight;
    legend('Location','northeast')
    legend('$\|V_1\|_p$','$\|V2\|_p$','$\|V3\|_p$','Interpreter','latex')
    title('p-norm of three random vectors')
    saveas(gcf,'p_norm_trend','epsc')