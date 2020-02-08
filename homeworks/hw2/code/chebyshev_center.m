% Compute the Chebyshev center of a polyhedron
% Boyd & Vandenberghe "Convex Optimization"
function [x_sol, r_sol] = chebyshev_center(A, b)
% The goal is to find the largest Euclidean ball (i.e. its center and
% radius) that lies in a polyhedron described by linear inequalites in this
% fashion: P = {x : a_i'*x <= b_i, i=1,...,m}

% Generate the data
% randn('state',0);
% m = 2; n = 10*m;
% A = randn(m,n);
% b = ones(n,1);
% b = ones(4,1);
% A = [2 2 -1 -1; 1 -1 2 -2];

[~,n]=size(A);

% Build and execute model
fprintf(1,'Computing Chebyshev center...');
cvx_begin
    variable r(1)
    variable x_c(2)
    maximize ( r )
    for k=1:n
        A(:, k)'*x_c + r*norm(A(:, k),2) <= b(k);
    end
cvx_end
fprintf(1,'Done! \n');
x_sol = x_c;
r_sol = r;

% Display results
fprintf(1,'The Chebyshev center coordinates are: \n');
disp(x_c);
fprintf(1,'The radius of the largest Euclidean ball is: \n');
disp(r);

% Generate the figure
x = linspace(-2,2);
for k=1:n
    plot(x, -x * A(1,k)./A(2,k) + b(k)./A(2,k),"b-");
    hold on
end
theta = 0:pi/100:2*pi;
plot( x_c(1) + r*cos(theta), x_c(2) + r*sin(theta), "r" );
plot(x_c(1),x_c(2),'b*');
xlabel("x_1")
ylabel("x_2")

txt = "# inequalities:" + num2str(n);
title({"Largest Euclidean ball lying in a 2D polyhedron", txt});
text(x_c(1), x_c(2), "\leftarrow  center")
axis([-1 1 -1 1])
axis equal
hold off
txt = "chebyshev_center_" + num2str(n);
saveas(gcf,txt,'epsc')
%saveas(gcf,txt,'png')


