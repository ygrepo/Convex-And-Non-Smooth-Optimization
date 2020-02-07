% Compute the Chebyshev center of a polyhedron
function [x_sol, r_sol, x, y] = chebyshev_center_with_norm(A, b, p)
% The goal is to find the largest Euclidean ball (i.e. its center and
% radius) that lies in a polyhedron described by linear inequalites in this
% fashion: P = {x : a_i'*x <= b_i, i=1,...,m}

rng('default')
format long g
[~,n]=size(A);

% Build and execute model
fprintf(1,'Computing Chebyshev center...');
cvx_begin
    variable r(1)
    variable x_c(2)
    maximize ( r )
    for k=1:n
        A(:, k)' * x_c + r * norm(A(:, k), p) <= b(k);
    end
cvx_end
fprintf(1,'Done! \n');
x_sol = x_c;
r_sol = r;


% Display results
fprintf(1,"The Chebyshev center coordinates are: \n");
disp(x_c);
txt = "Radius of largest 'scaled unit ball' using norm p:" + num2str(p) + "\n";
fprintf(1,txt);
disp(r);

% Generate the figure
x = linspace(-2,2);
for k=1:n
    plot(x, -x * A(1,k)./A(2,k) + b(k)./A(2,k),"b-");
    hold on
end

n_vecs = 100;
[x, y] = gen_random_vectors(n_vecs, p);
x = x.* r + x_c(1);
y = y.* r + x_c(2);

for i=1:n_vecs
    plot(x(i), y(i), "r." );
    hold on
end
plot(x_c(1),x_c(2),'b*');
xlabel("x-axis")
ylabel("y-axis")

txt1 = "# inequalities:" + num2str(n);
tx2 = "norm-p:" + num2str(p);
title({"Largest 'scaled unit ball' lying in a 2D polyhedron", txt1, tx2});
text(x_c(1), x_c(2), "\leftarrow  center")
axis([-1 1 -1 1])
axis equal
hold off


function [x, y] = gen_random_vectors(n, p)
    r = randn(n, 2); % Use a large n
    for i=1:n
        norm_r = norm(r(i,:), p);
        r(i, :) = r(i, :) ./ norm_r;
    end
    x = r(:, 1);
    y = r(:, 2);
    

