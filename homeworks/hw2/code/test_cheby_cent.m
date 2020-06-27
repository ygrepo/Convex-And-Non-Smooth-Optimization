%%
clear all; close all;
% cheby_test.m for part (a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CVX Example Code %%%%%%%%%%%%%%%%%%%%%%% %%%%%
% Boyd & Vandenberghe, "Convex Optimization"
% Joëlle Skaf - 08/16/05
% (a figure is generated)
%
% The goal is to find the largest Euclidean ball (i.e. its center and radius) 
% that lies in a polyhedron described by linear inequalites in this
% fashion: P = {x : a_i'*x <= b_i, i=1,...,m} where x is in R^2
% Generate the input data
a1 = [ 2; 1]; a2 = [ 2; -1]; a3 = [-1; 2]; a4 = [-1; -2]; b = ones(4,1);
% Create and solve the model
cvx_begin 
    variable r(1)
    variable x_c(2)
    maximize ( r )
    a1'*x_c + r*norm(a1,2) <= b(1); 
    a2'*x_c + r*norm(a2,2) <= b(2); 
    a3'*x_c + r*norm(a3,2) <= b(3); 
    a4'*x_c + r*norm(a4,2) <= b(4);
cvx_end
% Generate the figure
x = linspace(-2,2);
theta = 0:pi/100:2*pi;
plot( x, -x*a1(1)./a1(2) + b(1)./a1(2),'b-');
hold on
plot( x, -x*a2(1)./a2(2) + b(2)./a2(2),'b-');
plot( x, -x*a3(1)./a3(2) + b(3)./a3(2),'b-');
plot( x, -x*a4(1)./a4(2) + b(4)./a4(2),'b-');
plot( x_c(1) + r*cos(theta), x_c(2) + r*sin(theta), 'r'); 
plot(x_c(1),x_c(2),'k+')
xlabel('x_1')
ylabel('x_2')
title('Largest Euclidean ball lying in a 2D polyhedron');
axis([-1 1 -1 1])
axis equal 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % verification that cheby_center(A,b) does the same thing as CVX example
A = [a1 a2 a3 a4]; %put the same inequalities into matrix A figure(2);
[my_x_c, my_r] = cheby_cent(A,b,2);
disp(['difference between example and my x_c in 2-norm is:',num2str(norm(my_x_c-x_c))]);
disp(['difference between example and my r in magnitude is:',num2str(abs(my_r-r))]);


%%
% test on a different set of inequalities where the polyhedron is non 
% empty, namely the strip defined by x+y >= -1 and x+y <= 1, contains at
% least the l1 ball, known that the max radius is sqrt(2)/2;
A = [-1 1; -1 1]; b = [1, 1];
figure(3);
[new_x_c,new_r] = cheby_cent(A,b); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part (c)
% test for Chebyshev center with generic p-norm, 0 < p <= Inf a1 = [ 2; 1];
a2 = [ 2; -1];
a3 = [-1; 2];
a4 = [-1; -2];
b = ones(4,1);
A = [a1 a2 a3 a4]; %assemble matrix figure(1); 

%1-norm
[x_1, r_1] = cheby_cent(A,b,1);
figure(2);

%1.5-norm
[x_1p5, r_1p5] = cheby_cent(A,b,1.5);
figure(3); 

%Inf-norm
[x_inf,r_inf] = cheby_cent(A,b,Inf);

% part (d)
figure(4); 
%p-norm, p<1 
[x_p,r_p] = cheby_cent(A,b,7/10);

