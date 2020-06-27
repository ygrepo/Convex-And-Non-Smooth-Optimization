function [x_c,r] = cheby_cent(A,b,p)
% Modification on the CVX code example of finding the largest p-norm
% ball contained in in a polyhedron described by A'x <= b,
% Allows an arbitrary number of inequality constraints
%
% INPUT: %
%
% OUTPUT: x_c -- vector of size 2 x 1, center p-norm ball
% r -- radius of largest p-norm ball


if nargin < 3
    p = 2; %default is 2-norm
end

% Create and solve the model
cvx_begin
variable r(1)
variable x_c(2)
maximize ( r )
if p == Inf %p=Inf, conjugate is 1, cannot use p / (p-1) formula
    for i=1:size(A,2)
        a = A(:,i);
        a'*x_c + r*norm(a,1) <= b(i);
    end
else
    q = max(p,1); %if p < 1, solve with p=1
    for i=1:size(A,2)
        a = A(:,i);
        a'*x_c + r*norm(a,q/(q-1)) <= b(i);
    end
end
cvx_end

% Generate the figure
x = linspace(-2,2);
for i=1:size(A,2) %plot all the lines
    a = A(:,i);
    plot(x, -x*a(1)./a(2) + b(i)./a(2),'b-','LineWidth',1.5);
    hold on;
end
if p == Inf
    %plot the Inf-norm ball (square), all 4 sides
    xline = linspace(-r,r);
    yline = xline;
    plot(x_c(1)+xline, x_c(2)+r*ones(size(xline)),'r','LineWidth',1.5);
    hold on;
    plot(x_c(1)+xline, x_c(2)-r*ones(size(xline)),'r','LineWidth',1.5);
    hold on;
    plot(x_c(1)+r*ones(size(xline)),x_c(2) + yline,'r','LineWidth',1.5);
    hold on;
    plot(x_c(1)-r*ones(size(xline)),x_c(2) + yline,'r','LineWidth',1.5);
else
    %represent x^p + y^p = r^p as a graph for 0<x<r, rotate 4 times
    xp = linspace(0,r);
    rot = [0 -1; 1 0]; %rotation by pi/2
    y = abs((r^p - xp.^p).^(1/p));
    for i=0:3
        temp = (rot^i)*[xp;y];
        plot(x_c(1) + temp(1,:), x_c(2) + temp(2,:),'r','LineWidth',1.5);
        hold on;
    end
end
plot(x_c(1),x_c(2),'x','MarkerSize',8,'LineWidth',1.5)
xlabel('x_1')
ylabel('x_2')
title(['Largest ',num2str(p),'-norm ball lying in a 2D polyhedron']);
axis([-1 1 -1 1])
axis equal
set(gca,'FontSize',12);
end
