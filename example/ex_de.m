%% Diferencial Evolution - Example
% Antonio Horta Ribeiro - 2016
% Belo Horizonte - Brasil

clear all
clc
addpath('../')

% arg min f(x)
% subject to g(x)
f = @(x) peaks(x(1),x(2));
g = @(x) [(-x(1) + x(2) + 2.3); (x(1) + x(2) + 1.8)];

%% Differential Evolution(DE)
populSize = 14; % Population size.
nVars = 2; % Number of Variables.
upperBound = [3 3]; % Vector of upper bounds.
lowerBound = [-3 -3]; % Vector of lower bounds.

tic
[xBest, fBest, nGen, X, fx, gx] = de(f, nVars, lowerBound, upperBound, g);
toc

%% Plot 

x1 = -3:0.01:3;
x2 = -3:0.01:3;

[gridx1,gridx2]= meshgrid(x1,x2);

z = peaks(gridx1,gridx2);

figure(1)
hold on
[c, h] = contour(gridx1,gridx2,z, 10);
clabel(c,h);
colorbar
hold on
plot(x1,x1-2.3,'k--')
plot(x1,-x1 -1.8,'k--')

for k = 1:nGen
    % Plot distribution at each iter
    [c, h] = contour(gridx1,gridx2,z, 10);
    clabel(c,h);
    colorbar
    hold on
    for i = 1:populSize
    plot(X{k,i}(1),X{k,i}(2),'b*')
    hold on
    plot(x1,x1-2.3,'k--')
    plot(x1,-x1 -1.8,'k--')
    end
    hold off
    drawnow
end
