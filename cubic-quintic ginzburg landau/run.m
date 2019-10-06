clear all;
close all;
clc;

grid = cqglSolve;

figure(1);
colormap(cool);
s1 = surf(grid.t,grid.x,abs(grid.U) .^ 2,'FaceAlpha',1);
% s1.EdgeColor = 'none'; 
%title(sprintf('2D wave equation at t = %1.2f, con sigma = %1.2f y gamma = %1.2f',t(j),sigma, gamma),'Fontsize',11);
title("dark-dark solitary waves for eq(1) (Perturbed) A(x, t)");
xlabel('time'); ylabel('x'); zlabel("|A|^2");