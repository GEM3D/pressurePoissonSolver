% Script for testing fd2poisson over the square [a,b]x[a,b]
a = 0; 
b = 1;
k = 8;
m = 2^k;  % Number of interior grid points in one direction

f = @(x,y) -5*pi^2*sin(pi*x).*cos(2*pi*y);  % Laplacian(u) = f
g = @(x,y) sin(pi*x).*cos(2*pi*y);          % u = g on Boundary
uexact = @(x,y) g(x,y);                     % Exact solution is g.

% Compute and time the solution
tic
[usp,x,y] = fd2poissonsp(f,g,a,b,m,m/2);
gesp = toc;

tic
[ufft,x,y] = fd2poissonfft(f,g,a,b,m,m/2);
gefft = toc;

fprintf('Direct Gaussian elimination takes %d s\n',gesp);
fprintf('FFT solver takes %d s\n',gefft);
fprintf('Norm of the difference between the solutions is %1.4e\n',norm(usp(:)-ufft(:)));

%% Plot solution
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,ufft), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
title(strcat('Numerical Solution to Poisson Equation, m=',num2str(m)));

% Plot error
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,ufft-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
title(strcat('Error, m=',num2str(m)));
