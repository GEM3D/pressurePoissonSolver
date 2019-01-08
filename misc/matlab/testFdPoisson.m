%% Dirichlet test problem
% Script for testing fd2poisson over the square [a,b]x[a,b]
a = 0; 
b = 1;
k = 7;
m = 2^k;  % Number of interior grid points in one direction

f = @(x,y) -5*pi^2*sin(pi*x).*cos(2*pi*y);  % Laplacian(u) = f
g = @(x,y) sin(pi*x).*cos(2*pi*y);          % u = g on Boundary
uexact = @(x,y) g(x,y);                     % Exact solution is g.

% Compute and time the solution
tic
[usp,x,y] = fd2poissonsp(f,g,a,b,m,m);
gesp = toc;

tic
[ufft,x,y] = fd2poissonfft(f,g,[],[],a,b,m,m,1);
gefft = toc


%%

fprintf('Direct Gaussian elimination takes %d s\n',gesp);
fprintf('FFT solver takes %d s\n',gefft);
fprintf('Norm of the difference between the solutions is %1.4e\n',norm(usp(:)-ufft(:)));

figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,ufft), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
title(strcat('Numerical Solution to Poisson Equation, m=',num2str(m)));

% Plot error
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,ufft-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
title(strcat('Error, m=',num2str(m)));

%% Test of all 10 types of (mixed) boundary conditions
a = 0; 
b = 1;
k = 10;
m = 2^k;  % Number of interior grid points in one direction

f = @(x,y) -5*pi^2*sin(pi*(x+1/sqrt(99))).*cos(2*pi*(y-1/sqrt(101)));      % Laplacian(u) = f
dfun = @(x,y) sin(pi*(x+1/sqrt(99))).*cos(2*pi*(y-1/sqrt(101)));           % u = g on Boundary
uexact = @(x,y) dfun(x,y);                      % Exact solution is dfun.
nfunx = @(x,y) pi*cos(pi*(x+1/sqrt(99))).*cos(2*pi*(y-1/sqrt(101)));       % Neumann in x: du/dx 
nfuny = @(x,y) -2*pi*sin(pi*(x+1/sqrt(99))).*sin(2*pi*(y-1/sqrt(101)));    % Neumann in y: du/dy

% Dirichlet NSEW
[u,x,y] = fd2poissonfft(f,dfun,[],[],a,b,m,m,1);
fprintf('Error in Dirichlet NSEW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

%% Neumann NSEW
m = 2^11;
tic
[u,x,y] = fd2poissonfft(f,[],nfunx,nfuny,a,b,m,m,10);
fprintf('Error in Neumann NSEW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));
toc
%%

% Neumann SW Dirichlet NE
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,2);
fprintf('Error in Neumann SW Dirichlet NE solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann W Dirichlet NSE
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,3);
fprintf('Error in Neumann W Dirichlet NSE solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann NW Dirichlet SE
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,4);
fprintf('Error in Neumann NW Dirichlet SE solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann N Dirichlet SEW
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,5);
fprintf('Error in Neumann N Dirichlet SEW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann NE Dirichlet SW
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,6);
fprintf('Error in Neumann NE Dirichlet SW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann E Dirichlet NSW
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,7);
fprintf('Error in Neumann E Dirichlet NSW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann SE Dirichlet NW
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,8);
fprintf('Error in Neumann SE Dirichlet NW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

% Neumann S Dirichlet NEW
[u,x,y] = fd2poissonfft(f,dfun,nfunx,nfuny,a,b,m,m,9);
fprintf('Error in Neumann S Dirichlet NEW solution: %1.4e\n',norm(u(:)-uexact(x(:),y(:)),inf));

