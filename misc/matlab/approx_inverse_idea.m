%% Code for testing the approximate inverse idea.
n = 100;     % Size of the matrix
d = 0.85;   % Interval for the approximation is [0 d]
deg = 15;   % Degree of the best approximating polynomial.

% Construct matrix with random eigenvalues between (0,1).
rng(3242018)
ev = sort(rand(n,1));
A = diag(ev);
Q = rand(n); [Q,~]=qr(Q);
A = Q'*A*Q;

% Matrix we want to approximate the inverse of:
S = eye(n) - A;

% Construct the best approximating polynomial
f = chebfun(@(x) 1./(1-x),[0 d]);
p = minimax(f, deg);

% Get the coefficients for this polynomial in the monomial basis
x = chebpts(deg+1,[0 d]);
cm = flipud(polyfit(x,p(x),deg)');

% Get the coefficients for this polynomial in the Chebyshev basis
cc = chebcoeffs(p);

% Random vector to multiply the polynomials by: p_n(A)*b
rng(1399184);
b = rand(n,1);

% Evaluate the product using monomial basis
ym = polyPrecondMono(b,A,cm);
% Evaluate the product using Chebyshev basis
yc = polyPrecondCheb(b,A,cc,[0 d]);

% Check how close the two methods are (mathematically the difference should
% be zero).
norm(ym-yc,inf)

% Check how close they approximate the inverse of S times b:
yexact = S\b;
norm(ym-yexact,inf)
norm(yc-yexact,inf)

%% Testing out the schemes as preconditioners
% No preconditioning
[x,flag,relrest,itert,resvect] = gmres(S,b,50,1e-12,100);
% Chebyshev version
[xc,c_flag,c_relrest,c_itert,c_resvect] = gmres(S,b,50,1e-12,100,@polyPrecondCheb,[],[],A,cc,[0 d]);
% Monomial version
[xm,m_flag,m_relrest,m_itert,m_resvect] = gmres(S,b,50,1e-12,100,@polyPrecondCheb,[],[],A,cc,[0 d]);
semilogy(1:numel(c_resvect),c_resvect,'o-',...
         1:numel(m_resvect),m_resvect,'x-',...
         1:numel(resvect),resvect,'s-')

%% Actual test problem
% load 4_level_mesh_130_patches_16x16cells_schur;
d = 0.95;
f = chebfun(@(x) 1./(1-x),[0 d]);
deg = 15;
p = minimax(f, deg);
cc = chebcoeffs(p);

xpts = chebpts(deg+1,[0 d]);
c = flipud(polyfit(xpts,p(xpts),deg)');

n = size(S,1);
A = speye(n)-S;
rng(3242018)
b = rand(n,1);
tic
[x,flag,relrest,itert,resvect] = gmres(S,b,50,1e-12,100,@polyPrecondCheb,[],[],A,cc,[0 d]);
toc
tic
[xm,flagm,relrestm,itertm,resvectm] = gmres(S,b,50,1e-12,100,@polyPrecondMono,[],[],A,c);
toc
semilogy(resvect,'rx-')
hold on
semilogy(resvectm,'bo-')

