%% Best polynomial approximation preconditioner test.

% Load in the Schur complement matrix
load('../domain_decomp/matrix-files/4_level_mesh_130_patches_16x16cells_schur.mat');

% Schur complement has the form S = I-A.  We want to approximate inv(I-A)
% with a polynomial of A: p_n(A)
n = size(S,1);
A = speye(n)-S;

% Code for generating the best polynomial approximation to inv(I-A);
d = 0.95;
% % Uncomment to compute the coefficients on the fly
% deg = 15;
% f = chebfun(@(x) 1./(1-x),[0 d]);
% p = minimax(f, deg);
% % Coefficients in the Chebyshev basis
% cc = chebcoeffs(p);
% % Coefficients in the monomial basis
% xpts = chebpts(deg+1,[0 d]); cm = flipud(polyfit(xpts,p(xpts),deg)');

%
% Here are the coefficients from the above code
%
% Chebyshev
cc =[4.472135954953655e+00
     5.675247900481234e+00
     3.601012922685066e+00
     2.284885928634731e+00
     1.449787551186771e+00
     9.199076055378766e-01
     5.836924189936992e-01
     3.703598469934007e-01
     2.349977690621489e-01
     1.491089055767314e-01
     9.461139059090561e-02
     6.003206306517687e-02
     3.809106471898141e-02
     2.416923786484517e-02
     1.533567161022980e-02
     1.628851184599676e-02];
      
% Monomial
cm =[9.896646304586907e-01
     6.042710178142182e+00
    -4.072471800216493e+02
     1.304465012682393e+04
    -2.182169609100808e+05
     2.198110372567485e+06
    -1.444770634541706e+07
     6.509252624994410e+07
    -2.071251119656221e+08
     4.729408598029139e+08
    -7.779230221891359e+08
     9.137970461235683e+08
    -7.476486559170870e+08
     4.047130339719760e+08
    -1.302668757644882e+08
     1.887544576840549e+07];
 
% Create a vector to test solving x = S\b
rng(3242018)
b = rand(n,1);

% GMRES parameters
tol = 1e-12;
restart = 50;
maxiter = 100; 

% No preconditioning
tic
[x,flag,relrest,itert,resvect] = gmres(S,b,restart,tol,maxiter);
toc
% Preconditioning using Chebyshev basis
tic
[c_x,c_flag,c_relrest,c_itert,c_resvect] = gmres(S,b,50,1e-12,100,@polyPrecondCheb,[],[],A,cc,[0 d]);
toc
% Preconditioning using Monomial basis
tic
[m_x,m_flag,m_relrest,m_itert,m_resvect] = gmres(S,b,50,1e-12,100,@polyPrecondMono,[],[],A,cm);
toc
loglog(1:numel(c_resvect),c_resvect,'o-',...
         1:numel(m_resvect),m_resvect,'x-',...
         1:numel(resvect),resvect,'s-')
legend('Chebyshev','Monomial','No Preconditioning')

