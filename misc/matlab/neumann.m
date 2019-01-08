function neumann(method,sc)
% NEUMANN solves a 1d Neumann problem
%
% NEUMANN(METHOD,SOLVER_CHOICE) solves using a two step procesure (METHOD=1) or by
% setting an interior value (METHOD=2).   SOLVER_CHOICE is one of 'direct',
% 'pcg','bicgstab' or 'gmres'.
%


global A h solver_choice;

if (nargin < 2)
    sc = 'direct';
    if (nargin < 1)
        method = 1;    % 1 : Two step procedure;  2 : set interior value.
    end
end

plot_eig = false;    % Plot eigenvalues and return

% In the two step version, we can impose either mean on the boundary or the
% interior.  Choices : 
% 
%  enforce_mean_bc       : Enforce zero mean on the boundary
%  enforce_mean_domain   : Enforce zero mean in the domain.
% 

twostep_choice = 'enforce_mean_domain';   


solver_choice = sc;
tol = 1e-8;

% --------------------------------
% Physical parameters
% --------------------------------
u_true = @(x) cos(2*pi*x);
ux_true = @(x) -(2*pi)*sin(2*pi*x);
uxx_true = @(x) -(2*pi)^2*cos(2*pi*x);

% u_true = @(x) x - 0.5;
% ux_true = @(x) ones(size(x));
% uxx_true = @(x) zeros(size(x));

% Domain
ax = 0;
bx = 1;

% --------------------------------
% Numerical parameters
% --------------------------------

N = input('Input N : ');
h = (bx-ax)/N;

xe = linspace(0,1,N+1)';
xc = xe(1:end-1) + h/2;

uc = u_true(xc);
ul = u_true(xe(1));
ur = u_true(xe(end));
uxl = ux_true(xe(1));
uxr = ux_true(xe(end));

% Check compatibility
s = sum(uxx_true(xc));   % FTC!
tol_input = 1e-8;
if abs(s - (uxr-uxl)) > tol_input
    warning('Compatibility condition not satisfied');
end
if abs(ul + ur) > tol_input
    warning('The true solution on the boundary does not have zero mean')
end
if (abs(sum(uc)*h) > tol_input)
    warning('The true solution does not have zero mean');
end

switch method
    case 1
        % Two step procedure;  Enforce zero mean either on the boundary or 
        % in the domain.
        
        switch twostep_choice
            case 'enforce_mean_bc'
                enforce = @F;
            case 'enforce_mean_domain'
                enforce = @G;
        end
        
        % Set up matrix with Dirichlet condition at left edge; Neumann at
        % right.
        z = ones(N,1);
        A = spdiags([z -2*z z],-1:1,N,N);
        A(1,1) = -3;
        A(N,N) = -1;
        
        if (plot_eig)
            figure(2)
            clf;
            plot(eig(full(A)),'.','markersize',30);
            set(gca,'fontsize',16);
            hold on;
            plot(xlim,[0 0],'k--');
            plot([0 0],ylim,'k--');
            hold off;
            fprintf('Minimum eigenvalue : %12.4e\n',min(abs(eig(full(A)))));
            return
        end
        
        % Right hand side
        f = h^2*uxx_true(xc);
        
        % Solve for unknown value at left end point by solving F(g) = 0
        % where
        %          F(g) = a*g + b
        %
        
        tic;
        b = enforce(f,0,uxl,uxr);      % set g=0
        a = enforce(f,1,uxl,uxr) - b;  % set g=1
        g = -b/a;
        
        % Final solution
        [c,u,flag,relres,iter,resvec] = enforce(f,g,uxl,uxr);   % c should be 0 here
        t1 = toc;

    case 2
        % Neumann conditions at both ends with identity stencil at one
        % point
        
        z = ones(N,1);
        A = spdiags([z -2*z z],-1:1,N,N);
        A(1,1) = -1;
        A(N,N) = -1;
        
        if (plot_eig)
            figure(2)
            clf
            plot(sort(eig(full(A))),'b.','markersize',30);
            hold on;
            fprintf('Minimum eigenvalue : %12.4e\n',min(abs(eig(full(A)))));
        end
        
        % Set one point to identity stencil.  By setting to -2, 
        % we get a positive definite matrix (although PCG still doesn't
        % really work (???)
        A(N/2,(N/2-1):(N/2+1)) = [0 -2 0];
        
        if (plot_eig)
            plot(sort(eig(full(A))),'ro','markersize',7);
            plot(xlim,[0 0],'k--');
            plot([0 0],ylim,'k--');
            set(gca,'fontsize',16);
            hold off;
            legend({'Neumann', 'Modified Neumann'},'fontsize',16,'location','northwest')
            fprintf('Minimum eigenvalue : %12.4e\n',min(abs(eig(full(A)))));
            return
        end
        
        % Right hand side
        f = h^2*uxx_true(xc);
        
        % Set one point to identity stencil.
        f(N/2) = -2*u_true(xc(N/2));
        
        % BC
        bc = zeros(N,1);
        bc(1) = h*uxl;       % Dirichlet condition at left end
        bc(end) = -h*uxr;
        
        % Solver linear sysytem
        f = f + bc;
        maxit = 2*N;
        
        tic
        switch solver_choice
            case 'direct'
                u = A\f;
                flag = 0;
                relres = norm(A*u-f)/norm(f);
                iter = 1;
                resvec = relres;
                
            case 'bicgstab'
                [u,flag,relres,iter,resvec] = bicgstab(-A,-f,tol,maxit);
            case 'pcg'
                [u,flag,relres,iter,resvec] = pcg(-A,-f,tol,maxit);
            case 'gmres'
                restart = 20;
                maxit = 5*N;
                [u,flag,relres,iter,resvec] = gmres(-A,-f,restart,tol,maxit);                
        end
        t1 = toc;
    case 3
        % Mean value of the true solution; used for setting the constant in
        % the computed solution.
        meanSol = h*sum(u_true(xc));
        
        % Right hand side
        f = h^2*uxx_true(xc);
        tic
        fhat = dct2(f);
        den = -4*sin((0:N-1)'*pi/(2*N)).^2;
        phat = fhat./den;
        phat(1) = meanSol/h;  % Pick the constant to be the mean of the solution
        u = idct2(phat);
        t1 = toc;
        flag = 0;

        % Compute the residual
        z = ones(N,1);
        A = spdiags([z -2*z z],-1:1,N,N);
        A(1,1) = -1;
        A(N,N) = -1;
        relres = norm(A*u-f)/norm(f);
        iter = 1;
        resvec = relres;
end

% Plot the solution
figure(1);
clf
plot(xc,u,'k.','markersize',30);
hold on;

xf = linspace(0,1,500);
uf = u_true(xf);
plot(xf,uf,'r');
set(gca,'fontsize',16);


err = norm(u - u_true(xc),1)*h;

if (flag > 0)
    warning(sprintf('%s : Convergence not successful;  flag = %d\n',solver_choice,flag));
else
    fprintf(sprintf('%s : Successful convergence\n',solver_choice));
end
fprintf('\n');
fprintf('%15s %12d\n','N',N);
fprintf('%15s %12d\n','Flag',flag);
fprintf('%15s %12.4e\n','Rel. Residual',relres);
fprintf('%15s %12d\n','Iterations',round(iter(1)));  % BICGSTAB results non-integer value
fprintf('%15s %12.2e\n','Timing',t1);
fprintf('%15s %12.4e\n','Error',err);
fprintf('\n');

figure(3);
clf
semilogy(resvec,'.-','markersize',20);
set(gca,'fontsize',16);
title('Residual norm','fontsize',16);
xlabel('Iteration');
ylabel('Residual');

shg

end


function [c,u,flag,relres,iter,resvec] = F(f,g,uxl,uxr)

global A h solver_choice tol

N = length(A);

bc = zeros(N,1);
bc(1) = -2*g;       % Dirichlet condition at left end
bc(end) = -h*uxr;   % Neumann condition at the right end

% Solve problem with u=g at left end point; u_x = uxl at right end point
b = f + bc;

switch solver_choice
    case 'direct'
        u = A\b;
        flag = 0;
        relres = norm(A*u-b)/norm(b);
        iter = 1;
        resvec = relres;
        
    case 'bicgstab'
        [u,flag,relres,iter,resvec] = bicgstab(-A,-b,tol,2*N);
    case 'pcg'
        [u,flag,relres,iter,resvec] = pcg(-A,-b,tol,2*N);
    case 'gmres'
        restart = 20;
        maxit = 5*N;
        [u,flag,relres,iter,resvec] = gmres(-A,-b,restart,tol,maxit);
end

% Use Neumann condition at right end to get value at right end point
ur_ghost = u(end) + h*uxr;
ur = (u(end) + ur_ghost)/2;

% Impose Neumann condition at left edge as constraint :
%
%          (u(1) - ul_ghost)/h = uxl
%
% Use ur + ul = 0 to get g = ul = -ur and substitute into
% expression for Dirichlet condition
%
%             (u(1) + u_ghost)/2 = g
% to get

ul_ghost = 2*(-ur) - u(1);   % ul_ghost = 2*g - u(1)

% Final constraint condition.
c = (u(1) - ul_ghost)/h - uxl;


end

% Enforce zero mean in the domain
function [c,u,flag,relres,iter,resvec] = G(f,g,uxl,uxr)

global A h solver_choice tol

N = length(A);

bc = zeros(N,1);
bc(1) = -2*g;       % Dirichlet condition at left end
bc(end) = -h*uxr;   % Neumann condition at the right end

% Solve problem with u=g at left end point; u_x = uxl at right end point
b = f + bc;

switch solver_choice
    case 'direct'
        u = A\b;
        flag = 0;
        relres = norm(A*u-b)/norm(b);
        iter = 1;
        resvec = relres;
        
    case 'bicgstab'
        [u,flag,relres,iter,resvec] = bicgstab(-A,-b,tol,2*N);
    case 'pcg'
        [u,flag,relres,iter,resvec] = pcg(-A,-b,tol,2*N);
    case 'gmres'
        restart = 20;
        maxit = 5*N;
        [u,flag,relres,iter,resvec] = gmres(-A,-b,restart,tol,maxit);
end

% Use Neumann condition at right end to get value at right end point
ur_ghost = u(end) + h*uxr;
ur = (u(end) + ur_ghost)/2;

% Impose Neumann condition at left edge as constraint :
%
%          (u(1) - ul_ghost)/h = uxl
%
% Use sum(u) = 0 to get sum(u) = u(1) + (sum(u) - u(1)).  We want sum(u) = 0, 
% so we get u1 = -sum(u) + u(1).  Use this to enforce Neumann constraint
%
%             (u1 + ul_ghost)/2 = g
% to get

ul_ghost = 2*(-sum(u) + u(1)) - u(1);   % ul_ghost = 2*g - u(1)

% Final constraint condition.
c = (u(1) - ul_ghost)/h - uxl;


end


function coeffs = dct2(x,dim)
%  DCT2   Discrete Cosine Transform Type II computed using the fast Fourier Transform.
%     X = dct2(x) computes the Discrete Cosine Transform Type II (DCT-II) of the columns of X.
%
%     X = dct2(x,dim) computes the DCT along the dimension specified.
%     if dim = 1 (default) then the DCT is along the columns.
%     if dim = 2 then the DCT is along the rows.
%
%  See also idct2, dct, idct, dst, dst2, idst, idst2.

[m,n] = size(x);

if nargin == 1
    dim = 1;
end

% Trivial case (constant):
if ( m <= 1 )
    coeffs = x;
    return
end

if dim == 2
    j = (0:n-1);
    % Mirror the values for FFT:
    tmp = [x x(:,n:-1:1)];
    scale = n;
else
    j = (0:m-1)';
    % Mirror the values for FFT:
    tmp = [x; x(m:-1:1,:)];
    scale = m;
end

% Pre-compute the weight vector:
w = 2*exp(1i*j*pi/(2*scale));

% tmp = [x; x(m:-1:1,:)];
coeffs = ifft(tmp,[],dim);

% Truncate, flip the order, and multiply the weight vector:
coeffs = bsxfun(@times, w, coeffs(1:m,1:n));

% Scale the coefficient for the constant term:
coeffs = (scale)/2*coeffs;

% Post-process:
if ( isreal(x) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*x) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end

function coeffs = idct2(x,dim)
%  IDCT2   Inverse Discrete Cosine Transform Type II computed using the fast Fourier Transform.
%     X = idct2(x) computes the Inverse Discrete Cosine Transform Type II (DCTII) of the columns of X.
%
%     X = idct2(x,dim) computes the inverse DCT along the dimension specified.
%     if dim = 1 (default) then the inverse DCT is along the columns.
%     if dim = 2 then the inverse DCT is along the rows.
%
%  See also dct2.

[m,n] = size(x);

if nargin == 1
    dim = 1;
end

% Trivial case (constant):
if ( m <= 1 )
    coeffs = x;
    return
end

if dim == 2
    % Account for constant term
    x(:,1) = .5*x(:,1); 

    % Mirror the values for FFT:
    tmp = [x ones(m,1) x(:,end:-1:2)];
    scale = n;

    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*n-1)*pi/(2*n))/2);
    w([1, n+1]) = [2*w(1) 0];
    w(n+2:end) = -w(n+2:end);

else
    % Account for constant term
    x(1,:) = .5*x(1,:); 

    % Mirror the values for FFT:
    tmp = [x ; ones(1, n) ; x(end:-1:2,:)];
    scale = m;
    
    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*m-1)*pi/(2*m))/2).';
    w([1, m+1]) = [2*w(1); 0];
    w(m+2:end) = -w(m+2:end);
end

% Apply the weight vector:
tmp = bsxfun(@times, tmp, w);

coeffs = fft(tmp,[],dim);

% Truncate, flip the order, and scale:
% coeffs = (2/scale)*coeffs(j,k);
% Does fftshift work here? 
coeffs = (2/scale)*coeffs(1:m,1:n);

% Post-process:
if ( isreal(x) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*x) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end
