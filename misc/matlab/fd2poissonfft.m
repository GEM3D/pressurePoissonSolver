function [u,x,y] = fd2poissonfft(ffun,dfun,nfunx,nfuny,a,b,mx,my,bctype)
% Numerical approximation to Poisson's equation over the square [a,b]x[a,b]
% with boundary conditions specified by the 10 types indicated below:
%
%          N
%      --------
%     |        |
%   W |        | E
%     |        |
%      --------
%          S
%
% 1.  Dirichlet (NSEW)
% 2.  Neumann (SW) Dirichlet (NE)
% 3.  Neumann (W)  Dirichlet (NSE)
% 4.  Neumann (NW) Dirichlet (SE)
% 5.  Neumann (N)  Dirichlet (SEW)
% 6.  Neumann (NE) Dirichlet (SW)
% 7.  Neumann (E)  Dirichlet (NSW)
% 8.  Neumann (SE) Dirichlet (NW)
% 9.  Neumann (S)  Dirichlet (NEW)
% 10. Neumann (NSEW)
%
% Uses a second order finite differences on a cell-centered mesh with
% mx-by-my cells in the directions.  Solves the equations using the FFT
% (specifically the fast sine and fast cosine transforms - see
% http://en.wikipedia.org/wiki/Discrete_sine_transform).
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     dfun : the boundary function representing the Dirichlet B.C.
%    nfunx : the boundary function representing the Neumman B.C. in x-direction. 
%    nfuny : the boundary function representing the Neumman B.C. in y-direction. 
%      a,b : the interval defining the square
%       mx : the number of cells in x direction of the mesh.
%       my : the number of cells in y direction of the mesh.
%   bctype : the type of boundary conditions (1-10) as indicated above
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%

hx = (b-a)/mx;   % Mesh spacing
hy = (b-a)/my;

[x,y] = meshgrid(a+hx/2:hx:b-hx/2,a+hy/2:hy:b-hy/2); % Uniform mesh

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(ffun,x,y);

switch bctype
    % Dirichlet (NSEW)
    case 1
        % Compute boundary terms, south, north, east, west
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx);
        jy = (1:my)';
        fhat = dst2(dst2(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idst2(idst2(uhat,1),2);
        
    % Neumann (SW) Dirichlet (NE)
    case 2
        % Compute boundary terms, north, south, east, west
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ubs = feval(nfuny,x(1,:),y(1,:)-hy/2);  
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(nfunx,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) + ubs/hy;          % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) + ubw/hx;          % West

        jx = (1:mx)-1/2;
        jy = (1:my)'-1/2;
        fhat = dct4(dct4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idct4(idct4(uhat,1),2);

    % Neumann (W) Dirichlet (NSE)
    case 3
        % Compute boundary terms, north, south, east, west
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(nfunx,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) + ubw/hx;          % West

        jx = (1:mx)-1/2;
        jy = (1:my)';
        fhat = dct4(dst2(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idct4(idst2(uhat,1),2);

    % Neumann (NW) Dirichlet (SE)
    case 4
        % Compute boundary terms, north, south, east, west
        ubn = feval(nfuny,x(my,:),y(my,:)+hy/2); 
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(nfunx,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - ubn/hy;        % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) + ubw/hx;          % West

        jx = (1:mx)-1/2;
        jy = (1:my)'-1/2;
        fhat = dct4(dst4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idct4(idst4(uhat,1),2);
        
    % Neumann (N) Dirichlet (SEW)
    case 5
        % Compute boundary terms, north, south, east, west
        ubn = feval(nfuny,x(my,:),y(my,:)+hy/2); 
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - ubn/hy;        % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx);
        jy = (1:my)'-1/2;
        fhat = dst2(dst4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idst2(idst4(uhat,1),2);

    % Neumann (NE) Dirichlet (SW)
    case 6
        % Compute boundary terms, south, north, east, west
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ubn = feval(nfuny,x(my,:),y(my,:)+hy/2); 
        ube = feval(nfunx,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - ubn/hy;        % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - ube/hx;        % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx)-1/2;
        jy = (1:my)'-1/2;
        fhat = dst4(dst4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);
 
        uhat = fhat./den;
        u = idst4(idst4(uhat,1),2);

    % Neumann (E) Dirichlet (NSW)
    case 7
        % Compute boundary terms, north, south, east, west
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ubs = feval(dfun,x(1,:),y(1,:)-hy/2);  
        ube = feval(nfunx,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) - 2*ubs/hy^2;      % South
        f(:,mx) = f(:,mx) - ube/hx;        % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx)-1/2;
        jy = (1:my)';
        fhat = dst4(dst2(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idst4(idst2(uhat,1),2);
        
    % Neumann (SE) Dirichlet (NW)
    case 8
        % Compute boundary terms, south, north, east, west
        ubs = feval(nfuny,x(1,:),y(1,:)-hy/2);  
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ube = feval(nfunx,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) + ubs/hy;          % South
        f(:,mx) = f(:,mx) - ube/hx;        % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx)-1/2;
        jy = (1:my)'-1/2;
        fhat = dst4(dct4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);
 
        uhat = fhat./den;
        u = idst4(idct4(uhat,1),2);

    % Neumann (S) Dirichlet (NEW)
    case 9
        % Compute boundary terms, north, south, east, west
        ubn = feval(dfun,x(my,:),y(my,:)+hy/2); 
        ubs = feval(nfuny,x(1,:),y(1,:)-hy/2);  
        ube = feval(dfun,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(dfun,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
        f(1,:) = f(1,:) + ubs/hy;          % South
        f(:,mx) = f(:,mx) - 2*ube/hx^2;    % East
        f(:,1) = f(:,1) - 2*ubw/hx^2;      % West

        jx = (1:mx);
        jy = (1:my)'-1/2;
        fhat = dst2(dct4(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

        uhat = fhat./den;
        u = idst2(idct4(uhat,1),2);
        
    % Neumann (NSEW)
    case 10
        % Compute boundary terms, south, north, east, west
        ubs = feval(nfuny,x(1,:),y(1,:)-hy/2);  
        ubn = feval(nfuny,x(my,:),y(my,:)+hy/2); 
        ube = feval(nfunx,x(:,mx)+hx/2,y(:,mx));   
        ubw = feval(nfunx,x(:,1)-hx/2,y(:,1)); 

        % Adjust f for boundary terms
        f(my,:) = f(my,:) - ubn/hy;    % North
        f(1,:) = f(1,:) + ubs/hy;      % South
        f(:,mx) = f(:,mx) - ube/hx;    % East
        f(:,1) = f(:,1) + ubw/hx;      % West

        jx = (0:mx-1);
        jy = (0:my-1)';
        fhat = dct2(dct2(f,1),2);
        den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);
 
        uhat = fhat./den;
        uhat(1,1) = 0;   % Force the mean of the solution to be zero
        u = idct2(idct2(uhat,1),2);
end        
        
end


