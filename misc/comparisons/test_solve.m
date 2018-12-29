function test_solve()

close all

% Data from Gillman/Martisson paper
N_gaussian = [ 174720,  693504,  2763264,  11031552];  % Number of Gaussian pts.
N_tot      = [1815681, 7252225, 28987905, 115909633];  % Gaussian + Chebyshev

q = 21;          % Number of Gauss points on a side
p = 21;          % Number of Chebyshev nodes 
m = 0.5*(-1 + sqrt(1 + 2*N_gaussian/q));   % Number of leaf boxes on a side
L = log2(m);     % number of levels

% N_cheby = (# boxes)*(interior pts/box + pts on left/bottom edges/box) + 
%         (pts across top edge/domain) + (pts across right edge/domain) +
%         (pt at upper right corner/domain).
N_cheby = ((p-1)*2.^L).^2 + 2*(p-1)*2.^L + 1;
N_check = N_cheby + N_gaussian - N_tot;  
if (N_check ~= 0)
    error('Point count in HPS method is not correct');
end

Nvec_hps = N_cheby;   % Compare to FD discretizations, below.

t_hps_build = [91.68, 371.15, 1661.97, 6894.31];  % Build Ainv operator?
t_hps_solve = [0.34, 1.803, 6.97, 30.67];   % Build Ainv?
t_hps_apply = [0.035,0.104, 0.168, 0.367];  % Apply Ainv? 

t_hps{1} = t_hps_build;
t_hps{2} = t_hps_solve + t_hps_apply;

% ---------------------
% Timing data (matlab)
% ---------------------

load time_matlab_data;     % t_matlab_build, t_matlab_solve

Nvec_matlab = [40000, 160000, 640000, 2560000];    % == N^2

t_matlab{1} = t_matlab_build;
t_matlab{2} = t_matlab_solve;

t_matlab_gm{1} = [2.46e-1, 1.29, 6.87, 49.86];
t_matlab_gm{2} = [5.32e-3, 2.74e-2, 1.33e-1, 6.98e-1]; 


% ---------------
% From CBMS code
% ---------------

% nleaf = 5
Nvec_smatrix = [1600, 6400, 25600, 102400, 409600, 1638400];  % == N^2

t_smatrix{1} = [0.06, 0.09, 0.58, 5.11, 77.1, 1861.73];
t_smatrix{2} = [0.01, 0.02, 0.07, 0.26,  1.03,    6.50];

%{
% nleaf = 8
Nvec_smatrix = [4096, 16384, 65536, 262144, 1048576];
t_smatrix{1} = [0.04, 0.20, 1.51, 19.63, 313.64];
t_smatrix{2} = [0.01, 0.02, 0.09,  0.35,   1.50];
%}

% nleaf = 16; L = 2,3,4,5,6,7
% Memory for last solve : 19574421/1024^2 = 18GB - so timing is
% terrible.
Nvec_smatrix = [4096, 16384, 65536, 262144, 1048576, 4194304];
t_smatrix{1} = [0.07, 0.23, 1.19, 9.95, 119.48, 2810.22];
t_smatrix{2} = [0.01, 0.03, 0.05, 0.19, 0.80, 177.75];

%{
% nleaf incerases; keep level fixed at 4
Nvec_smatrix = [4096, 16384, 65536, 262144, 1048576];
t_smatrix{1} = [0.13, 0.21, 1.19, 11.86, 394.56];
t_smatrix{2} = [0.03, 0.03, 0.05,  0.30, 304.16];
%}


% ---------------
% Plot build data16384
% ---------------

Nvec_all = sort([Nvec_matlab Nvec_hps Nvec_smatrix 2e9]);


hd = zeros(1,4);
for i = 1:2   % Build and solve
    figure(i);
    clf
    
    % HPS
    loglog(Nvec_hps,t_hps{i},'k.','markersize',30);
    p_hps{i} = polyfit(log(Nvec_hps),log(t_hps{i}),1);
    hold on;
    hd(1) = loglog(Nvec_all,exp(polyval(p_hps{i},log(Nvec_all))),'r-','linewidth',2);
    lstr{1} = sprintf('HPS (slope = %5.2f)',p_hps{i}(1));


    % Matlab
    hold on;
    p_lu{i} = polyfit(log(Nvec_matlab),log(t_matlab{i}),1);
    hd(2) = loglog(Nvec_all,exp(polyval(p_lu{i},log(Nvec_all))),'b','linewidth',2);
    lstr{2} = sprintf('LU (slope = %5.2f)',p_lu{i}(1));

    % Smatrix (from CBMS code main01.m)    
    idstart = 1;
    idend = 5;
    loglog(Nvec_smatrix,t_smatrix{i},'ok','markersize',15,'linewidth',2);
    p_smatrix{i} = polyfit(log(Nvec_smatrix(idstart:5)),...
        log(t_smatrix{i}(idstart:5)),1);
    hd(3) = loglog(Nvec_all,exp(polyval(p_smatrix{i},log(Nvec_all))),...
        'm','linewidth',2);
    lstr{3} = sprintf('S-matrix (slope = %5.2f)',p_smatrix{i}(1));
    
    % Matlab data
    hd(4) = loglog(Nvec_matlab,t_matlab{i},'pk','markersize',15,'linewidth',2);
    lstr{4} = 'LU (matlab)';

    
    
    % Reference data
    hd(5) = loglog(Nvec_matlab, t_matlab_gm{i},'ro','markersize',15,'linewidth',2);
    lstr{5} = 'LU (GM)';
        
    
    lh = legend(hd,lstr);
    set(lh,'fontsize',16,'location','northwest');
    set(lh,'AutoUpdate','off');

    set(gca,'fontsize',16);
end

figure(1)
title('Build time','fontsize',18);
xlabel('N','fontsize',16);
ylabel('t','fontsize',16);

figure(2)
title('Solve time','fontsize',18);
xlabel('N','fontsize',16);
ylabel('t (s)','fontsize',16);


shg


% ------------------
% Print out results
% ------------------

set(lh,'AutoUpdate','off');
levels = L;

l = double('-')*ones(1,55);
header_str = sprintf('%16s %12s %12s %12s','N (L)','HPS','LU (est)','S-matrix');

%%
% Build times
fprintf('%s\n',l);
fprintf('%s (q = p = %d; N = # leaf boxes)\n','Build times',q);
fprintf('%s\n',header_str);
fprintf('%s\n',l);
for i = 1:length(N_gaussian)
    N = Nvec_hps(i);
    lu_est = exp(polyval(p_lu{1},log(N)));
    smatrix_est = exp(polyval(p_smatrix{1},log(N)));
    fprintf('%12d (%d) %12.2f %12.2f %12.2f\n',N,L(i),t_hps_build(i),lu_est,smatrix_est);
    figure(1);
    hold on;
%     loglog(N,lu_est,'kp','markersize',15);
end
fprintf('%s\n',l);

%%
% Solve times
fprintf('%s\n',l);
fprintf('%s (q = p = %d; N = # leaf boxes)\n','Solve times',q);
fprintf('%s\n',header_str);
fprintf('%s\n',l);
for i = 1:length(N_gaussian)
    N = Nvec_hps(i);
    lu_est = exp(polyval(p_lu{2},log(N)));
    smatrix_est = exp(polyval(p_smatrix{2},log(N)));
    fprintf('%12d (%d) %12.2f %12.2f %12.2f\n',N,L(i),t_hps_solve(i),lu_est,smatrix_est);
    figure(2);
    hold on;
%     loglog(N,lu_est,'kp','markersize',15);
end
fprintf('%s\n',l);

end



  

