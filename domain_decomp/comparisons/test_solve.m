function test_solve()

close all

% Data from Gillman/Martisson paper
N_gaussian = [174720, 693504, 2763264, 11031552];  % Number of Gaussian pts.

q = 21;                 % Number of Gauss points on a side
p = 0.5*(-1 + sqrt(1 + 2*N_gaussian/q));   % Number of leaf boxes on a side
L = log2(p);     % number of levels
Nvec_hps = p.*2.^L;  % Number of total points on a side.

t_hps_build = [91.68, 371.15, 1661.97, 6894.31];
t_hps_solve = [0.34, 1.803, 6.97, 30.67];

t_hps{1} = t_hps_build;
t_hps{2} = t_hps_solve;



% ---------------------
% Timing data (matlab)
% ---------------------

load time_matlab_data;     % t_matlab_build, t_matlab_solve

Nvec_matlab = sqrt([40000, 160000, 640000, 2560000]);

t_matlab{1} = t_matlab_build;
t_matlab{2} = t_matlab_solve;

t_matlab_gm{1} = [2.46e-1, 1.29, 6.87, 49.86];
t_matlab_gm{2} = [5.32e-3, 2.74e-2, 1.33e-1, 6.98e-1]; 

% ---------------
% Plot build data
% ---------------

hd = zeros(1,4);
for i = 1:2
    figure(i);
    clf
    loglog(Nvec_hps,t_hps{i},'k.','markersize',30);
    p_hps{i} = polyfit(log(Nvec_hps),log(t_hps{i}),1);
    hold on;
    hd(1) = loglog(Nvec_hps,exp(polyval(p_hps{i},log(Nvec_hps))),'r-','linewidth',2);
    lstr{1} = sprintf('HPS (slope = %5.2f)',p_hps{i}(1));


    % Matlab
    hd(3) = loglog(Nvec_matlab,t_matlab{i},'pk','markersize',15,'linewidth',2);
    hold on;
    lstr{3} = 'LU (matlab)';

    p_lu{i} = polyfit(log(Nvec_matlab),log(t_matlab{i}),1);

    Nv2 = [Nvec_matlab Nvec_hps];
    hd(2) = loglog(Nv2,exp(polyval(p_lu{i},log(Nv2))),'b','linewidth',2);
    lstr{2} = sprintf('LU (slope = %5.2f)',p_lu{i}(1));
    
    % Reference data
    hd(4) = loglog(Nvec_matlab, t_matlab_gm{i},'ro','markersize',15,'linewidth',2);
    lstr{4} = 'LU (GM)';
        
    
    lh = legend(hd,lstr);
    set(lh,'fontsize',16,'location','northwest');

    set(gca,'fontsize',16);
end

figure(1)
title('Build time','fontsize',18);
xlabel('N (leaf boxes per side)','fontsize',16);
ylabel('t','fontsize',16);

figure(2)
title('Solve time','fontsize',18);
xlabel('N (leaf boxes per side)','fontsize',16);
ylabel('t (s)','fontsize',16);


shg


% ------------------
% Print out results
% ------------------

set(lh,'AutoUpdate','off');
levels = L;

l = double('-')*ones(1,38);
header_str = sprintf('%12s %12s %12s','N (L)','HPS','LU (est)');

%%
% Build times
fprintf('%s\n',l);
fprintf('%s (q = %d; N = # leaf boxes)\n','Build times',q);
fprintf('%s\n',header_str);
fprintf('%s\n',l);
for i = 1:length(N_gaussian)
    N = Nvec_hps(i);
    lu_est = exp(polyval(p_lu{1},log(N)));
    fprintf('%8d (%d) %12.2f %12.2f\n',N,L(i),t_hps_build(i),lu_est);
    figure(1);
    hold on;
    loglog(N,lu_est,'kp','markersize',15);
end
fprintf('%s\n',l);

%%
% Solve times
fprintf('%s\n',l);
fprintf('%s (q = %d; N = # leaf boxes)\n','Solve times',q);
fprintf('%s\n',header_str);
fprintf('%s\n',l);
for i = 1:length(N_gaussian)
    N = Nvec_hps(i);
    lu_est = exp(polyval(p_lu{2},log(N)));
    fprintf('%8d (%d) %12.2f %12.2f\n',N,L(i),t_hps_solve(i),lu_est);
    figure(2);
    hold on;
    loglog(N,lu_est,'kp','markersize',15);
end
fprintf('%s\n',l);

end



  

