function time_matlab()

% ------------------
% Plot Matlab data
% ------------------
N0 = 100;
N1 = 400;
Nvec_matlab = round(logspace(log10(N0),log10(N1),5));

% numbers compare well with G/M paper.
Nvec_matlab = sqrt([40000, 160000, 640000, 2560000]);  

t_matlab_build = zeros(size(Nvec_matlab));
t_matlab_solve = zeros(size(Nvec_matlab));

for i = 1:length(Nvec_matlab)
    N = Nvec_matlab(i);
    
    b = rand(N^2,1);
    
    tic;   
    A = create_matrix(N);
    [L,U,P,Q] = lu(A);
    t_matlab_build(i) = toc;
    
    tic;
    y = L\b;
    x = U\y;
    
    t_matlab_solve(i) = toc;
    
end

save time_matlab_data t_matlab_build t_matlab_solve;

end
