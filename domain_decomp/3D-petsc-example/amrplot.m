function [] = amrplot( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
grid = mmread(filename);
fine_grid = mmread(strcat(filename,'.amr'));
n = size(grid,1);
x = 1:n;
x = x/n;
x = x-1/(2*n);
y = x;
fine_x = 1:2*n;
fine_x = fine_x/(2*n);
fine_x = fine_x-1/(4*n);
fine_y = fine_x;
fine_x = fine_x+1;
figure;
surf(x,y,grid,'EdgeColor','none');
hold on;
surf(fine_x,fine_y,fine_grid,'EdgeColor','none');

end

