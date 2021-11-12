function [x_grid,dx] = Generate_x_grid(x_b,x_e,N)
% Generate grid
%   This m-file assures consistency between grids for basis function and
%   estimations
dx = (x_e-x_b)/(N-1); % N and not (N-1) because N-1 results in edge sensors to be off the grid             
x_grid = x_b:dx:x_e;  
end


% PROBLEM: even or uneven gridpoints N serious errors in relationship to even or
% uneven sensors and number of basis-functions???