function [inp_mat] = Basis_function(x_meas,x,xsim,Dsim,Vsim,Ksim,Psim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Generate B-spline for basis function

degree  = 3;
% choose control points B-spline
c = linspace(x_meas(1),x_meas(end),9);

% Generate intial guesss
Dc = interp1(xsim,Dsim,c,'pchip');
Vc = interp1(xsim,Vsim,c,'pchip');
Kc = interp1(xsim,Ksim,c,'pchip');
Pc = interp1(xsim,Psim,c,'pchip');

Dc = ones(size(interp1(xsim,Dsim,c,'pchip')));
Pc = ones(size(interp1(xsim,Psim,c,'pchip')))*5;

[inp_mat.BP.D,inp_mat.gamma_i.D,~] = Bspline_generation(c,Dc,x,degree); % Diffusion coefficient
[inp_mat.BP.V,inp_mat.gamma_i.V,~] = Bspline_generation(c,Vc,x,degree); % Convective velocity
[inp_mat.BP.K,inp_mat.gamma_i.K,~] = Bspline_generation(c,Kc,x,degree); % Damping coefficient K_inv
[inp_mat.BP.P,inp_mat.gamma_i.P,~] = Bspline_generation(c,Pc,x,degree); % Power deposition profile

%% Generate polynomial approximation

Bpoly = [x.^5;x.^4;x.^3;x.^2;x.^1;x.^0]; inp_mat.BP.D = Bpoly.'; inp_mat.gamma_i.D = 1*ones(size(inp_mat.BP.D,2),1);
inp_mat.gamma_i.D = [0;0;0;0;0;1];

% Bpoly = [x.^1;x.^0]; inp_mat.BP.D = Bpoly.'; inp_mat.gamma_i.D = 10*ones(size(inp_mat.BP.D,2),1);
% inp_mat.gamma_i.D = [0;10];
%% ++++++++++++++++++ Finite difference++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end

