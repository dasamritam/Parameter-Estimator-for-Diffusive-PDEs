function L = Lmatrix_generation(N,dx)
% function [L_chi,L_V,L_tau] = Lmatrix_generation(N,dx)
% This function generates the L-matrix components of the finite difference
% approximation dT/dt = (chi*L_chi + V*L_V + L_*tau_inv)*T(x) 
% of the partial differential equation without boundary conditions:
% dT/dt = chi*dT^2/dx^2 + V*dT/dx + tau*T

v = ones(N,1);

% L_chi construction
L.D = (1/dx^2)*spdiags([v -2*v v],-1:1,N,N).';

% L_V construction
L.V = (0.5/dx)*spdiags([v 0*v -v],-1:1,N,N).';

% L_tau construction
L.K = -spdiags(v,0,N,N).';


end

