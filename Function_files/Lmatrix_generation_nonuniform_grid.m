function L = Lmatrix_generation_nonuniform_grid(N,dx)
% function [L_chi,L_V,L_tau] = Lmatrix_generation(N,dx)
% This function generates the L-matrix components of the finite difference
% approximation dT/dt = (chi*L_chi + V*L_V + L_*tau_inv)*T(x) 
% of the partial differential equation without boundary conditions:
% dT/dt = chi*dT^2/dx^2 + V*dT/dx + tau*T

v = ones(N,1);
w = ((dx.dx_forward-dx.dx_back)./(dx.dx_forward+dx.dx_back)).';
% L_chi construction
L.D = sparse(diag(2./(dx.dx_forward.^2 + dx.dx_back.^2)))*spdiags([v+w -2*v v-w],-1:1,N,N).';

% L_V construction
L.V = sparse(diag(1./(dx.dx_forward + dx.dx_back)))*spdiags([v 0*v -v],-1:1,N,N).';

% L_tau construction
L.K = -spdiags(v,0,N,N).';


end

