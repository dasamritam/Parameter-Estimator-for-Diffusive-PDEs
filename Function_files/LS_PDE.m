function [gamma_best, CostBest, CostIdeal, CovP, inp_mat] = LS_PDE(inp_mat,omega,output,input,x_meas,options)
% Generate non-recursive components of optimization algorithm
%   Detailed explanation goes here


%% Define which gamma's need to be estimated
if inp_mat.cases.D_on; gamma.D = inp_mat.gamma_i.D; else gamma.D = []; end
if inp_mat.cases.V_on; gamma.V = inp_mat.gamma_i.V; else gamma.V = []; end
if inp_mat.cases.K_on; gamma.K = inp_mat.gamma_i.K; else gamma.K = []; end
if inp_mat.cases.P_on; gamma.P = inp_mat.gamma_i.P; else gamma.P = []; end


%% Build discretization grid
[x_grid,dx] = Generate_x_grid(x_meas(1),x_meas(end),options.FD.N);

%% Generate non-recursive (permanent) state-space matrices
v_bc = 2:length(x_grid)-1; % determines size of A, B, C matrices
lv_bc = length(v_bc);

% Sensor locations on grid
[C,S] = C_matrix_generation(x_grid,x_meas(2:end-1),options.FD.N); % warning sensors are moved to grid points
inp_mat.C = C(:,v_bc); % define C-matrix

% Make complex vector s=0+1i
inp_mat.di = spdiags(1i*ones(lv_bc,1),0,lv_bc,lv_bc);

% Define Laplacian matrix in terms of D, V, and K
inp_mat.L = Lmatrix_generation(lv_bc,dx);

% store variables
inp_mat.dx = dx;
inp_mat.v_bc = v_bc;
inp_mat.S = S;

%% Optimizations algorithm
[gamma_best, CostBest, CovP] = LM_optimization(inp_mat,gamma,omega,output,input,options);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Check cost function for simulated profiles
inp_mat_ideal = inp_mat; inp_mat_ideal.cases.D_on = 0; inp_mat.cases.V_on = 0;
inp_mat_ideal.cases.K_on = 0; inp_mat_ideal.cases.P_on = 0; 

[e_ideal] = Cost_function_Jacobian_ideal(inp_mat_ideal,gamma,omega,output,input);
CostIdeal = e_ideal'*e_ideal;



end

