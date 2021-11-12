% Author: Amritam Das
% Last Upadte: 17-09-2018

%% VARIABLES WE NEED TO CHOOSE/ TUNE. NO SYSTEMATIC WAY TO DESIGN THE ALGO. IDEPENDENT OF THESE.  
% 1. Number of sensors 
% 2. Number of control points in bspline
% 3. Choice of Basis
% 4. Number of Basis

%% The code still suffers when the grid for Data Generator and estimator model are chosen to be different.
%% Some how the assymetric gridding does not solve the problem completely.
%% Assymetric griding does make the result and adjustment in the jacobian does improve the results. 

clear;
close all;
clc

mainfolder = cd; % make string with main folder code
addpath(genpath(mainfolder)); % add all subfolders in path

set(0,'defaultaxesfontsize',14);
set(0,'defaultlinelinewidth',1.5);
%% ++++++++++++++++++++++++ Assignment of the logic arguments +++++++++++++++++++++++++++ %%

% Use the same grid from the data generator

simulation_grid_YES = 0;
orthonormalization_YES = 0;


%% ++++++++++++++++++ Generate measurement data +++++++++++++++++++++++++ %%

% This generates the required ip/op data for the estimation.

[Dsim,Vsim,Ksim,Psim,xsim,dxsim,x_meas,omega_all,Input,Output] = data_generator;

%% ++++++++++++++++++++ Setting up the Estimation Problem ++++++++++++++++++++++++++++++ %

%% Part A:  Descrtization grid for the MODEL.

% Intial Grid with extremum sensors

[x_old, x_b, x_e, N] = initial_model_grid(x_meas, dxsim, simulation_grid_YES);

% Include the missing sensors in the grid

[x, dx, N,miss] = add_missing_sensors(x_old, x_meas);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate B-spline for basis function BF

degree  = 3;
% choose control points B-spline
c = linspace(x_b,x_e,13);

% Generate intial guesss
Dc = ones(size(interp1(xsim,Dsim,c,'pchip')));
Pc = ones(size(interp1(xsim,Psim,c,'pchip')));
Vc = ones(size(interp1(xsim,Vsim,c,'pchip')));
Kc = ones(size(interp1(xsim,Ksim,c,'pchip')));

[inp_mat.BP.D_nonorth,inp_mat.gamma_i.D_nonorth,ass.D] = Bspline_generation(c,Dc,x,degree); % Diffusion coefficient
[inp_mat.BP.V_nonorth,inp_mat.gamma_i.V_nonorth,ass.V] = Bspline_generation(c,Vc,x,degree); % Convective velocity
[inp_mat.BP.K_nonorth,inp_mat.gamma_i.K_nonorth,ass.K] = Bspline_generation(c,Kc,x,degree); % Damping coefficient K_inv
[inp_mat.BP.P_nonorth,inp_mat.gamma_i.P_nonorth,ass.P] = Bspline_generation(c,Pc,x,degree); % Power deposition profile

%% Generate polynomial approximation

Dpoly = [x.^7;x.^6;x.^5;x.^4;x.^3;x.^2;x.^1;x.^0]; inp_mat.BP.D_nonorth = Dpoly.'; inp_mat.gamma_i.D_nonorth = 0.1*ones(size(inp_mat.BP.D_nonorth,2),1);
inp_mat.gamma_i.D_nonorth = [0;0;0;0;0;0;0;1];

Vpoly = [x.^5;x.^4;x.^3;x.^2;x.^1;x.^0]; inp_mat.BP.V_nonorth = Vpoly.'; inp_mat.gamma_i.V_nonorth = 0.1*ones(size(inp_mat.BP.V_nonorth,2),1);
inp_mat.gamma_i.V_nonorth = [0;0;0;0;0;1];

Kpoly = [x.^2;x.^1;x.^0]; inp_mat.BP.K_nonorth = Kpoly.'; inp_mat.gamma_i.K_nonorth = 0.1*ones(size(inp_mat.BP.K_nonorth,2),1);
inp_mat.gamma_i.K_nonorth = [0;0;1];

%% Orthonormalization

if orthonormalization_YES == 1
    inp_mat.BP.D = orthonormalization(inp_mat.BP.D_nonorth,Dc); 
    inp_mat.BP.V = orthonormalization(inp_mat.BP.V_nonorth,Vc); 
    inp_mat.BP.K = orthonormalization(inp_mat.BP.K_nonorth,Kc);
    inp_mat.BP.P = orthonormalization(inp_mat.BP.P_nonorth,Pc);
else
    inp_mat.BP.D = inp_mat.BP.D_nonorth; 
    inp_mat.BP.V = inp_mat.BP.V_nonorth; 
    inp_mat.BP.K = inp_mat.BP.K_nonorth;
    inp_mat.BP.P = inp_mat.BP.P_nonorth;
end

inp_mat.gamma_i.D = inp_mat.gamma_i.D_nonorth; 
inp_mat.gamma_i.V = inp_mat.gamma_i.V_nonorth; 
inp_mat.gamma_i.K = inp_mat.gamma_i.K_nonorth;
inp_mat.gamma_i.P = inp_mat.gamma_i.P_nonorth;
%% ++++++++++++++++++ Finite difference++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate intial guesss
Din = interp1(xsim,Dsim,x,'pchip');
Vin = interp1(xsim,Vsim,x,'pchip');
Kin = interp1(xsim,Ksim,x,'pchip');
Pin = interp1(xsim,Psim,x,'pchip');

% Store original profiles
inp_mat.coeff.D = Din; inp_mat.coeff.V = Vin; inp_mat.coeff.K = Kin; inp_mat.coeff.P = Pin;


% Boundary conditions are measurements

xbc_meas = x_meas(2:end-1);

xdom = [x_meas(1),xbc_meas,x_meas(end)];

% 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting the Model for Estimation

Unknowns = 'DVKP';


v_bc = 2:N-1; % determines size of A, B, C matrices
lv_bc = length(v_bc);
%% Sensor locations on grid
[C,S] = C_matrix_generation(x,xbc_meas,N); % this is questionable
inp_mat.C = C(:,v_bc); 
%% Make complex vector
inp_mat.di = spdiags(1i*ones(lv_bc,1),0,lv_bc,lv_bc);

%% Define Laplacian matrix in terms of D, V, and K_inv
inp_mat.L = Lmatrix_generation_nonuniform_grid(lv_bc,dx);

inp_mat.cases = Cases_variables_estimated(Unknowns);
% ==================================================================
% ==================================================================
%% +++++++++++++++++++ Algorithm starts here +++++++++++++++++++
if inp_mat.cases.D_on; gamma.D = inp_mat.gamma_i.D;
else gamma.D = []; end
    
if inp_mat.cases.V_on; gamma.V = inp_mat.gamma_i.V;
else gamma.V = []; end

if inp_mat.cases.K_on; gamma.K = inp_mat.gamma_i.K;
else gamma.K = []; end

if inp_mat.cases.P_on; gamma.P = inp_mat.gamma_i.P;
else gamma.P = []; end


%% Select parameters to be estimated


inp_mat.dx = dx;
inp_mat.v_bc = v_bc;




%% Least-squares scheme
options.maxNbfOfSteps = 300;

[gamma_best, CostBest, CovP] = LM_optimization(inp_mat,gamma,omega_all,Output,Input,options);




inp_mat_ideal = inp_mat; inp_mat_ideal.cases.D_on = 0; inp_mat_ideal.cases.V_on = 0;
inp_mat_ideal.cases.K_on = 0; inp_mat_ideal.cases.P_on = 0; 


gamma_ideal.D = []; gamma_ideal.P = [];gamma_ideal.V = [];gamma_ideal.K = [];
[e] = Cost_DVKP_Bspline_ideal(inp_mat_ideal,gamma_ideal,omega_all,Output,Input,xdom);


CostBest
CostIdeal = e'*e

%%
for jj = 1:length(x_meas); indm(jj) = find(x_meas(jj) == xsim); end %% find locations

%% ++++++++++++++++++++++++++++++++++ Visualization ++++++++++++++++++++ %%
close all;

set(gcf,'Renderer','painters');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',20,'interpreter','latex')
set(findall(figureHandle,'type','axes'),'Color','w','fontSize',20)
set(findall(figureHandle,'type','axes'),'TickLabelInterpreter','latex')

if inp_mat.cases.D_on  
subplot(3,1,[1,2])
Dest = inp_mat.BP.D*gamma_best.D;
Dinitial = inp_mat.BP.D*gamma.D;

hold on
plot(xsim,Dsim,'-','color',0.4*[1,1,1])
plot(x,Dinitial,'k--')
plot(x(1:5:end),Dest(1:5:end),'r+')
stem(x_meas,Dsim(indm),'o','color',[0,0,1])
box on
hold off
grid on
ylabel('$D\left(x\right)$','Interpreter','latex')
l1D = legend('Simulated $D\left(x\right)$',...
             'Initial guess $D\left(x\right)$',...
             'Estimated $D\left(x\right)$',...
             'Sensor locations','location','northwest');
set(l1D,'Interpreter','latex')

subplot(3,1,[3])
hold on
Dsimerror = interp1(xsim,Dsim,x).';
error_relD = (Dsimerror-Dest)./Dest;
plot(x,error_relD,'k')
hold off
ylabel('$\varepsilon_{rel}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
grid on
box on

end


%%
if inp_mat.cases.V_on
figure
subplot(3,1,[1,2])
Vest = inp_mat.BP.V*gamma_best.V;
Vinitial = inp_mat.BP.V*gamma.V;

hold on
plot(xsim,Vsim,'-','color',0.4*[1,1,1])
plot(x,Vinitial,'k--')
plot(x(1:5:end),Vest(1:5:end),'r+')
stem(x_meas,Vsim(indm),'o','color',[0,0,1])
box on
hold off
grid on
ylabel('$U\left(x\right)$','Interpreter','latex')
l1D = legend('Simulated $U\left(x\right)$',...
             'Initial guess $U\left(x\right)$',...
             'Estimated $U\left(x\right)$',...
             'Sensor locations','location','northwest');
set(l1D,'Interpreter','latex')
axis([xsim(1),xsim(end),min(Vsim),1.1*max(Vsim)])
subplot(3,1,[3])
hold on
Vsimerror = interp1(xsim,Vsim,x).';
error_relV = (Vsimerror-Vest)./Vest;
plot(x,error_relV,'k')
hold off
ylabel('$\varepsilon_{rel}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
axis([xsim(1),xsim(end),1.1*min(error_relV),1.1*max(error_relV)])
grid on
box on

end

% % savePDF(figureHandle.Position([3]),figureHandle.Position([4]))
%%
if inp_mat.cases.K_on
figure
subplot(3,1,[1,2])
Kest = inp_mat.BP.K*gamma_best.K;
Kinitial = inp_mat.BP.K*gamma.K;

hold on
plot(xsim,Ksim,'-','color',0.4*[1,1,1])
plot(x,Kinitial,'k--')
plot(x(1:5:end),Kest(1:5:end),'r+')
h = stem(x_meas,Ksim(indm),'o','color',[0,0,1]);
set(h,'BaseValue',min(Ksim));
box on
hold off
grid on
ylabel('$K\left(x\right)$','Interpreter','latex')
l1D = legend('Simulated $K\left(x\right)$',...
             'Initial guess $K\left(x\right)$',...
             'Estimated $K\left(x\right)$',...
             'Sensor locations','location','northwest');
set(l1D,'Interpreter','latex')

subplot(3,1,[3])
hold on
Ksimerror = interp1(xsim,Ksim,x).';
error_relK = (Ksimerror-Kest)./Kest;
plot(x,error_relK,'k')
hold off
ylabel('$\varepsilon_{rel}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
axis([xsim(1),xsim(end),1.1*min(error_relK),1.1*max(error_relK)])
grid on
box on

end


%%

if inp_mat.cases.P_on
Pest = inp_mat.BP.P*gamma_best.P;
Pinitial = inp_mat.BP.P*gamma.P;
figure
subplot(3,1,[1,2])
hold on
plot(xsim,Psim,'-','color',0.4*[1,1,1])
plot(x,Pinitial,'k--')
plot(x(1:5:end),Pest(1:5:end),'r+')
stem(x_meas,Psim(indm),'o','color',[0,0,1])
grid on

hold off
box on
ylabel('$P\left(x\right)$','Interpreter','latex')
l1P = legend('Simulated $P\left(x\right)$',...
             'Initial guess $P\left(x\right)$',...
             'Estimated $P\left(x\right)$',...
             'Sensor locations'); %,'location','northwest'
set(l1P,'Interpreter','latex')
subplot(3,1,[3])
Psimerror = interp1(xsim,Psim,x).';
error_relP = (Psimerror-Pest)./Pest;
plot(x,error_relP,'k')
ylabel('$\varepsilon_{rel}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
axis([xsim(1),xsim(end),1.1*min(error_relP),1.1*max(error_relP)])
grid on
box on
end
%% 

% % savePDF(figureHandle.Position([3]),figureHandle.Position([4]))
return

%% ++++++++++++++++++++++++++ Test Jacobian +++++++++++++++++++++++++++ %%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

[Gout,Jout] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);

dtheta = 1e-8; 

%% Diffusion
if inp_mat.cases.D_on;
Jnum1D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,1);
Jnum2D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,2);
Jnum3D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,3);
JD_error = [dtheta,norm(Jnum1D-Jout.Jd(:,:,1)),norm(Jnum2D-Jout.Jd(:,:,2)),norm(Jnum3D-Jout.Jd(:,:,3))]
end
%% Convection
if inp_mat.cases.V_on; 
Jnum1V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,1);
Jnum2V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,2);
Jnum3V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,3);
JV_error = [dtheta,norm(Jnum1V-Jout.Jv(:,:,1)),norm(Jnum2V-Jout.Jv(:,:,2)),norm(Jnum3V-Jout.Jv(:,:,3))]
end
%% Damping
if inp_mat.cases.K_on; 
Jnum1K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,1);
Jnum2K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,2);
Jnum3K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,3);
JK_error = [dtheta,norm(Jnum1K-Jout.Jk(:,:,1)),norm(Jnum2K-Jout.Jk(:,:,2)),norm(Jnum3K-Jout.Jk(:,:,3))]
end
%% P
if inp_mat.cases.P_on; 
JnumP = Test_Jacobian_P(inp_mat,gamma,omega,dtheta,1);
JP_error = [dtheta,norm(JnumP-Jout.Jp(:,:,1))]
end

%% % Test full Jacobian
[E,JE] =  Cost_function_Jacobian(inp_mat,gamma,omega_all,Output,[Input,0*Input,0*Input],xdom);

[Jnum] = Test_Jacobian_full(inp_mat,gamma,omega_all,Output,[Input,0*Input,0*Input],xdom,dtheta);

Jfull_error = norm(JE-Jnum)

%% Figure presentation
figure(123)
hold on
plot(xsim,inp_mat.BP.P*inp_mat.gamma_i.P,'k.-')
for i = 1:length(gamma.P);
    plot(xsim,inp_mat.BP.P(:,i).*gamma.P(i))
end
plot(xsim,inp_mat.BP.P*inp_mat.gamma_i.P,'k--')
hold off
box on
grid on
xlabel('x')
ylabel('P(x)')
legend('P(x_i)','B-splines')

