clear all;
close all;
clc

mainfolder = cd; % make string with main folder code
addpath(genpath(mainfolder)); % add all subfolders in path

set(0,'defaultaxesfontsize',14);
set(0,'defaultlinelinewidth',1.5);


%% ++++++++++++++++++ Generate measurement data +++++++++++++++++++++++++ %%
%% Simulate temperature measurements

%% Input parameters
omega = 2*pi*25; x_b = 0; x_e = 1; 
Nsim = 500; dxsim = (x_e-x_b)/(Nsim-1);
xsim = x_b:dxsim:x_e; 

% Measurement points
NS = 16; % number of sensors;
DXS = round((x_e-x_b)/(NS)/(dxsim)); % sensor distance
x_meas = xsim(2*DXS:DXS:(NS-2)*DXS);
int_dom = 2*DXS:1:(NS-2)*DXS;

% input power definition
Ptot = 0.7; sigma = 0.1; xdep = 0.5; MW2keVs = 1;%6.24e21/2.1e19; % Factor to convert to keV   
Pdep = 0.2+MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep).^2/sigma.^2);

% Values on grid x
Dsim = 10*xsim.^3-xsim+5;
Vsim = 0*xsim;
Ksim = 0*ones(size(Dsim));
Psim = Pdep;

dt= 1e-4; t = 0:dt:0.04-dt; u = square(2*pi*25*t,70); U = fft(u)/length(u);
f = 1:4;
omega_all = 2*pi*25*f;

[Gh,Rp,y0,profiles] = SlabFD_v2(Dsim,Vsim,Ksim,omega_all,xsim,x_meas,Psim,x_b,x_e,Nsim);
Pomega = U(f+1);
Theta = Gh.*repmat(Pomega(:),1,size(Gh,2));



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ++++++++++++++++++++ Input to the code ++++++++++++++++++++++++++++++ %%

% Define integration grid
options.FD.N = 300; % number of grid-points
[x,dx] = Generate_x_grid(x_meas(1),x_meas(end),options.FD.N);
options.maxNbfOfSteps = 300;

% Define input and output measurements
input  = [Pomega(:),Theta(:,1),Theta(:,end)];
output =  Theta(:,2:end-1);

% Define basis functions
[inp_mat] = Basis_function(x_meas,x,xsim,Dsim,Vsim,Ksim,Psim); % nonsense file needs to be worked out

% Define initial guess by resampling on to grid
Din = interp1(xsim,Dsim,x,'pchip');
Vin = interp1(xsim,Vsim,x,'pchip');
Kin = interp1(xsim,Ksim,x,'pchip');
Pin = interp1(xsim,Psim,x,'pchip');

inp_mat.coeff.D = Din; inp_mat.coeff.V = Vin; inp_mat.coeff.K = Kin; inp_mat.coeff.P = Pin; % Store original profiles

Unknowns = 'DP'; 
inp_mat.cases = Cases_variables_estimated(Unknowns);

%% %%%%%%%%%%%%%%%%%%%% Least-squares scheme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gamma_best, CostBest, CostIdeal, CovP, inp_mat] = LS_PDE(inp_mat,omega,output,input,x_meas,options);

CostBest
CostIdeal

%% %%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1:length(x_meas); indm(jj) = find(x_meas(jj) == xsim); end %% find locations

close all
if inp_mat.cases.D_on
    
   
figure
subplot(3,1,[1,2])
Dest = inp_mat.BP.D*gamma_best.D;
Dinitial = inp_mat.BP.D*inp_mat.gamma_i.D;

hold on
plot(xsim,Dsim,'-','color',0.4*[1,1,1])
plot(x,Dinitial,'k--')
plot(x(1:5:end),Dest(1:5:end),'k.')
plot(x_meas,Dsim(indm),'o','color',0.4*[1,1,1])
plot(x(1:5:end),Dest(1:5:end),'k.')
box on
hold off
grid on
ylabel('$D\left(x\right)$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
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
%l2D = legend('orientation','horizontal','location','northwest');
%set(l2D,'Interpreter','latex')
axis([xsim(1),xsim(end),1.1*min(error_relD),1.1*max(error_relD)])
grid on
box on

end
if inp_mat.cases.V_on
figure
plot(rpd,Vsim,'k',...
     x,inp_mat.BP.V*gamma.V,'k.',...
     x,inp_mat.BP.V*gamma_best.V,'r--')
end

if inp_mat.cases.K_on
subplot(2,2,3)
plot(x,K,'k',...
     x,inp_mat.BP.K*gamma.K,'k.',...
     x,inp_mat.BP.K*gamma_best.K,'r--',...
     c,BPc.K*gamma_best.K,'bo')
end
if inp_mat.cases.P_on
Pest = inp_mat.BP.P*gamma_best.P;
Pinitial = inp_mat.BP.P*inp_mat.gamma_i.P;
figure
subplot(3,1,[1,2])
hold on
plot(xsim,Psim,'-','color',0.4*[1,1,1])
plot(x,Pinitial,'k--')
plot(x(1:5:end),Pest(1:5:end),'k.')
plot(x_meas,Psim(indm),'o','color',0.4*[1,1,1])
plot(x(1:5:end),Pest(1:5:end),'k.')
grid on

hold off
box on
ylabel('$P\left(x\right)$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
l1P = legend('Simulated $P\left(x\right)$',...
             'Initial guess $P\left(x\right)$',...
             'Estimated $P\left(x\right)$',...
             'Sensor locations','location','northwest');
set(l1P,'Interpreter','latex')
subplot(3,1,[3])
Psimerror = interp1(xsim,Psim,x).';
error_relP = (Psimerror-Pest)./Dest;
plot(x,error_relP,'k')
ylabel('$\varepsilon_{rel}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
%l2D = legend('orientation','horizontal','location','northwest');
%set(l2D,'Interpreter','latex')
axis([xsim(1),xsim(end),1.1*min(error_relP),1.1*max(error_relP)])
grid on
box on
end

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

