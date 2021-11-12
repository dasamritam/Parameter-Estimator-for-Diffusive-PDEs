function [Dsim,Vsim,Ksim,Psim,xsim,dxsim,x_meas,omega_all,Input,Output] = data_generator
%% Simulate temperature measurements

% Output variables::

% 1. Dsim,Vsim,Ksim,Psim
% 2. xsim,dxsim,x_meas
% 3. omega_all
% 4. Input
% 5. Output

%% Input parameters
omega = 2*pi*25; x_b = 0; x_e = 1; 
Nsim = 301; dxsim = (x_e-x_b)/(Nsim-1);
xsim = x_b:dxsim:x_e; 

% Measurement points
NS = 10; 
DXS = round((x_e-x_b)/(NS)/(dxsim)); % sensor distance
x_meas = xsim(2*DXS:DXS:(NS-1)*DXS);
% x_meas = round(x_meas,8);
% x_meas = [x_b x_meas x_e];
Sensor_number = length(x_meas); % number of sensors;

% input power definition
Ptot = 0.7; sigma = 0.1; xdep = 0.35; MW2keVs = 1;%6.24e21/2.1e19; % Factor to convert to keV   
Pdep = 0.2+MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep).^2/sigma.^2)...
    +0.8*MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-0.6).^2/sigma.^2);

% Values on grid x
Dsim = 5*xsim.^3-0.005*xsim+5;
Vsim = 15*xsim.^2-0.005;
Ksim = -3*xsim;%xsim*ones(size(Dsim));
Psim = Pdep;




dt= 1e-4; t = 0:dt:0.04-dt; u = square(2*pi*25*t,70); U = fft(u)/length(u);
f = 1:4;
omega_all = 2*pi*25*f;

[Gh,Rp,y0,profiles] = SlabFD_v2(Dsim,Vsim,Ksim,omega_all,xsim,x_meas,Psim,x_b,x_e,Nsim);
Pomega = U(f+1);
Theta = Gh.*repmat(Pomega(:),1,size(Gh,2));

% Define input
Input  = [Pomega(:),Theta(:,1),Theta(:,end)];
Output =  Theta(:,2:end-1);
