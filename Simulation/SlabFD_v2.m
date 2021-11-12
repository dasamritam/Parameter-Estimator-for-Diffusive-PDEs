



function [Ghnum,Rp,y0,profiles] = SlabFD_v2(Chi_e,V_e,tau_e,omega,x,X,Pdep,x_b,x_e,N)

% This function calculates the transfer functions from P -> Te(X) of the PDE:
  % dT/dt = chi*dT^2/dx^2 + V*dT/dx + tau*T + Pdep
%   X = location of the sensors
%   xe = location of the outer edge
%   N = number of discretization points
%   rp = radial coord. of the profile given
%   rpd = radial coord. of the deposition profile
%   chi = diffusion coefficient/profile
%   V  = convection coefficient/profile
%   tau = damping coefficient/profile

%% +++++++++++++++++++ Generate discretization +++++++++

dx = (x_e-x_b)/(N-1);
x = x_b:dx:x_e;

% Check if profiles are given and interpolate

% if length(rp) == 1;
%     Chi_e = chi+0*x;
%     V_e = V+0*x;
%     tau_inv = tau+0*x;
% elseif norm(rp - x)==0;
%     Chi_e = chi;
%     V_e = V;
%     tau_inv = tau;
% else
%     Chi_e = interp1(rp,chi,x,'cubic');
%     V_e = interp1(rp,V,x,'cubic');
%     tau_inv = interp1(rp,tau,x,'cubic');
% end


% interpolate source
% if norm(rpd - x)==0;
%     P = Pdep;
% else P = interp1(rpd,Pdep,x,'cubic');
% end

%% +++++++++++++++++++ A-matrix ++++++++++++++++++++++++
% Generate coefficients    
a_r = Chi_e/(dx^2)+(V_e)/(2*dx);
b_r = -2*Chi_e/(dx^2)-tau_e;
c_r = Chi_e/(dx^2)-(V_e)/(2*dx);

% Generate matrix
A = spdiags([a_r' b_r' c_r'],-1:1,N,N).';

% Boundary conditions: gpunt= (g(x+h)-g(x-h))/2h=0
A(1,2)   = a_r(1)+c_r(1); 

%% Input matrix
B(:,1) = Pdep; % deposition profile
B(N,2) = 0; % Temperature at boundary (not relevant for TF)

%% +++++++++++++++++++ C-matrix ++++++++++++++++++++++++

% Select source and measurement locations
Sensor = interp1(x,x,X,'nearest');
for i = length(Sensor):-1:1;
    indSf(i) = find(x==Sensor(i));
end
if indSf(end) == length(x);
    indSf = indSf(1:end-1);
end

% Measurement locations
Rp = x(indSf); % Measurement positions

C = zeros(length(Sensor),N); j = 1;
for i = indSf 
    C(j,i) = 1;
    j=1+j;
end
C = sparse(C);
% D = zeros(length(vector),2);
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++ Calculate transfer function++++++++++++++++++++

di = spdiags(1i*ones(N,1),0,N,N);
for j = length(omega):-1:1
    sI_A = omega(j)*di-A;
    Ghnum(j,:) = (C/sI_A)*B(:,1);
end

% Steady-state value
y0 = -C/A*B;

profiles = [x(:),Chi_e(:),V_e(:),tau_e(:),Pdep(:)];
end

