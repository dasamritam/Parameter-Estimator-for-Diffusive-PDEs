function [G,J] = G_and_J_G_nonuniform_grid(inp_mat,gamma,omega)
% Calculates transfer function and corresponding Jacobian of transport
% model where chi has been approximated using B-splines
%   L is the laplacian matrix
%   B is the B matrix of a state-space model
%   C is the C matrix of a state-space model
%   BP is B-spline matrix
%   gamma are the coefficients to de adjusted
%   omega is the frequency to be considered
%   di is a predefined matrix
% %%%%
%   G transfer function
%   J dG/dgamma_i Jacobian
v_bc = inp_mat.v_bc;

%% Check which variables should be iterated

if inp_mat.cases.D_on 
   D = inp_mat.BP.D*gamma.D;
else D = inp_mat.coeff.D; 
end
    
if inp_mat.cases.V_on 
   V = inp_mat.BP.V*gamma.V;
else V = inp_mat.coeff.V; 
end

if inp_mat.cases.K_on 
   K = inp_mat.BP.K*gamma.K;
else K = inp_mat.coeff.K; 
end

if inp_mat.cases.P_on 
   P = inp_mat.BP.P*gamma.P;
else P = inp_mat.coeff.P; 
end

%% Construct A matrix from components
A = diag(D(v_bc))*inp_mat.L.D...
  + diag(V(v_bc))*inp_mat.L.V...
  + diag(K(v_bc))*inp_mat.L.K;

%% Construct B matrix from components
B(:,1) = P(v_bc);

B(1,2) = 2*D(2)/(inp_mat.dx.dx_forward(1)^2 + inp_mat.dx.dx_back(1)^2)...
    *(1+ (inp_mat.dx.dx_forward(1)-inp_mat.dx.dx_back(1))/(inp_mat.dx.dx_forward(1)+inp_mat.dx.dx_back(1)))...
    -(1/(inp_mat.dx.dx_forward(1)+inp_mat.dx.dx_back(1))) *V(2);

B(end,3) = 2*D(end-1)/(inp_mat.dx.dx_forward(end)^2 + inp_mat.dx.dx_back(end)^2)...
    *(1- (inp_mat.dx.dx_forward(end)-inp_mat.dx.dx_back(end))/(inp_mat.dx.dx_forward(end)+inp_mat.dx.dx_back(end)))...
    +(1/(inp_mat.dx.dx_forward(end)+inp_mat.dx.dx_back(end)))*V(end-1);

% B(1,2) = D(1)/inp_mat.dx.^2 + 0.5/inp_mat.dx*V(1); % wrong second row not first row
% B(end,3) = D(end)/inp_mat.dx.^2 - 0.5/inp_mat.dx*V(end); % not expected like this

%% Calculate transfer function
sI_A = omega*inp_mat.di-A; %(sI-A)

C_inv_sI_A = inp_mat.C/sI_A; % C*inv(sI-A)

G = C_inv_sI_A*B; % Calculate transfer function

%% +++++++++++ Calculate corresponding Jacobian +++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Predefine inverse matrix calculations
AtotsB = sI_A\B; % inv(sI-A)*B


%% +++++++++++++++++ D derivative +++++++++++++++++++++++++
if inp_mat.cases.D_on
    
    Linv_D = inp_mat.L.D*AtotsB;    % L_D*inv(sI-A)*B

    for i = length(gamma.D):-1:1 % this could maybe improved
        Jd(:,i,:) = C_inv_sI_A*diag(inp_mat.BP.D(v_bc,i))*Linv_D;
    end

    % if b is boundary conidtions
    
    % derivative of B-matrix column 2   
    J2bd = C_inv_sI_A(:,1)*inp_mat.BP.D(2,:).*2/(inp_mat.dx.dx_forward(1)^2 + inp_mat.dx.dx_back(1)^2)...
    .*(1+ (inp_mat.dx.dx_forward(1)-inp_mat.dx.dx_back(1))/(inp_mat.dx.dx_forward(1)+inp_mat.dx.dx_back(1)));

     % derivative of B-matrix column 3
    J3bd = C_inv_sI_A(:,end)*inp_mat.BP.D(end-1,:).*2/(inp_mat.dx.dx_forward(end)^2 + inp_mat.dx.dx_back(end)^2)...
    *(1- (inp_mat.dx.dx_forward(end)-inp_mat.dx.dx_back(end))/(inp_mat.dx.dx_forward(end)+inp_mat.dx.dx_back(end)));

    Jd(:,:,2) = Jd(:,:,2) + J2bd; % Combine derivatives 
    Jd(:,:,3) = Jd(:,:,3) + J3bd; % Combine derivatives
    
    J.Jd = Jd;
else Jd = [];
end

%% +++++++++++++++++ V derivative +++++++++++++++++++++++++
if inp_mat.cases.V_on

    Linv_V = inp_mat.L.V*AtotsB;    % L_V*inv(sI-A)*B

    for i = length(gamma.V):-1:1 % this could maybe improved
        Jv(:,i,:) = C_inv_sI_A*diag(inp_mat.BP.V(v_bc,i))*Linv_V;
    end
    
    % if b is boundary conditions
    
    % derivative of B-matrix column 2 
    J2bv = C_inv_sI_A(:,1)*-1*inp_mat.BP.V(2,:)*(1/(inp_mat.dx.dx_forward(1)+inp_mat.dx.dx_back(1)));  
    
    % derivative of B-matrix column 3
    J3bv = C_inv_sI_A(:,end)*1*inp_mat.BP.V(end-1,:)*(1/(inp_mat.dx.dx_forward(end)+inp_mat.dx.dx_back(end))); 

    Jv(:,:,2) = Jv(:,:,2) + J2bv; % Combine derivatives 
    Jv(:,:,3) = Jv(:,:,3) + J3bv; % Combine derivatives
    
    J.Jv = Jv;
else Jv = [];
end
%% +++++++++++++++++ K derivative +++++++++++++++++++++++++
if inp_mat.cases.K_on
    
    Linv_K = inp_mat.L.K*AtotsB;    % % L_V*inv(sI-A)*B

    for i = length(gamma.K):-1:1 % this could maybe improved
        Jk(:,i,:) = C_inv_sI_A*diag(inp_mat.BP.K(v_bc,i))*Linv_K; %
    end % B is independent of K
    
    J.Jk = Jk;
else Jk = [];
end

%% +++++++++++++++++ P derivative +++++++++++++++++++++++++
% dA/dp = 0; dG/dp = C*(sI-A)^{-1}*dB/dp %
if inp_mat.cases.P_on
    Jp = zeros(size(C_inv_sI_A,1),length(gamma.P),size(B,2)); % make same size vectors
    for i = length(gamma.P):-1:1 % this could maybe improved
        Jp(:,i,1) = C_inv_sI_A*(inp_mat.BP.P(v_bc,i));
    end
    
    J.Jp = Jp;
else Jp = [];
end
%% Combine
J.all = [Jd,Jv,Jk,Jp];

% reshape in J form
J.Ju = reshape(J.all,size(J.all,1)*size(J.all,2),size(J.all,3)); % reshape 


end

% Try to simplify Jacobian but failed

% A1 = Jcb;
% v = Linv_D(:,1);
% 
% A2 = reshape(Jcb, [], numel(v));      % flatten first two dimensions
% B2 =  A2*v;
% B3  = reshape(B2, size(A1, 1), []);

% AA1 = repmat(Linv_D(:,1),1,length(gamma.D));
% AA  = reshape(AA1,[1,1,length(gamma.D)]);
% size(repmat(Linv_D(:,1),1,length(gamma.D)))
% J11 = bsxfun(@times,Jcb,AA);
