function [e,J] = Cost_function_BsplineDVKP(inp_mat,gamma,omega_all,Y,U,xdom)
%This function calculates the weighted cost function using the transfer
%function and its Jacobian
%   Detailed explanation goes here (uitschrijven)



%% Calculate cost-function/Jacobian multiple frequencies where
% e = Y-G*U;
% de/dg = -dG/dg*U 

for k = length(omega_all):-1:1  
    [G,J_G] = Bspline_G_and_J_G_case(inp_mat,gamma,omega_all(k));
    % Cost function
    enw(:,k) = Y(k,:).'-G*U(k,:).';
    % Jacobian
    Jg = J_G.all;
    JgB = -reshape(Jg,size(Jg,1)*size(Jg,2),size(Jg,3))*U(k,:).'; % reshape for input
    Jnw(:,k,:) = reshape(JgB,size(Jg,1),size(Jg,2)); % reshape back measurements x unknowns    
end

%% Weighting
% Spatial weighting
deltax = xdom(2:end)-xdom(1:end-1); wx = 1./deltax(1:size(Jnw,1)); % Riemann approximation
%W = ones(size(W)); % One could also program simpson rule by weights

% Frequency weighting
wf = ones(size(omega_all)); % sigma e
W = wx.'*wf; 

%% Weight output
% Weight cost function
e = W.*enw;

% Weight Jacobian
Wfull = repmat(W,[1,1,size(Jnw,3)]); % check weights later on
J = Wfull.*Jnw;

end 




%% Compute reshape
% A = J_G.all
% B = reshape(A,size(A,1)*size(A,2),size(A,3))*U(1,:).'
% C = reshape(B,size(A,1),size(A,2));
% 
% C2 = J_G.all(:,:,1)*U(1,1).'+J_G.all(:,:,2)*U(1,2).'+J_G.all(:,:,3)*U(1,3).';