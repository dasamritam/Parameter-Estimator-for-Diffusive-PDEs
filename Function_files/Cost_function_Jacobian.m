function [e,J] = Cost_function_Jacobian(inp_mat,gamma,omega,Y,U)
%This function calculates the weighted cost function using the transfer
%function and its Jacobian
%   Detailed explanation goes here (uitschrijven)


%% Calculate cost-function/Jacobian multiple frequencies where
% e = Y-G*U;
% de/dg = -dG/dg*U 
[W,H] = LS_weighting(Y, inp_mat);


for k = length(omega):-1:1  
%    [G,J_G] = G_and_J_G(inp_mat,gamma,omega(k));
    [G,J_G] = G_and_J_G_nonuniform_grid(inp_mat,gamma,omega(k));
    % Cost function
    r1(:,k) = Y(k,:).'-G*U(k,:).'; % non-regularized part
    r2(:,k) = G*U(k,:).';          % regularized part

    ew = r1;%[W.*r1;H.*r2];
    
    
    %% Jacobian
    JgB = -J_G.Ju*U(k,:).'; % reshape for input
    Jr1(:,k,:) = reshape(JgB,size(J_G.all,1),size(J_G.all,2));  % non-regularized part and reshape back measurements x unknowns    
    Jr2(:,k,:) = reshape(-JgB,size(J_G.all,1),size(J_G.all,2)); % regularized part and reshape back measurements x unknowns    
    
    Jw = Jr1;% [W.*Jr1;H.*Jr2]
        
end

%% Reshape output     V = sum(sum(|e(meas,omega)|^2)) = sum(|e|^2(:))
e = ew(:);
%% Reshape Jacobian V_J = J(meas,omega,unknowns)
J = reshape(Jw,size(Jw,1)*size(Jw,2),size(Jw,3));

end 





%% Compute reshape
% A = J_G.all
% B = reshape(A,size(A,1)*size(A,2),size(A,3))*U(1,:).'
% C = reshape(B,size(A,1),size(A,2));
% 
% C2 = J_G.all(:,:,1)*U(1,1).'+J_G.all(:,:,2)*U(1,2).'+J_G.all(:,:,3)*U(1,3).';