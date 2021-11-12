function [e_ideal] = Cost_function_Jacobian_ideal(inp_mat,gamma,omega,Y,U)
%This function calculates the weighted cost function using the transfer
%function and its Jacobian
%   Detailed explanation goes here (uitschrijven)


%% Calculate cost-function/Jacobian multiple frequencies where
% e = Y-G*U;
% de/dg = -dG/dg*U 
[W,H] = LS_weighting(Y, inp_mat);


for k = length(omega):-1:1  
   % [G,~] = G_and_J_G(inp_mat,gamma,omega(k));
    [G,~] = G_and_J_G_nonuniform_grid(inp_mat,gamma,omega(k));
    % Cost function
    r1(:,k) = Y(k,:).'-G*U(k,:).'; % non-regularized part
    r2(:,k) = G*U(k,:).';          % regularized part

    ew = [W.*r1;H.*r2];
    
 
        
end

%% Reshape output     V = sum(sum(|e(meas,omega)|^2)) = sum(|e|^2(:))
e_ideal = ew(:);

end 
