function [e,J] = Cost_function_BsplineDVKP(inp_mat,gamma,omega_all,Y,U)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
keyboard
for j = size(Y,2):-1:1
    for k = length(omega_all):-1:1
        [G,J_G] = Bspline_G_and_J_G_case(inp_mat,gamma,omega_all(k));
        e(j,k) = Y(k,j)-G(j,:)*U(k,:).';
        %J = J_G.Jd*U(k,:).';
        JGG(:,k) = reshape(J_G.Jd,18*11,3)*U(k,:).'
    end
end
JGGG = reshape(JGG,18,11,4);
JGGGG = reshape(JGGG,72,11);

J_G.Jd*U.'

JGG*U.'

%     for k = length(omega_all):-1:1
%         [G,J_G] = Bspline_G_and_J_G_case(inp_mat,gamma,omega_all(k));
%         e( = Y(:)-G*U(omega,:);
%         J = J_G*U;
%     end

%  e = Y-G*U;
% dedgamma = -dGdgamma*U 

% Build gamma

% Build J

end

