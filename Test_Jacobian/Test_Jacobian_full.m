function [Jnum] = Test_Jacobian_full(inp_mat,gamma,omega_all,Y,U,xdom,dtheta)

% Determine number of unknowns and position in the total vector
gl_d = length(gamma.D); gl_v = length(gamma.V);
gl_k = length(gamma.K); gl_p = length(gamma.P);

% Total vector theta positions of coefficients
ind_d = 1:gl_d; ind_v = gl_d+(1:gl_v);
ind_k = gl_d+gl_v+(1:gl_k); ind_p = gl_d+gl_v+gl_k+(1:gl_p);

% Make full vector
theta = [gamma.D;gamma.V;gamma.K;gamma.P];


for i=length(theta):-1:1
    % left side
    theta_l = theta;
    theta_l(i) = theta(i)*(1-dtheta);
    
    gamma_l.D = theta_l(ind_d); gamma_l.V = theta_l(ind_v); 
    gamma_l.K = theta_l(ind_k); gamma_l.P = theta_l(ind_p);
    
    [eo1,~] = Cost_DVKP_Bspline(inp_mat,gamma_l,omega_all,Y,U,xdom);
    
    % right side
    theta_r = theta;
    theta_r(i) = theta(i)*(1+dtheta);
    
    gamma_r.D = theta_r(ind_d); gamma_r.V = theta_r(ind_v); 
    gamma_r.K = theta_r(ind_k); gamma_r.P = theta_r(ind_p);
    
    [eo2,~] = Cost_DVKP_Bspline(inp_mat,gamma_r,omega_all,Y,U,xdom);
    
    % Calculate Jacobian
    Jnum(:,i)=(eo2-eo1)/(2*dtheta*theta(i));
    %Jnum(:,idxFree(k)) = J(1:lri,idxFree(k))+J(lri+1:end,idxFree(k))*1i;
end   


end


















% 
% 
% 
% P = [gamma.D;gamma.V;gamma.K;gamma.P];%gamma;
% 
% lri = size(inp_mat.C,1);
% Jnum = zeros(lri,length(P));
% 
% 
% idxFree = 1:length(P);
% 
% for k=length(idxFree):-1:1
%     
%     % profielen
%     theta = P;
%     thetad = P;
%     
%     % left side
%     thetad(idxFree(k)) = thetad(idxFree(k))*(1-dtheta);
%     % right side
%     theta(idxFree(k))  = theta(idxFree(k))*(1+dtheta);
% 
%       
%     % thetad (left_side)
%     gamma = thetad;
%     [e1,~] = Cost_function_BsplineDVKP(inp_mat,gamma,omega_all,Y,U);
%     eo1 =[real(e1);imag(e1)];
%     
%   
%     % theta (right_side)
%     gamma = theta;
%     [e2,~] = Cost_function_BsplineDVKP(inp_mat,gamma,omega_all,Y,U);
%     eo2=[real(e2);imag(e2)];
% 
%     % Calculate Jacobian
% %     J(:,idxFree(k))=(eo2-eo1)/(2*dtheta*theta(idxFree(k)));
% %     Jnum(:,idxFree(k)) = J(1:lri,idxFree(k))+J(lri+1:end,idxFree(k))*1i;
% end
% end