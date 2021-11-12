function [Jnum] = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,Bnumber)


P = gamma.K;

lri = size(inp_mat.C,1);
Jnum = zeros(lri,length(P));


idxFree = 1:length(gamma.K);

for k=length(idxFree):-1:1
    
    % profielen
    theta = P;
    thetad = P;
    
    % left side
    thetad(idxFree(k)) = thetad(idxFree(k))*(1-dtheta);
    % right side
    theta(idxFree(k))  = theta(idxFree(k))*(1+dtheta);

      
    % thetad (left_side)
    gamma.K = thetad;
    [Gd,~] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);
    e1 = Gd(:,Bnumber);
    eo1 =[real(e1);imag(e1)];
    

    
    % theta (right_side)
    gamma.K = theta;
    [G,~] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);
    e2 = G(:,Bnumber);

    eo2=[real(e2);imag(e2)];

    % Calculate Jacobian
    J(:,idxFree(k))=(eo2-eo1)/(2*dtheta*theta(idxFree(k)));
    Jnum(:,idxFree(k)) = J(1:lri,idxFree(k))+J(lri+1:end,idxFree(k))*1i;
end
end