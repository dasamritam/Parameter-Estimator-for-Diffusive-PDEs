function [gamma, CostBest, CovP] = LM_optimization_Bspline(inp_mat,gamma,omega_all,Y,U,xdom,options)

% Copyright (c) 2017 Matthijs van Berkel, Gerd Vandersteen
% All rights reserved.
% Software can be used freely for non-commercial applications only. 
% If used in publications please reference 

% Compute the MLE error vector and the Jaobian matrix

%if (nargin < 7)
    idxFree=[];
%end
if isempty(idxFree)
    idxFree = 1:length(gamma);
end

[e,J] = Cost_DVKP_Bspline(inp_mat,gamma,omega_all,Y,U,xdom);
CostBest = real(e'*e);
            lambda0 = norm(J);

counter = 0;
lambda = 0;

while ((counter < options.maxNbfOfSteps) & (lambda < 1e10*lambda0));
    % While the cost in decreasing, perform a Newton-Gauss
    % optimization step.
    counter = counter + 1;
    
    % Compute the necessary step.
    if (lambda ~= 0)
        % Levenberg-Marquart
        dP = - [real(J(:,idxFree)); imag(J(:,idxFree)); lambda*eye(length(idxFree))]\[real(e); imag(e); zeros(length(idxFree), 1)];
    else
        % Gauss-Newton P = P0 - (Jr)^(-1)*e(P)
        dP = - [real(J(:,idxFree)); imag(J(:,idxFree))]\[real(e); imag(e)];
    end
    
    Gamma_new = gamma;
    Gamma_new(idxFree) = gamma(idxFree)+dP;
    
    % Evaluate the new optimized point.
    [eNew,JNew] = Cost_DVKP_Bspline(inp_mat,gamma,omega_all,Y,U,xdom);
    
    %    NewCost = real(NewJ(:,end)'*NewJ(:,end));
    CostNew = real(eNew'*eNew);
 
    % If the step is succesful: take the step and decrease lambda
    if (CostNew < CostBest)
        J = JNew;
        e = eNew;
        gamma = Gamma_new;
        CostBest = CostNew;
        lambda = lambda / 3;
    else
        % If the step is not succesful: use lamdba to controle the stepsize
        if (lambda == 0)
            lambda0 = norm(JNew);
            lambda = lambda0;
        else
            lambda = lambda * 10;
        end
    end
end


J_RI = [real(J(:,idxFree)); imag(J(:,idxFree))];
CovP = zeros(length(gamma),length(gamma));
CovP(idxFree,idxFree) = 0.5*inv(J_RI'*J_RI);

fprintf('counter = %d ,lambda = %d ,cond = %d, \n',counter,lambda,cond(J));

