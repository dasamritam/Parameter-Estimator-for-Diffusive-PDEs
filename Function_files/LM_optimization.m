function [gamma, CostBest, CovP] = LM_optimization(inp_mat,gamma,omega,Y,U,options)

% Copyright (c) 2017 Matthijs van Berkel, Gerd Vandersteen
% All rights reserved.
% Software can be used freely for non-commercial applications only. 
% If used in publications please reference 

% Compute the MLE error vector and the Jaobian matrix

%% Predefine variables
% Determine number of unknowns and position in the total vector
gl_d = length(gamma.D); gl_v = length(gamma.V);
gl_k = length(gamma.K); gl_p = length(gamma.P);

% Total vector theta positions of coefficients theta = [gamma.D;gamma.V;gamma.K;gamma.P];
ind_d = 1:gl_d; ind_v = gl_d+(1:gl_v);
ind_k = gl_d+gl_v+(1:gl_k); ind_p = gl_d+gl_v+gl_k+(1:gl_p);

% Variable to manually control which parameters which are going to change
lfree = gl_d+gl_v+gl_k+gl_p;
idxFree = 1:(gl_d+gl_v+gl_k+gl_p); %default length theta

%% Initial estimate cost

[e,J] = Cost_function_Jacobian(inp_mat,gamma,omega,Y,U);
CostBest = real(e'*e);
            lambda0 = norm(J);

counter = 0;
lambda = 0;
alpha = 0;

while ((counter < options.maxNbfOfSteps) & (lambda < 1e10*lambda0));
    % While the cost in decreasing, perform a Newton-Gauss
    % optimization step.
    counter = counter + 1;

    % Compute the necessary step.

    if (lambda ~= 0)
        % Levenberg-Marquart
        dtheta = - [real(J(:,idxFree)); imag(J(:,idxFree)); lambda*eye(length(idxFree))]\[real(e); imag(e); zeros(length(idxFree), 1)];

%         JJ = [real(J(:,idxFree)); imag(J(:,idxFree))];
%         dtheta2 = -(JJ'*JJ/lambda+eye(lfree)+alpha*eye(lfree))\(JJ'*[real(e); imag(e)]);
%         norm(dtheta2-dtheta)
    else
        % Gauss-Newton P = P0 - (Jr)^(-1)*e(P)
        dtheta = - [real(J(:,idxFree)); imag(J(:,idxFree))]\[real(e); imag(e)];
%         JJ = [real(J(:,idxFree)); imag(J(:,idxFree))];
%         dtheta2 = -(JJ'*JJ+alpha*eye(lfree))\JJ'*[real(e); imag(e)];
%         norm(dtheta2-dtheta)
    end
    
    % Define step
    theta = ([gamma.D;gamma.V;gamma.K;gamma.P]); 
    theta_new = theta; % initialize for idxFree case
    theta_new(idxFree) = theta(idxFree)+dtheta; % calculate step-size
    
    % Transform theta_new to gamma_new
    gamma_new.D = theta_new(ind_d); gamma_new.V = theta_new(ind_v); 
    gamma_new.K = theta_new(ind_k); gamma_new.P = theta_new(ind_p);
    
    % Evaluate the new optimized point.
    [eNew,JNew] = Cost_function_Jacobian(inp_mat,gamma_new,omega,Y,U);
    
    % NewCost = real(NewJ(:,end)'*NewJ(:,end));
    CostNew = real(eNew'*eNew);
 
    % If the step is succesful: take the step and decrease lambda
    if (CostNew < CostBest)
        J = JNew;
        e = eNew;
        gamma = gamma_new;
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

CovP =0;
% J_RI = [real(J(:,idxFree)); imag(J(:,idxFree))];
% CovP = zeros(length(gamma),length(gamma));
% CovP(idxFree,idxFree) = 0.5*inv(J_RI'*J_RI);

fprintf('counter = %d ,lambda = %d ,cond = %d, \n',counter,lambda,cond(J));

