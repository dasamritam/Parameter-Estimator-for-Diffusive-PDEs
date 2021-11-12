function [W,L] = LS_weighting(Y,inp_mat)
%%%%%%%% Least Square Weight on Points and on Frequency Lines %%%%%%%
 
% % Update: 23-08-2018
 
%% Spatial weighting
% % For weighting the error around the first two and last two 
% % boundary points.

wx = ones(size(Y,2),1)'; % intiate the spatial weighting matrix 
%wx(1:2) = 1/1e2; % Weighting the first two points
%wx(end-1:end) = 1/1e2; % Weighting the last two points
 
%% Frequency weighting
% % For weighting the error around different frequency lines
% % We keep equal weight for all frequency lines.
 
wf = ones(1,size(Y,1)); % intiate the frequency weighting matrix
 
%% Spatial Weight .* Frequency Weight
W = wx.'*wf; 

%% Regularization weight
L = 1e-10*wx.'*wf; % bullshit

end

