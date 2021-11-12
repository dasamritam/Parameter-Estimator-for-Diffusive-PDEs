function gamma = theta2gamma(theta,gamma)
%This function transforms gamma to theta
%   Detailed explanation goes here
gamma.D = theta(1:length(gamma.D));
gamma.V = theta(length(gamma.D)+1:length(gamma.D)+length(gamma.V));
gamma.K = theta(length(gamma.V)+1:length(gamma.V)+length(gamma.V));
theta = [gamma.D,gamma.V,gamma.K,gamma.P];

end

