function [Q,gamma_i,BPc] = Bspline_generation_orthonormal(c,theta,x,degree)
% function [BP,gamma,BPc] = Bspline_generation(c,y,x,degree)
% generates the b-spline matrix and calculates the initial
% fit to the initial guess of the transport
% coefficient  using least-squares fit: gamma = BP\theta
% 
% Input
%   c       control point locations uniformly spaced (C)
%   theta   value of transport coefficient at c, ordinates of control
%           points
%   x       discretized locations (N)
%   degree  B-spline degree (cubic is degree = 4)
%
% Output
%   BP      matrix with B-spline values on x-grid
%   gamma_i B-spline coefficients, i.e., y = BP*gamma
%   BPc     matrix with B-spline values only on c-grid

% Use original basis: function [BP, gamma_i, BPc]
% Use unweighted orthonormal basis: function [Q, gamma_i, BPc] and choose w
% = ones(size(A,1),1)
% Use weighted orthonormal basis: function [Q, gamma_i, BPc] and choose w
% as sigmoid function


nr_extra_knots = degree+1; % default
delta_right = (c(end)-c(end-1))*1;
knots_extra_right = linspace(c(end)+delta_right,c(end)+...
                           nr_extra_knots*delta_right,nr_extra_knots);
                       
knots=[c knots_extra_right]; % knots extra right don't understand ask Peter
    
deltac = diff(knots); % Calculate the difference between each knot
knots  = knots-(degree+1)/2.*[deltac deltac(end)];
%% Generate sparse grid B-spline to calculate initial gamma
BPc = 0;
j=0;    
for s_point = c;
    j=j+1;
    for i=length(theta):-1:1,
        BPc(j,i) = BBspline(degree,i,s_point,knots);    
    end       
end,
    
gamma_i = BPc\theta(:)*0;

%% Generate B-splines on finite difference grid for evaluation

BP = 0;
j=0;    
for x_point = x
    j=j+1;
    for i = length(theta):-1:1,
        BP(j,i) = BBspline(degree,i,x_point,knots);
    end       
end
    
%% Orthonormalization
% Initialize
A = BP;


% Weighting for WGS orthonormalization
w = ones(size(A,1),1);


% Weight Gram-Schmidt orthonormalization
for j  = 1:size(A,2)
    V(:,j) = A(:,j);
    for i = 1:j-1
        R(i,j) = sum(V(:,i).*A(:,j).*w)/sum(V(:,i).*V(:,i).*w);
        V(:,j) = V(:,j) - R(i,j)*V(:,i);
    end
    R(j,j) = sqrt(sum(w.*V(:,j).^2));
    Q(:,j) = V(:,j)/R(j,j);
end
