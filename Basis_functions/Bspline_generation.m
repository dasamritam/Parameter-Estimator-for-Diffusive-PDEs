function [BP,gamma_i,BPc] = Bspline_generation(c,theta,x,degree)
% function [BP,gamma,BPc] = Bspline_generation(c,y,x,degree)
% generates the b-spline matrix and calculates the initial
% fit to the initial guess of the transport
% coefficient  using least-squares fit: gamma = BP\theta
% 
% Input
%   c       control point locations uniformly spaced (C)
%   theta   value of transport coefficient at c
%   x       discretized locations (N)
%   degree  B-spline degree (cubic is degree = 4)
% Output
%   BP      matrix with B-spline values on x-grid
%   gamma_i B-spline coefficients, i.e., y = BP*gamma
%   BPc     matrix with B-spline values only on c-grid

nr_extra_knots = degree+1; % default
delta_right = (c(end)-c(end-1))*1;
knots_extra_right = linspace(c(end)+delta_right,c(end)+...
                           nr_extra_knots*delta_right,nr_extra_knots);
                       
knots=[c knots_extra_right]; % knots extra right don't understand ask Peter
    
deltac = diff(knots);
knots  = knots-(degree+1)/2.*[deltac deltac(end)];

%% Generate sparse grid B-spline to calculate initial gamma
BPc = 0;
j=0;    
for s_point = c
    j=j+1;
    for i=length(theta):-1:1
        BPc(j,i) = BBspline(degree,i,s_point,knots);    
    end       
end
    
gamma_i = BPc\theta(:);

%% Generate B-splines on finite difference grid for evaluation

BP = 0;
j=0;    
for x_point = x
    j=j+1;
    for i = length(theta):-1:1
        BP(j,i) = BBspline(degree,i,x_point,knots);
    end       
end
%BP = struct('BP',BP);
%% For testing only
% %%
% clear all; clc;
% x = 0:1e-3:2.2; 
% s = linspace(0,2.2,10);
% y = -s.^3+s+5;
% degree = 4;
% 
% y1=BP*gamma; 
% 
% figure
% plot(x,BP,'r',s,BPc,'ko')
% 
% figure
% plot(x,y1,'k',s,BPc*gamma,'o-')
