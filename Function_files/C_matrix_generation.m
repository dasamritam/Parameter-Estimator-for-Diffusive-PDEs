function [C,S] = C_matrix_generation(x,x_meas,N)
%Construct matrix with measurement locations for full x (needs to be
%reduced depending on BC.
%   Measurement are approximated to grid instead of making grid based on
%   measurement locations
%   S are the sensor locations
%% Exclude sensors which are being used as boundary conditions
S_length = length(x_meas); % number of sensors

%% Find sensor locations
Sensor = interp1(x,x,x_meas,'nearest'); % select nearest sensor

for i = S_length:-1:1
    indSf(i) = find(x==Sensor(i));
end

S = x(indSf); % grid locations closest to the sensor locations 
if norm(x_meas-S) > 1e-2; warning('error on measurement location compared to interpolated larger than 1%, increase N'); end 

% generate C matrix
C = zeros(S_length,N); j = 1;
for i = indSf 
    C(j,i) = 1;
    j=1+j;
end

C = sparse(C);

%% Reduce size of matrix depending on boundary conditions
% if bcl == 1 && bcr == 1; % Sensor boundaries
%     Cred = C(2:end-1,2:end-1);
% elseif bcl == 0 && bcr == 1; % Neumann left
%     Cred = C(1:end-1,1:end-1);
% elseif bcl == 1 && bcr == 0; % Neumann right
%     Cred = C(2:end,2:end);
% else Cred = C; % Neumann right and left
% end
% 
%    
% if bcl == 1 && bcr == 1;
%     si = S_length-1:-1:2;
% elseif bcl == 1 && bcr == 0;
%     si = S_length:-1:2;
% elseif bcl == 0 && bcr == 1;
%     si = S_length-1:-1:1;
% else si = S_length:-1:1;
% end



end

