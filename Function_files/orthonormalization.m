function [Q] = orthonormalization(BP, theta)
%% Orthonormalization
% Initialize
A = BP;

% % A=inp_mat.BP.D_nonorth;
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

