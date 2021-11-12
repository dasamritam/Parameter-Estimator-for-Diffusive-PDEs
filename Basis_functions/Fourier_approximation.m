function [Q,cof ] = Fourier_approximation(x,theta,degree)
% This function compute the fourier series according to the chosen degree
% x: the spatial domain
% theta: parameter profile over the spatial domain x
% degree: the estimation degree we want to use

% Use original basis: function [F, cof]
% Use unweighted orthonormal basis: function [Q, cof] and choose w
% = ones(size(A,1),1)
% Use weighted orthonormal basis: function [Q, cof] and choose w
% as sigmoid function

F = ones(length(x), degree);
scale = (x(end)-x(1))/(2*pi);% scale it to the spatial domain [0.1,0.9]

% Basis
% eg: If degree = 5, basis functions are 1, cos(x), sin(x), cos(2x), sin(2x)
%     If degree = 4, basis functions are 1, cos(x), sin(x), cos(2x)
for i = 2:2:degree
    F(:,i) = cos(1/scale * (i/2) * (x - 0.5));
end

for i = 3:2:degree
    F(:,i) = sin(1/scale * fix(i/2) *(x - 0.5));
end

% theta = theta';
%  cof = (theta\F)';
%  cof = 0.5*ones(7,1)+3;
 
%% Initial 
A = F;
%% Classical Gram-Schmidt 

% for j  = 1:size(A,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) =Q(:,i)'*A(:,j);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end

%% Weight Gram-Schmidt orthonormalization


w = ones(size(F,1),1);
%  w(1:50) = 1e2;
%  w(350:400) = 1e2;

% Weighting function chosen as a combination of two sigmoid functions
% scale2 = (x(end) - x(1))/10; % scale it to the spatial domain [0.1,0.9]
% w = (-dsigmf(1/scale2*(x-0.1),[5 2 5 8]) + 1)*99 + 1; % check the dsigmf function on mathwork website
% w = w';

% for j  = 1:size(F,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) = sum(v.*Q(:,i).*w)/sum(Q(:,i).*Q(:,i).*w);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
% %     R(j,j) = sqrt(sum(w.*v.^2));
%     Q(:,j) = v/R(j,j); 
% %     Q(:,j) = QQ(:,j)/(sqrt(sum(w .* QQ(:,j).^2)));

% end
for j  = 1:size(F,2)
    V(:,j) = A(:,j);
    for i = 1:j-1
        R(i,j) = sum(V(:,i).*A(:,j).*w)/sum(V(:,i).*V(:,i).*w);
        V(:,j) = V(:,j) - R(i,j)*V(:,i);
    end
%     R(j,j) = norm(V(:,j));
    R(j,j) = sqrt(sum(w.*V(:,j).^2));
    Q(:,j) = V(:,j)/R(j,j);
%     Q(:,j) = QQ(:,j)/(sqrt(sum(w .* QQ(:,j).^2)));
end

 cof = Q\(theta)'*0;% coefficients

end

