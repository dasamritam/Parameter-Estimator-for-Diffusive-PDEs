function [P,Bpoly] = Polynomial_generation(x,Pin,degree)
%This function calculates polynomial approximation of unknown parameter
%   Input x and the quantity Pin (D,V,K,P)
%   Degree of polynomial to be estimated


for i          = (degree+1):-1:1;
    Bpoly(:,i) = (x(:)).^(i-1);
end

P = Bpoly\Pin(:);


end

