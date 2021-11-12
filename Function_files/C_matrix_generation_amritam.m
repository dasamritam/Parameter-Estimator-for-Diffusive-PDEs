function [C] = C_matrix_generation_amritam(x,xbc_meas,N)
x_indomain = x(2:N-1);
C = sparse(zeros(length(xbc_meas),length(x_indomain)));
[~,index_measured,~] = intersect(x_indomain, xbc_meas,'sorted');
for op=1:length(index_measured)
    C(op,index_measured(op)) = 1;
    op = op+1;
end
end

