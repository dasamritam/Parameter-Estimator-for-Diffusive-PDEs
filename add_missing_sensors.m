function [x, dx, N, miss] = add_missing_sensors(x_old, x_meas)

% Step 1:: Find the missing points in the initial grid.

      missing_index=1;
      miss = [];

      for i=1:length(x_meas)
          if isempty(find(round(x_old,8)==round(x_meas(i),8)))== 1 
             miss(missing_index) = i;
             missing_index=missing_index+1;
          end
      end

% Step 2:: Add the missing ponits in order.
       
      if isempty(miss)==0
         x_new = [x_old x_meas(miss)];
         x = round(sort(x_new),8);
      else
         x =round(x_old,8);
      end

% Calculating Non-uniform grid length
        
      index=1;
      dx_new = [];
    for j = 2:length(x)-1
        dx_new.dx_forward(index) = round(x(j+1)-x(j),5);
        dx_new.dx_back(index) = round(x(j)-x(j-1),5);
        index = index+1;
    end

N = length(x);
dx = dx_new;