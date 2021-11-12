function [x_old, x_b, x_e, N] = initial_model_grid(x_meas, dxsim, simulation_grid_YES)

  % The main purpose of differentiating two different gridding is to ...
  % enusure that the x_meas are included in the grid. If you uniform ...
  % spacing sensors the simulation_grid_YES == 1 is applicable. ...
  % Otherwise, use the later. 
  
if simulation_grid_YES == 1
    %%%%% Same grid length for data generator and the model.
    dx = dxsim; 
    x_b = x_meas(1); x_e  = x_meas(end);
    x_old = x_b:dx:x_e;
    N = length(x_old);
else
    %%%%% You have the freedom to choose your own grid. 
    N_initial = 301;
    x_b = x_meas(1); x_e  = x_meas(end); 
    % N and not (N-1) because N-1 results in edge sensors to be off ...
    % the grid, whereas in the case of N they are on the grid.
    dx = (x_e-x_b)/(N_initial);                      
    x_old = x_b:dx:x_e;
    N = N_initial;
end