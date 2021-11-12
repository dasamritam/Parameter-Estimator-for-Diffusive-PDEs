function [cases] = Cases_variables_estimated(FreeVariables)
% Sets which variables have to be estimated
%   Various if statements throughout the code dependend on the choice.
%   These are set in this function

%FreeVariables = 'alll';

str_error = sprintf(...
'FreeVariables = %s is not a proper definition of what has to be estimated:\nplease select any combination of DVKP e.g. DV to estimate D and V',FreeVariables);

switch FreeVariables
    % All variables        
    case {'all','DVKP','All','ALL'} % all variables
        cases.D_on = 1; cases.V_on = 1; cases.K_on = 1; cases.P_on = 1;
        disp('all are estimated')
    
    % Three variables
    case {'DVK'}
         cases.D_on = 1; cases.V_on = 1; cases.K_on = 1; cases.P_on = 0;
    case {'DVP'}
         cases.D_on = 1; cases.V_on = 1; cases.K_on = 0; cases.P_on = 1;
    case {'DKP'}
         cases.D_on = 1; cases.V_on = 0; cases.K_on = 1; cases.P_on = 1;
    case {'VKP','VPK','PKV','KPV','KVP','PVK'}
         cases.D_on = 0; cases.V_on = 1; cases.K_on = 1; cases.P_on = 1;

   
    % Two variables           
    case {'DV','VD'}
         cases.D_on = 1; cases.V_on = 1; cases.K_on = 0; cases.P_on = 0;
    case {'DK','KD'}
         cases.D_on = 1; cases.V_on = 0; cases.K_on = 1; cases.P_on = 0;
    case {'DP','PD'}
         cases.D_on = 1; cases.V_on = 0; cases.K_on = 0; cases.P_on = 1;
    case {'VK','KV'}
         cases.D_on = 1; cases.V_on = 0; cases.K_on = 1; cases.P_on = 0;
    case {'VP','PV'}
         cases.D_on = 0; cases.V_on = 1; cases.K_on = 0; cases.P_on = 1;
    case {'KP','PK'}
         cases.D_on = 0; cases.V_on = 0; cases.K_on = 1; cases.P_on = 1;

    % One variables           
    case {'D'}
         cases.D_on = 1; cases.V_on = 0; cases.K_on = 0; cases.P_on = 0;
    case {'V'}
         cases.D_on = 0; cases.V_on = 1; cases.K_on = 0; cases.P_on = 0;
    case {'K'}
         cases.D_on = 0; cases.V_on = 0; cases.K_on = 1; cases.P_on = 0;
    case {'P'}
         cases.D_on = 0; cases.V_on = 0; cases.K_on = 0; cases.P_on = 1;
    otherwise
        error(str_error)
end



