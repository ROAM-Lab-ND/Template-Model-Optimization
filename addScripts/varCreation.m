%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% varCreation %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to create CasADi variables and track them in the
% appropriate indices.

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   problem - Structure used to pass necessary information to NLP solver
%   name - Variable name to add to problem struct
%   num - Variable index number
%   len - Size of variable to add to problem struct
%   first - Initial guess of variable value to seed NLP solver
%   lb - Lower bound value permissible for variable
%   ub - Upper bound value permissible for variable

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   problem - Structure used to pass necessary information to NLP solver

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, problem] =...
    varCreation(varS, problem, name, num, len, first, lb, ub) 

    import casadi.* 
        
    % Check variable name for storing variable indices in correct place
    if contains(name, 'theta')

        inds = 'ang';

    elseif contains(name, 'kSLagDot')||contains(name, 'kDLagDot')

        inds = 'kLagDot';
        
    elseif contains(name, 'kSLeadDot')||contains(name, 'kDLeadDot')
        
        inds = 'kLeadDot';

    elseif contains(name, 'foot')

        inds = 'foot';
        
    elseif contains(name, 'z')
        
        inds = 'state';
        
    elseif contains(name, 'state') && ~contains(name, '_')
        
        inds = 'state';
        
    elseif contains(name, 'state') && contains(name, '_')
        
        inds = 'stateFull';
        
    elseif contains(name, 'u')
        
        inds = 'u';
        
    elseif contains(name, 'tf')
        
        inds = 'time';
        
    elseif contains(name, 'tDat')
        
        inds = 'tHum';
        
    elseif contains(name, 'rVPP')
        
        inds = 'rVPP';

    else

        inds = 'gen';

    end
       
    % Call into addVariable to generate and store variable properly    
    [problem, varS.var.([name num2str(num)]), varS.inds.(inds)(end+1,:)] =...
        addVariable(problem, [name num2str(num)], len, first, lb, ub);    

end