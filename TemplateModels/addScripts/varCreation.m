%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% varCreation %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to create CasADi variables and track them in the
% appropriate indices.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem] = varCreation(vS, problem, name, num, len, first, lb, ub) 

    import casadi.* 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%   
        
    [problem, vS.var.([name num2str(num)]), vS.inds.(inds)(end+1,:)] =...
        addVariable(problem, [name num2str(num)], len, first, lb, ub);    

end