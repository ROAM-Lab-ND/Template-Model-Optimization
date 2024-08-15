%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% colProp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 December 2021
% Last Updated: 8 December 2021

% This function is used to build the collocation propogation and call to
% the dynamics stored in a CasADi-built function.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem, stateKEnd, Q, R] = colProp(vS, dynStruct, problem, params)

    import casadi.*
    
    k = params{1};
    stateK = params{2};
    P = params{3};
    
    if contains(vS.params.phase, 'DS')
        
        U = params{4};
        
    end
    
    % Double Support Collocation
    stateKEnd = dynStruct.contEqCoeff(1)*stateK;
    Q = 0;
    R = 0;
    
    for j = 1:vS.params.M
        
            [vS, problem] = varCreation(vS, problem, ['state'...
                    num2str(k) '_'], j, length(stateK),...
                    vS.init.state(k,:)', vS.lb.state, vS.ub.state);
            
    end
    
    DT = 1/vS.params.N;
    
    for j = 1:vS.params.M
        
        % State Derivative at Collocation Point
        stateKDot = dynStruct.colEqCoeff(1, j+1)*stateK;
        
        for r = 1:vS.params.M
            
            stateKDot = stateKDot + dynStruct.colEqCoeff(r+1, j+1)*...
                vS.var.(['state' num2str(k) '_' num2str(r)]);
            
        end
        
        if contains(vS.params.phase, 'DS')
        
            [FJ, QJ, RJ] = dynStruct.fD(vS.var.(['state' num2str(k) '_' num2str(j)]), P, U);
        
        else
            
            [FJ, QJ, RJ] = dynStruct.fS(vS.var.(['state' num2str(k) '_' num2str(j)]), P);
            
        end
        
        [problem] = addConstraint(problem,...
            DT*FJ - stateKDot, zeros(length(stateK),1),...
            zeros(length(stateK),1));
        
        stateKEnd = stateKEnd + dynStruct.contEqCoeff(j+1)*...
            vS.var.(['state' num2str(k) '_' num2str(j)]);
        
        Q = Q + dynStruct.quadFuncCoeff(j+1)*QJ*DT;
        R = R + dynStruct.quadFuncCoeff(j+1)*RJ*DT;
        
    end
    
end