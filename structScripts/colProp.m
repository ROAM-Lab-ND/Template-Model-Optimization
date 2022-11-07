%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% colProp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 December 2021
% Last Updated: 8 December 2021

% This function is used to build the collocation propogation and call to
% the dynamics stored in a CasADi-built function object.

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   dynS - Structure used for dynamics storage
%   problem - Structure used to pass necessary information to NLP solver
%   params - Cell containing variable information for current state

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   problem - Structure used to pass necessary information to NLP solver
%   stateKEnd - Resulting propogation of state variables at current step
%   Q - resulting objective cost of current step
%   R - resulting objective cost of current step

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, problem, stateKEnd, Q, R] =...
    colProp(varS, dynS, problem, params)

    % Import CasADi toolbox
    import casadi.*
    
    % Grab shoot count
    k = params{1};
    
    % Grab state variables
    stateK = params{2};
    
    % Grab non-state variables
    P = params{3};
    
    % Check if in double support or single support
    if contains(varS.params.phase, 'DS')
        
        % Grab input variable
        U = params{4};
        
    end
    
    % Initialize the propogation of the state variables
    stateKEnd = dynS.contEqCoeff(1)*stateK;
    
    % Initialize cost variables
    Q = 0;
    R = 0;
    
    % For each finite element
    for j = 1:varS.params.M
        
        % Create state variable
        [varS, problem] = varCreation(varS, problem, ['state'...
            num2str(k) '_'], j, length(stateK),...
            varS.init.state(k,:)', varS.lb.state, varS.ub.state);
            
    end
    
    % Calculate time step
    DT = 1/varS.params.N;
    
    % For each finite element
    for j = 1:varS.params.M
        
        % Initialize state ROC
        stateKDot = dynS.colEqCoeff(1, j+1)*stateK;
        
        % For each finite element
        for r = 1:varS.params.M
            
            % Calculate state ROC
            stateKDot = stateKDot + dynS.colEqCoeff(r+1, j+1)*...
                varS.var.(['state' num2str(k) '_' num2str(r)]);
            
        end
        
        % Check if in double support or single support
        if contains(varS.params.phase, 'DS')
        
            % Call to CasADi function for state and cost
            [FJ, QJ, RJ] =...
                dynS.fD(varS.var.(['state' num2str(k) '_' num2str(j)]),...
                P, U);
        
        else
            
            % Call to CasADi function for state and cost
            [FJ, QJ, RJ] =...
                dynS.fS(varS.var.(['state' num2str(k) '_' num2str(j)]), P);
            
        end
        
        % Constraint: State ROC and state dynamics from CasADi must align
        [problem] = addConstraint(problem,...
            DT*FJ - stateKDot, zeros(length(stateK),1),...
            zeros(length(stateK),1));
        
        % Calculate state propogation
        stateKEnd = stateKEnd + dynS.contEqCoeff(j+1)*...
            varS.var.(['state' num2str(k) '_' num2str(j)]);
        
        % Calculate running costs
        Q = Q + dynS.quadFuncCoeff(j+1)*QJ*DT;
        R = R + dynS.quadFuncCoeff(j+1)*RJ*DT;
        
    end
    
end