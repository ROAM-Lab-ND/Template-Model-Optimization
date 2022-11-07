%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% nlpBSLIPSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to propagate the BSLIP model through the Single
% Support (SS) phase.

%   varS - Structure used for variable and post-optimization storage
%   dynS - Structure used for dynamic storage
%   problem - Structure used to pass necessary information to NLP solver
%   i - Iterator tracking number of steps taken by model
%   stateK - State vector at current time instance

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   problem - Structure used to pass necessary information to NLP solver
%   stateEnd - State vector after propogating final time instance of phase

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, problem, stateEnd] = ...
    nlpBSLIPSS(varS, problem, dynS, i, stateK)
    
    import casadi.*
    
    % Iterate through phase for desired number of shooting elements
    for k = varS.params.shootCount:...
            (varS.params.shootCount + varS.params.N - 1)
        
        % Store current state information for potential future use
        stateEnd = stateK;
        
        % Set phase indicator to single support
        varS.params.phase = 'SS';
        
        % Update time tracking variable
        tDat = varS.params.timeTrack +...
            ((k - (2*i - 1)*varS.params.N)/varS.params.N)*...
            varS.var.(['tfS' num2str(i)]);
        
        % Store non-state/non-control parameters in P vector
        P = [varS.var.(['kSLagDot' num2str(varS.params.inputCountS)]);...
            varS.var.(['kSLeadDot' num2str(varS.params.inputCountS)]);...
            varS.var.(['footSStance' num2str(i)]);...
            varS.var.(['tfS' num2str(i)]); tDat; varS.var.len0];
        
        % Propogate state based on multishooting or collocation
        if strcmp(varS.params.method, 'Multishooting')
            
            % Pass state and parameter variables to integrator
            Fk = dynS.FS('state0', stateK, 'p', P);        

            % Store propogated state for future use
            stateEnd = Fk.states;

            % Update running cost 
            problem.cost = problem.cost + Fk.qs;
            problem.resVert = problem.resVert + Fk.rs;
            
        elseif strcmp(varS.params.method, 'Collocation')
            
            % Store parameters in vector for collocation propogation
            params = {k , stateK , P, []};
            
            % Pass params to integrator
            [varS, problem, stateEnd, Q, R] =...
                colProp(varS, dynS, problem, params);

            % Update running cost
            problem.cost = problem.cost + Q;
            problem.resVert = problem.resVert + R;
            
        end
        
        % Iterate single support counter
        varS.params.inputCountS = varS.params.inputCountS + 1;
        
        % Check if current step is the final iteration for this phase
        if (k ~= (varS.params.shootCount + varS.params.N - 1))
            
            % Create state variables for shoot k+1
            [varS, problem] = varCreation(varS, problem, 'state', k+1, 6,...
                varS.init.state(k+1,:)', varS.lb.state, varS.ub.state);
            stateK = varS.var.(['state' num2str(k+1)]);

            % Constraint: State at (k+1) must start where state at k ends
            [problem] = addConstraint(problem, stateEnd - stateK,...
                zeros(6,1), zeros(6,1));
            
            % Create control variables for shoot k+1
            [varS, problem] = varCreation(varS, problem, 'kSLagDot',...
                varS.params.inputCountS, 1,...
                varS.init.kLagDot(k+1), varS.lb.kDot, varS.ub.kDot);

            [varS, problem] = varCreation(varS, problem, 'kSLeadDot',...
                varS.params.inputCountS, 1,...
                varS.init.kLeadDot(k+1), varS.lb.kDot, varS.ub.kDot);
            
            % Check if model is using constant or varying leg stiffness
            if ~strcmp(varS.params.springType, 'Constant')
                
                % Check if within two shooting elements of end of phase
                if (k ~= (varS.params.shootCount + varS.params.N - 2))

                    % Constraint: Lag leg stiffness ROC within bounds
                    [problem] = addConstraint(problem,...
                        varS.var.(['kSLagDot'...
                        num2str(varS.params.inputCountS)]) -...
                        varS.var.(['kSLagDot'...
                        num2str(varS.params.inputCountS-1)]),...
                        -varS.params.stiffVaryMax,...
                        varS.params.stiffVaryMax);

                end
                
                % Constraint: Lead leg stiffness ROC within bounds
                [problem] = addConstraint(problem,...
                    varS.var.(['kSLeadDot'...
                    num2str(varS.params.inputCountS)]) -...
                    varS.var.(['kSLeadDot'...
                    num2str(varS.params.inputCountS-1)]),...
                    -varS.params.stiffVaryMax, varS.params.stiffVaryMax);

            end
            
        end

    end
    
    % Update overall instance count to k+1
    varS.params.shootCount = k + 1;
    
end