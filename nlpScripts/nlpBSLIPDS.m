%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% nlpBSLIPDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 17 November 2021

% This function is used to propagate the BSLIP model through the Double
% Support (DS) phase.

% INPUTS:
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

function [varS, problem, stateEnd] =...
    nlpBSLIPDS(varS, dynS, problem, i, stateK)
    
    import casadi.*

    % Iterate through phase for desired number of shooting elements
    for k = varS.params.shootCount:...
            (varS.params.shootCount + varS.params.N - 1)
        
        % Store current state information for potential future use
        stateEnd = stateK;
        
        % Set phase indicator to double support
        varS.params.phase = 'DS';
        
        % Create input control variable u (Currently forced to 0)
        [varS, problem] = varCreation(varS, problem, 'u',...
            varS.params.inputCountD, 1,...
            varS.init.u(k-((i-1)*varS.params.N)), varS.lb.u, varS.ub.u);
        
        % Update time tracking variable
        tDat = varS.params.timeTrack +...
            (((k) - 2*(i-1)*varS.params.N)/varS.params.N)*...
            varS.var.(['tfD' num2str(i)]);

        % Store non-state/non-control parameters in P vector
        P = [varS.var.(['state0D' num2str(i)]);...
            varS.var.(['thetaD' num2str(i)]);...
            varS.var.(['kDLagDot' num2str(varS.params.inputCountD)]);...
            varS.var.(['kDLeadDot' num2str(varS.params.inputCountD)]);...
            varS.var.(['footDLag' num2str(i)]);...
            varS.var.(['tfD' num2str(i)]); tDat; varS.var.len0];
        
        % Propogate state based on multishooting or collocation
        if strcmp(varS.params.method, 'Multishooting')
            
            % Pass state, parameter, and control variables to integrator
            Fk = dynS.FD('state0', stateK, 'p', P, 'u',...
                varS.var.(['u' num2str(varS.params.inputCountD)]));
            
            % Store propogated state for future use
            stateEnd = Fk.stated;

            % Update running cost
            problem.cost = problem.cost + Fk.qd;
            problem.resVert = problem.resVert + Fk.rd;
            
        elseif strcmp(varS.params.method, 'Collocation')
            
            % Store parameters in vector for collocation propogation
            params = {k, stateK, P,...
                varS.var.(['u' num2str(varS.params.inputCountD)])};
            
            % Pass params to integrator
            [varS, problem, stateEnd, Q, R] =...
                colProp(varS, dynS, problem, params);
        
            % Update running cost
            problem.cost = problem.cost + Q;
            problem.resVert = problem.resVert + R;
            
        end
        
        % Iterate double support counter
        varS.params.inputCountD = varS.params.inputCountD + 1;
        
        % Check if current step is the final iteration for this phase
        if (k ~= (varS.params.shootCount + varS.params.N - 1))
        
            % Create state variables for shoot k+1
            [varS, problem] = varCreation(varS, problem, 'state', k+1,...
                6, varS.init.state(k+1,:)', varS.lb.state, varS.ub.state);
            
            % Store state in new variable for easier constraint coding
            stateK = varS.var.(['state' num2str(k+1)]);

            % Constraint: State at (k+1) must start where state at k ends
            [problem] = addConstraint(problem, stateEnd - stateK,...
                zeros(6,1), zeros(6,1));
            
            % Create control variables for shoot k+1
            [varS, problem] = varCreation(varS, problem, 'kDLagDot',...
                varS.params.inputCountD, 1,...
                varS.init.kLagDot(k+1), varS.lb.kDot, varS.ub.kDot);
            
            [varS, problem] = varCreation(varS, problem, 'kDLeadDot',...
                varS.params.inputCountD, 1,...
                varS.init.kLeadDot(k+1), varS.lb.kDot, varS.ub.kDot);
            
            % Store hip position
            xHip = stateK(1);
            yHip = stateK(2);
            
            % Calculate lag leg length
            lenLag = sqrt((xHip -...
                varS.var.(['footDLag' num2str(i)]))^2 + yHip^2);
            
            % Constraint: Force along lag leg must never be less than 0
            [problem] = addConstraint(problem,...
                1000*(stateK(5)*(varS.var.len0 - lenLag) +...
                varS.var.(['u' num2str(varS.params.inputCountD-1)])*...
                yHip/(lenLag)), 0, inf);
            
            % Check if model is using constant or varying leg stiffness
            if ~strcmp(varS.params.springType, 'Constant')
                
                % Constraints: Leg stiffness ROC must be within bounds
                [problem] = addConstraint(problem,...
                    varS.var.(['kDLagDot'...
                    num2str(varS.params.inputCountD)]) -...
                    varS.var.(['kDLagDot'...
                    num2str(varS.params.inputCountD-1)]),...
                    -varS.params.stiffVaryMax, varS.params.stiffVaryMax);
                
                [problem] = addConstraint(problem,...
                    varS.var.(['kDLeadDot'...
                    num2str(varS.params.inputCountD)]) -...
                    varS.var.(['kDLeadDot'...
                    num2str(varS.params.inputCountD-1)]),...
                    -varS.params.stiffVaryMax, varS.params.stiffVaryMax);

            end
            
        end
        
    end
    
    % Update overall instance count to k+1
    varS.params.shootCount = k + 1;

end