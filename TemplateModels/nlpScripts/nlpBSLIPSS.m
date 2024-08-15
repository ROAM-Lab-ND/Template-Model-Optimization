%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% nlpBSLIPSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to propagate the BSLIP model through the Single
% Support (SS) phase.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem, stateEnd] = nlpBSLIPSS(vS, problem, dynStruct, i, stateK)
    
    import casadi.*
    
    % Single Support Phase
    for k = vS.params.shootCount:(vS.params.shootCount + vS.params.N - 1)
        
        statePrev = stateK;
        vS.params.phase = 'SS';
        
%         % Iterate input counter and create variable for input control at
%         % time step k
%         vS.params.inputCountS = vS.params.inputCountS + 1;
%         [vS, problem] = varCreation(vS, problem, 'kSLagDot',...
%             vS.params.inputCountS, 1,...
%             vS.init.kLagDot(k+1), vS.lb.kDot, vS.ub.kDot);
%         [vS, problem] = varCreation(vS, problem, 'kSLeadDot',...
%             vS.params.inputCountS, 1,...
%             vS.init.kLeadDot(k+1), vS.lb.kDot, vS.ub.kDot);
        
        tDat = vS.params.timeTrack +...
            ((k - (2*i - 1)*vS.params.N)/vS.params.N)*...
            vS.var.(['tfS' num2str(i)]);
        
        % Store non-state/non-control parameters in P vector for passing to
        % RK4 integrator
        P = [vS.var.(['kSLagDot' num2str(vS.params.inputCountS)]);...
            vS.var.(['kSLeadDot' num2str(vS.params.inputCountS)]);...
            vS.var.(['footSStance' num2str(i)]);...
            vS.var.(['tfS' num2str(i)]); tDat; vS.var.len0];
        
        if strcmp(vS.params.method, 'Multishooting')
            
            % Pass time, state, parameter, and control variables to integrator
            Fk = dynStruct.FS('state0', stateK, 'p', P);        

            % Store result of state after propogating shoot through integrator
            stateEnd = Fk.states;

            % Add resulting objective cost from shoot k to total cost 
            problem.cost = problem.cost + Fk.qs;
            problem.resVert = problem.resVert + Fk.rs;
            
        elseif strcmp(vS.params.method, 'Collocation')

            params = {k , stateK , P, []};
        
            [vS, problem, stateEnd, Q, R] = colProp(vS, dynStruct, problem, params);

            % Add resulting objective cost from shoot k to total cost 
            problem.cost = problem.cost + Q;
            problem.resVert = problem.resVert + R;
            
        end
        
        % Iterate input counter and create variable for input control at
        % time step k
        vS.params.inputCountS = vS.params.inputCountS + 1;
        
        
        if (k ~= (vS.params.shootCount + vS.params.N - 1))
            
            
            % Create state variables for shoot k+1
            [vS, problem] = varCreation(vS, problem, 'state', k+1, 6,...
                vS.init.state(k+1,:)', vS.lb.state, vS.ub.state);
            stateK = vS.var.(['state' num2str(k+1)]);

            % Constraint: State at (k+1) must start where state at k ends
            [problem] = addConstraint(problem, stateEnd - stateK,...
                zeros(6,1), zeros(6,1));
            
            if (k == (vS.params.shootCount + vS.params.N - 2))
            
                [vS, problem] = varCreation(vS, problem, 'kSLagDot',...
                    vS.params.inputCountS, 1,...
                    vS.init.kLagDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);
                [vS, problem] = varCreation(vS, problem, 'kSLeadDot',...
                    vS.params.inputCountS, 1,...
                    vS.init.kLeadDot(k+1), vS.lb.kDotSS, vS.ub.kDotSS);
                
            else
                
                [vS, problem] = varCreation(vS, problem, 'kSLagDot',...
                    vS.params.inputCountS, 1,...
                    vS.init.kLagDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);
                [vS, problem] = varCreation(vS, problem, 'kSLeadDot',...
                    vS.params.inputCountS, 1,...
                    vS.init.kLeadDot(k+1), vS.lb.kDotSS, vS.ub.kDotSS);
                
            end
        
            if ~strcmp(vS.params.springType, 'Constant')

                if (k ~= (vS.params.shootCount + vS.params.N - 2))

                    [problem] = addConstraint(problem,...
                        vS.var.(['kSLagDot' num2str(vS.params.inputCountS)]) -...
                        vS.var.(['kSLagDot' num2str(vS.params.inputCountS-1)]),...
                        -vS.params.stiffVaryMax, vS.params.stiffVaryMax);
%                     [problem] = addConstraint(problem,...
%                         vS.var.(['kSLeadDot' num2str(vS.params.inputCountS)]) -...
%                         vS.var.(['kSLeadDot' num2str(vS.params.inputCountS-1)]),...
%                         -vS.params.stiffVaryMax, vS.params.stiffVaryMax);

                end
                
                [problem] = addConstraint(problem,...
                    vS.var.(['kSLeadDot' num2str(vS.params.inputCountS)]) -...
                    vS.var.(['kSLeadDot' num2str(vS.params.inputCountS-1)]),...
                    -vS.params.stiffVaryMax, vS.params.stiffVaryMax);

%             else
% 
%                 [problem] = addConstraint(problem,...
%                     vS.var.(['kSLagDot' num2str(vS.params.inputCountS)]) -...
%                     vS.var.(['kDLagDot' num2str(vS.params.inputCountD)]),...
%                     -vS.params.stiffVaryMax, vS.params.stiffVaryMax);
%                 [problem] = addConstraint(problem,...
%                     vS.var.(['kSLeadDot' num2str(vS.params.inputCountS)]) -...
%                     vS.var.(['kDLeadDot' num2str(vS.params.inputCountD)]),...
%                     -vS.params.stiffVaryMax, vS.params.stiffVaryMax);

            end
            
        end
        
%         if k == (vS.params.shootCount + vS.params.N - 1)
%             
%             [problem] = addConstraint(problem,...
%                 stateK(5) - stateEnd(6), -vS.params.stiffPhaseMax, vS.params.stiffPhaseMax);
%             [problem] = addConstraint(problem,...
%                 stateK(6) - stateEnd(5), -vS.params.stiffPhaseMax, vS.params.stiffPhaseMax);
%             
%         end

    end
    
    vS.params.shootCount = k + 1;
    
end