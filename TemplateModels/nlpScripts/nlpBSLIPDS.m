%% %%%%%%%%%%%%%%%%%%%%%%%%% structInitBSLIP %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 17 November 2021

% This function is used to propagate the BSLIP model through the Double
% Support (DS) phase. Note that force input is currently only allowed
% during DS to mimick the act of ankle push off in bipedal walking.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem, stateEnd] = nlpBSLIPDS(vS, dynStruct, problem, i, stateK)
    
    import casadi.*

    % Double Support Phase
    for k = vS.params.shootCount:(vS.params.shootCount + vS.params.N - 1)
        
        statePrev = stateK;
        vS.params.phase = 'DS';
        
        [vS, problem] = varCreation(vS, problem, 'u',...
            vS.params.inputCountD, 1,...
            vS.init.u(k-((i-1)*vS.params.N)), vS.lb.u, vS.ub.u);
        
%         [vS, problem] = varCreation(vS, problem, 'kDLagDot',...
%             vS.params.inputCountD, 1,...
%             vS.init.kLagDot(k+1), vS.lb.kDot, vS.ub.kDot);
%         [vS, problem] = varCreation(vS, problem, 'kDLeadDot',...
%             vS.params.inputCountD, 1,...
%             vS.init.kLeadDot(k+1), vS.lb.kDot, vS.ub.kDot);
        
        tDat = vS.params.timeTrack +...
            (((k) - 2*(i-1)*vS.params.N)/vS.params.N)*...
            vS.var.(['tfD' num2str(i)]);

        % Store non-state/non-control parameters in P vector for passing to
        % RK4 integrator
        P = [vS.var.(['state0D' num2str(i)]);...
            vS.var.(['thetaD' num2str(i)]);...
            vS.var.(['kDLagDot' num2str(vS.params.inputCountD)]);...
            vS.var.(['kDLeadDot' num2str(vS.params.inputCountD)]);...
            vS.var.(['footDLag' num2str(i)]);...
            vS.var.(['tfD' num2str(i)]); tDat; vS.var.len0];
        
        if strcmp(vS.params.method, 'Multishooting')
            
            % Pass time, state, parameter, and control variables to integrator
            Fk = dynStruct.FD('state0',stateK,'p',P,'u',...
                vS.var.(['u' num2str(vS.params.inputCountD)]));
            
            % Store result of state after propogating shoot through integrator
            stateEnd = Fk.stated;

            % Add resulting objective cost from shoot k to total cost 
            problem.cost = problem.cost + Fk.qd;
            problem.resVert = problem.resVert + Fk.rd;
            
        elseif strcmp(vS.params.method, 'Collocation')

            params = {k , stateK , P , vS.var.(['u' num2str(vS.params.inputCountD)])};

            [vS, problem, stateEnd, Q, R] = colProp(vS, dynStruct, problem, params);
        
            % Add resulting objective cost from shoot k to total cost 
            problem.cost = problem.cost + Q;
            problem.resVert = problem.resVert + R;
            
        end
        
        % Iterate input counter and create variable for input control at
        % time step k
        vS.params.inputCountD = vS.params.inputCountD + 1;
        
        if (k ~= (vS.params.shootCount + vS.params.N - 1))
        
            % Create state variables for shoot k+1
            [vS, problem] = varCreation(vS, problem, 'state', k+1, 6,...
                vS.init.state(k+1,:)', vS.lb.state, vS.ub.state);
            
            stateK = vS.var.(['state' num2str(k+1)]);

            % Constraint: State at (k+1) must start where state at k ends
            [problem] = addConstraint(problem, stateEnd - stateK,...
                zeros(6,1), zeros(6,1));
            
            [vS, problem] = varCreation(vS, problem, 'kDLagDot',...
                vS.params.inputCountD, 1,...
                vS.init.kLagDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);
            [vS, problem] = varCreation(vS, problem, 'kDLeadDot',...
                vS.params.inputCountD, 1,...
                vS.init.kLeadDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);

        
            % Constraint: X-velocity of state must not overshoot desired final
            % velocity by more than epsilon
    %         [problem] = addConstraint(problem,...
    %             stateK(4) - vS.var.state0D1(4) - vS.params.vDiff,...
    %             vS.lb.eps, vS.ub.eps);

            xHip = stateK(1);
            yHip = stateK(2);

            lenLag = sqrt((xHip - vS.var.(['footDLag' num2str(i)]))^2 +...
                yHip^2);

            [problem] = addConstraint(problem,...
                1000*(stateK(5)*(vS.var.len0 - lenLag) +...
                vS.var.(['u' num2str(vS.params.inputCountD-1)])*...
                yHip/(lenLag)), 0, inf);

            if ~strcmp(vS.params.springType, 'Constant')

                [problem] = addConstraint(problem,...
                    vS.var.(['kDLagDot' num2str(vS.params.inputCountD)]) -...
                    vS.var.(['kDLagDot' num2str(vS.params.inputCountD-1)]),...
                    -vS.params.stiffVaryMax, vS.params.stiffVaryMax);
                [problem] = addConstraint(problem,...
                    vS.var.(['kDLeadDot' num2str(vS.params.inputCountD)]) -...
                    vS.var.(['kDLeadDot' num2str(vS.params.inputCountD-1)]),...
                    -vS.params.stiffVaryMax, vS.params.stiffVaryMax);

%                 [problem] = addConstraint(problem,...
%                     vS.var.(['kDLagDot' num2str(vS.params.inputCountD)]) -...
%                     vS.var.(['kSLeadDot' num2str(vS.params.inputCountS)]),...
%                     -vS.params.stiffVaryMax, vS.params.stiffVaryMax);
%                 [problem] = addConstraint(problem,...
%                     vS.var.(['kDLeadDot' num2str(vS.params.inputCountD)]) -...
%                     vS.var.(['kSLagDot' num2str(vS.params.inputCountS)]),...
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