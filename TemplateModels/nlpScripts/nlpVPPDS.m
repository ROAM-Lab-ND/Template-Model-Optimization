%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nlpVPPDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 2 September 2021
% Last Updated: 12 November 2021

% This function is used to simulate and build the double support phase of
% the VPP model

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem, stateEnd] = nlpVPPDS(vS, dynStruct, problem, i, stateK)
    
    import casadi.*

    % Double Support Phase
    for k = vS.params.shootCount:(vS.params.shootCount + vS.params.N - 1)
        
        statePrev = stateK;
        vS.params.phase = 'DS';
        
        [vS, problem] = varCreation(vS, problem, 'u',...
            vS.params.inputCountD, 1,...
            vS.init.u(k-((i-1)*vS.params.N)), vS.lb.u, vS.ub.u);
        
        tDat = vS.params.timeTrack +...
            (((k) - 2*(i-1)*vS.params.N)/vS.params.N)*...
            vS.var.(['tfD' num2str(i)]);

        % Store non-state/non-control parameters in P vector for passing to
        % Collocation propagation
        P = [vS.var.(['state0D' num2str(i)]);...
            vS.var.(['thetaD' num2str(i)]);...
            vS.var.(['kDLagDot' num2str(vS.params.inputCountD)]);...
            vS.var.(['kDLeadDot' num2str(vS.params.inputCountD)]);...
            vS.var.(['footDLag' num2str(i)]);...
            vS.var.(['tfD' num2str(i)]); tDat; vS.var.len0];
        
        if strcmp(vS.params.vppType, 'Varying')
            
            P = [P; vS.var.(['rVPPD' num2str(i)])];
        
        end
        
        % Call to either CasADi function or Collocation matrices
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
            [vS, problem] = varCreation(vS, problem, 'state', k+1, 8,...
                vS.init.state(k+1,:)', vS.lb.state, vS.ub.state);

            stateK = vS.var.(['state' num2str(k+1)]);
            
            [problem] = addConstraint(problem, stateEnd - stateK,...
                zeros(8,1), zeros(8,1));
            
            [vS, problem] = varCreation(vS, problem, 'kDLagDot',...
                vS.params.inputCountD, 1,...
                vS.init.kLagDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);
            [vS, problem] = varCreation(vS, problem, 'kDLeadDot',...
                vS.params.inputCountD, 1,...
                vS.init.kLeadDot(k+1), vS.lb.kDotDS, vS.ub.kDotDS);
            
            xHip = stateK(1) - vS.params.rH*sin(stateK(3));
            yHip = stateK(2) - vS.params.rH*cos(stateK(3));

            lenLag = sqrt((xHip - vS.var.(['footDLag' num2str(i)]))^2 +...
                yHip^2);

            [problem] = addConstraint(problem,...
                1000*(stateK(7)*(vS.var.len0 - lenLag) +...
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

            end
            
        end
        
    end
    
    vS.params.shootCount = k + 1;

end