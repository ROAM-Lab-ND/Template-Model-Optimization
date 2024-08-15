%% %%%%%%%%%%%%%%%%%%%%%% templateModel_HumanFit %%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 17 December 2021
% Last Updated: 15 November 2023

% This code uses the CasADi framework to model a walking gait based on
% the Virtual Pivot Point (VPP) or Bipedal Spring Loaded Inverted Pendulum
% (BSLIP) template.

function templateModelOptCRC(dataStruct, tempParamsSetup, tempParamsVar)

    % Import CasADi functions
    import casadi.*

    % Initialize variable and dynamics structures
    varStruct = [];
    dynStruct = [];

    % If using collocation, initialize collocation points
    if strcmp(tempParamsVar{16}, 'Collocation')

        % Initialize dynStruct struct
        dynStruct = colPoints('radau', tempParamsVar{9});

    end
    
    % Temporarily store current optimization run and max number of runs
    tempJ = 1;
    tempRuns = tempParamsSetup{2};

    % Run optimization framework desired number of times
    while (tempJ <= tempRuns)
        
        %% %%%%%%%%%%%%%%% Dynamics and Structure Initialization %%%%%%%%%%%%%%% %%

        % Initial storage in varStruct
        [varStruct, problem] = varInit(varStruct, dynStruct,...
            dataStruct, tempParamsVar);

        % Storage of dynamics and CasADi functions in dynStruct struct
        [dynStruct, problem] = dynInit(varStruct, dynStruct, problem);

        % Initial state variable boundaries based on model type
        switch varStruct.params.model

            case 'VPP'

                % Initial bounds
                tempInitLB = [0; varStruct.lb.state(2); varStruct.lb.state(3);...
                    varStruct.params.xInit-0.2; varStruct.lb.state(5:6);...
                    10; varStruct.lb.state(8)];

                tempInitUB = [0; varStruct.ub.state(2); varStruct.ub.state(3);...
                    varStruct.params.xInit+0.2; -0.01; varStruct.ub.state(6:8)];

            case 'BSLIP'

                % Initial bounds
                tempInitLB = [0; varStruct.lb.state(2);...
                    varStruct.params.xInit-0.2; varStruct.lb.state(4);...
                    10; varStruct.lb.state(6)];

                tempInitUB = [0; varStruct.ub.state(2);...
                    varStruct.params.xInit+0.2; -0.01; varStruct.ub.state(5:6)];

            otherwise

                tempMsg = 'Invalid template model chosen';
                error(tempMsg);

        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%% NLP Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        % Loop for progressing through collocation/shooting elements
        for tempi = 1:varStruct.params.steps

            %%%%%%%%%%%%%%%%%% Double Support (DS) Phase %%%%%%%%%%%%%%%%%%
            
            % Create state variable
            if (tempi == 1)

                [varStruct, problem] = varCreation(varStruct, problem, 'state',...
                    varStruct.params.shootCount,...
                    length(varStruct.init.state(1,:)),...
                    varStruct.init.state(1,:)', tempInitLB, tempInitUB);

                % Constraint: Start lag leg foot at x<0m
                [problem] = addConstraint(problem, varStruct.var.footDLag1,...
                    varStruct.lb.foot, 0);

                if strcmp(varStruct.params.springType, 'Constant')

                    [problem] = addConstraint(problem,...
                        varStruct.var.(['state'...
                        num2str(varStruct.params.shootCount)])(end) -...
                        varStruct.var.(['state'...
                        num2str(varStruct.params.shootCount)])(end-1), -0.001, 0.001);

                end

            else

                [varStruct, problem] = varCreation(varStruct, problem, 'state',...
                    varStruct.params.shootCount,...
                    length(varStruct.init.state(1,:)),...
                    varStruct.init.state(varStruct.params.shootCount,:)',...
                    varStruct.lb.stateDS, varStruct.ub.stateDS);

            end

            % Constraint: Next SS lead leg becomes DS lag leg
            [problem] = addConstraint(problem,...
                varStruct.var.(['footDLag' num2str(tempi+1)]) -...
                varStruct.var.(['footSStance' num2str(tempi)]), 0, 0);

            % Store state in temporary variable
            tempStateK = varStruct.var.(['state'...
                num2str(varStruct.params.shootCount)]);

            % Establish location of hip based on model
            if strcmp(varStruct.params.model, 'VPP')

                % Calculate hip location
                tempXHip = tempStateK(1) - varStruct.params.rH*sin(tempStateK(3));
                tempZHip = tempStateK(2) - varStruct.params.rH*cos(tempStateK(3));

            elseif strcmp(varStruct.params.model, 'BSLIP')

                % Calculate hip location
                tempXHip = tempStateK(1);
                tempZHip = tempStateK(2);

            end

            % Store DS initial state variable
            varStruct.var.(['state0D' num2str(tempi)]) = tempStateK;

            [varStruct, problem] = varCreation(varStruct, problem, 'kDLagDot',...
                varStruct.params.inputCountD, 1,...
                varStruct.init.kLagDot(varStruct.params.inputCountD),...
                varStruct.lb.kDotDS, varStruct.ub.kDotDS);

            [varStruct, problem] = varCreation(varStruct, problem, 'kDLeadDot',...
                varStruct.params.inputCountD, 1,...
                varStruct.init.kLagDot(varStruct.params.inputCountD),...
                varStruct.lb.kDotDS, varStruct.ub.kDotDS);

            if (tempi ~= 1)

                if ~strcmp(varStruct.params.springType, 'Constant')

                    [problem] = addConstraint(problem, tempStateEnd(end) - tempStateK(end-1),...
                        0, 0);

                    [problem] = addConstraint(problem, tempStateEnd(end-1) - tempStateK(end),...
                        0, 0);

                    [problem] = addConstraint(problem, tempStateEnd(1:end-2) - tempStateK(1:end-2),...
                        zeros(length(varStruct.init.state(1,:))-2,1),...
                        zeros(length(varStruct.init.state(1,:))-2,1));

                    [problem] = addConstraint(problem,...
                        varStruct.var.(['kDLagDot' num2str(varStruct.params.inputCountD)]) -...
                        varStruct.var.(['kSLeadDot' num2str(varStruct.params.inputCountS-1)]),...
                        -varStruct.params.stiffVaryMax, varStruct.params.stiffVaryMax);

                else

                    [problem] = addConstraint(problem, tempStateEnd - tempStateK,...
                        zeros(length(varStruct.init.state(1,:)),1),...
                        zeros(length(varStruct.init.state(1,:)),1));

                end

                % Constraint: vertical velocity negative
                [problem] = addConstraint(problem, tempStateK((length(tempStateK)/2 + 1)), -0.5, -0.005);

                % Constraint: Next DS starts when lag leg touches down
                [problem] = addConstraint(problem, tempZHip -...
                    varStruct.var.len0*...
                    cos(varStruct.var.(['thetaD' num2str(tempi+1)])), 0, 0);

            end

            % Constraint: Start step at instant of touchdown
            [problem] = addConstraint(problem,...
                varStruct.var.len0 -...
                sqrt((tempZHip*...
                tan(varStruct.var.(['thetaD' num2str(tempi)])))^2 +...
                tempZHip^2), 0, 0);

            % Constraint: SS lead leg x position should be at initial touchdown of DS lead leg
            [problem] = addConstraint(problem,...
                varStruct.var.(['footSStance' num2str(tempi)]) -...
                tempXHip -...
                tempZHip.*tan(varStruct.var.(['thetaD' num2str(tempi)])),...
                0, 0);

            % Call into DS dynamics based on model
            switch varStruct.params.model

                case 'VPP'

                    [varStruct, problem, tempStateEnd] =...
                        nlpVPPDS(varStruct, dynStruct, problem, tempi, tempStateK);

                case 'BSLIP'

                    [varStruct, problem, tempStateEnd] =...
                        nlpBSLIPDS(varStruct, dynStruct, problem, tempi, tempStateK);

            end

            %%%%%%%%%%%%%%%%%% Single Support (SS) Phase %%%%%%%%%%%%%%%%%%
            
            % Update time tracking variable
            varStruct.params.timeTrack = varStruct.params.timeTrack +...
                varStruct.var.(['tfD' num2str(tempi)]);

            % Create state variable
            [varStruct, problem] = varCreation(varStruct, problem, 'state',...
                varStruct.params.shootCount,...
                length(varStruct.init.state(1,:)),...
                varStruct.init.state(varStruct.params.shootCount,:)',...
                varStruct.lb.state, varStruct.ub.state);

            % Store state in temporary variable
            tempStateK = varStruct.var.(['state'...
                num2str(varStruct.params.shootCount)]);

            [problem] = addConstraint(problem, tempStateEnd - tempStateK,...
                zeros(length(varStruct.init.state(1,:)),1),...
                zeros(length(varStruct.init.state(1,:)),1));

            [varStruct, problem] = varCreation(varStruct, problem, 'kSLagDot',...
                varStruct.params.inputCountS, 1,...
                varStruct.init.kLagDot(varStruct.params.inputCountS),...
                varStruct.lb.kDotSS, varStruct.ub.kDotSS);

            [varStruct, problem] = varCreation(varStruct, problem, 'kSLeadDot',...
                varStruct.params.inputCountS, 1,...
                varStruct.init.kLeadDot(varStruct.params.inputCountS),...
                varStruct.lb.kDotSS, varStruct.ub.kDotSS);

            if ~strcmp(varStruct.params.springType, 'Constant')

                [problem] = addConstraint(problem,...
                    varStruct.var.(['kSLagDot' num2str(varStruct.params.inputCountS)]) -...
                    varStruct.var.(['kDLagDot' num2str(varStruct.params.inputCountD-1)]),...
                    -varStruct.params.stiffVaryMax, varStruct.params.stiffVaryMax);

                [problem] = addConstraint(problem,...
                    varStruct.var.(['kSLeadDot' num2str(varStruct.params.inputCountS)]) -...
                    varStruct.var.(['kDLeadDot' num2str(varStruct.params.inputCountD-1)]),...
                    -varStruct.params.stiffVaryMax, varStruct.params.stiffVaryMax);

            end

            % Establish location of hip based on model
            if strcmp(varStruct.params.model, 'VPP')

                % Calculate hip location
                tempXHip = tempStateK(1) - varStruct.params.rH*sin(tempStateK(3));
                tempZHip = tempStateK(2) - varStruct.params.rH*cos(tempStateK(3));

            elseif strcmp(varStruct.params.model, 'BSLIP')

                % Calculate hip location
                tempXHip = tempStateK(1);
                tempZHip = tempStateK(2);

            end

            % Constraint: SS starts when lag leg reaches rest length
            [problem] = addConstraint(problem,...
                varStruct.var.len0 -...
                sqrt((tempXHip -...
                varStruct.var.(['footDLag' num2str(tempi)])).^2 +...
                tempZHip.^2), 0, 0);

            % Constraint: vertical velocity should be positive
            [problem] = addConstraint(problem, tempStateK((length(tempStateK)/2 + 1)), 0.005, 0.5);

            % Call into SS dynamics based on model
            switch varStruct.params.model

                case 'VPP'

                    [varStruct, problem, tempStateEnd] =...
                        nlpVPPSS(varStruct, problem, dynStruct, tempi, tempStateK);

                case 'BSLIP'

                    [varStruct, problem, tempStateEnd] =...
                        nlpBSLIPSS(varStruct, problem, dynStruct, tempi, tempStateK);

            end

            % Update time tracking variable
            varStruct.params.timeTrack = varStruct.params.timeTrack +...
                varStruct.var.(['tfS' num2str(tempi)]);

        end

        % Constraint: Time tracking variable should be close to human data time
        [problem] = addConstraint(problem, varStruct.params.timeTrack,...
            dataStruct.data.timeGait - 0.005,...
            dataStruct.data.timeGait + 0.005);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLP Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        % Call to solver, update flag to denote first run is complete
        [varStruct, dataStruct] = nlpSolver(varStruct, dataStruct,...
            dynStruct, problem);
        varStruct.params.run1 = false;
        tempRun1 = false;
        tempJ = tempJ + 1;

        % If objective cost is not low enough and not at run limit, rerun
        % by seeding with previous optimization result
        if (varStruct.optims.obj > 1e-4) && (tempRuns < 5)

            tempRuns = tempRuns + 1;
            varStruct.params.runs = tempRuns;
            
            tempInputCountD = 1;
            tempInputCountS = 1;
            tempShootCount = 1;

            tempParamsVar{10} = tempRun1;

        end

    end

    % Remove CasADi variable field to save on space
    varStruct = rmfield(varStruct, 'var');
    varStruct.params = rmfield(varStruct.params, 'timeTrack');

    % Save data
    saveData(varStruct, dataStruct)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Post Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    % Clear ALL temporary variables if flag is set
    if tempParamsSetup{3}

        clear temp*

    end

end