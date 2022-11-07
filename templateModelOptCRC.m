%% %%%%%%%%%%%%%%%%%%%%%% templateModel_HumanFit %%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 17 December 2021
% Last Updated: 19 May 2022

% This code uses the CasADi framework to model a walking gait based on
% the Virtual Pivot Point (VPP) or Bipedal Spring Loaded Inverted Pendulum 
% (BSLIP) template.

% INPUTS:
%   datS - Structure used for subject data storage
%   tempParamsSetup - Cell of optimization configuration options
%   tempParamsVar - Cell of variable optimization options

% OUTPUTS:
%   NONE

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function templateModelOptCRC(datS, tempParamsSetup, tempParamsVar)

    % Import CasADi toolbox
    import casadi.*

    % Initialize variable and dynamic structs
    varS = [];
    dynS = [];

    % Check if using collocation or multishooting
    if strcmp(tempParamsVar{16}, 'Collocation')

        % Grab collocation information
        dynS = colPoints('radau', tempParamsVar{9});

    end

    % Create run tracking variables
    tempJ = 1;
    tempRuns = tempParamsSetup{2};

    % While current run iterator is below desired number of runs
    while (tempJ <= tempRuns)
        
        %%%%%%%%%%%%%% Dynamics and Structure Initialization %%%%%%%%%%%%%%

        % Build variable structure
        [varS, problem] = varInit(varS, dynS, datS, tempParamsVar);

        % Build dynamics structure
        [dynS, problem] = dynInit(varS, dynS, problem);

        % Check if using VPP or BSLIP template model
        switch varS.params.model

            case 'VPP'

                % Define initial state value bounds
                tempInitLB = [0; varS.lb.state(2); varS.lb.state(3);...
                    varS.params.xInit-0.2; varS.lb.state(5:8)];

                tempInitUB = [0; varS.ub.state(2); varS.ub.state(3);...
                    varS.params.xInit+0.2; -0.01; varS.ub.state(6:8)];

            case 'BSLIP'

                % Define initial state value bounds
                tempInitLB = [0; varS.lb.state(2);...
                    varS.params.xInit-0.2; varS.lb.state(4:6)];

                tempInitUB = [0; varS.ub.state(2);...
                    varS.params.xInit+0.2; -0.01; varS.ub.state(5:6)];

            otherwise

                tempMsg = 'Invalid template model chosen';
                error(tempMsg);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%% NLP Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%

        % For each step to be taken by template model
        for tempi = 1:varS.params.steps

            %%%%%%%%%%%%%%%%%%%% Double Support Phase %%%%%%%%%%%%%%%%%%%%%
            
            % Check if current step is first step
            if (tempi == 1)

                % Create initial state variables
                [varS, problem] = varCreation(varS, problem, 'state',...
                    varS.params.shootCount,...
                    length(varS.init.state(1,:)),...
                    varS.init.state(1,:)', tempInitLB, tempInitUB);

                % Constraint: Start with lag leg behind origin
                [problem] = addConstraint(problem, varS.var.footDLag1,...
                    varS.lb.foot, 0);

                % Check if using constant or varying spring stiffness
                if strcmp(varS.params.springType, 'Constant')

                    % Constraint: Stiffness must be same between legs
                    [problem] = addConstraint(problem,...
                        varS.var.(['state'...
                        num2str(varS.params.shootCount)])(end) -...
                        varS.var.(['state'...
                        num2str(varS.params.shootCount)])(end-1),...
                        -0.001, 0.001);

                end

            else

                % Create state variables
                [varS, problem] = varCreation(varS, problem, 'state',...
                    varS.params.shootCount,...
                    length(varS.init.state(1,:)),...
                    varS.init.state(varS.params.shootCount,:)',...
                    varS.lb.state, varS.ub.state);

            end

            % Constraint: DS lag leg must align with SS leg
            [problem] = addConstraint(problem,...
                varS.var.(['footDLag' num2str(tempi+1)]) -...
                varS.var.(['footSStance' num2str(tempi)]), 0, 0);

            % Grab current state for easier coding
            tempStateK = varS.var.(['state'...
                num2str(varS.params.shootCount)]);

            % Check if using VPP or BSLIP model
            if strcmp(varS.params.model, 'VPP')

                % Calculate hip location
                tempXHip = tempStateK(1) -...
                    varS.params.rH*sin(tempStateK(3));
                tempZHip = tempStateK(2) -...
                    varS.params.rH*cos(tempStateK(3));

            elseif strcmp(varS.params.model, 'BSLIP')

                % Calculate hip location
                tempXHip = tempStateK(1);
                tempZHip = tempStateK(2);

            end

            % Store DS state variable as initial state for DS
            varS.var.(['state0D' num2str(tempi)]) = tempStateK;

            % Create leg stiffness ROC variables
            [varS, problem] = varCreation(varS, problem, 'kDLagDot',...
                varS.params.inputCountD, 1,...
                varS.init.kLagDot(varS.params.inputCountD),...
                varS.lb.kDot, varS.ub.kDot);

            [varS, problem] = varCreation(varS, problem, 'kDLeadDot',...
                varS.params.inputCountD, 1,...
                varS.init.kLagDot(varS.params.inputCountD),...
                varS.lb.kDot, varS.ub.kDot);

            % Check if on first step of gait cycle
            if (tempi ~= 1)

                % Check if using constant or varying spring stiffness
                if ~strcmp(varS.params.springType, 'Constant')

                    % Constraint: Leg stiffnesses must equal each other
                    [problem] = addConstraint(problem,...
                        tempStateEnd(end) - tempStateK(end-1), 0, 0);

                    [problem] = addConstraint(problem,...
                        tempStateEnd(end-1) - tempStateK(end), 0, 0);

                    % Constraint: State propogation must be continuous
                    [problem] = addConstraint(problem,...
                        tempStateEnd(1:end-2) - tempStateK(1:end-2),...
                        zeros(length(varS.init.state(1,:))-2,1),...
                        zeros(length(varS.init.state(1,:))-2,1));

                    % Constraint: Leg stiffness ROC must maintain bounds
                    % between phases (SS leg becomes DS lag leg)
                    [problem] = addConstraint(problem,...
                        varS.var.(['kDLagDot'...
                        num2str(varS.params.inputCountD)]) -...
                        varS.var.(['kSLeadDot'...
                        num2str(varS.params.inputCountS-1)]),...
                        -varS.params.stiffVaryMax,...
                        varS.params.stiffVaryMax);

                else

                    % Constraint: State propogation must be continuous
                    [problem] = addConstraint(problem,...
                        tempStateEnd - tempStateK,...
                        zeros(length(varS.init.state(1,:)),1),...
                        zeros(length(varS.init.state(1,:)),1));

                end

                % Constraint: Vertical CoM velocity must be negative
                [problem] = addConstraint(problem,...
                    tempStateK((length(tempStateK)/2 + 1)), -0.5, -0.005);

                % Constraint: End SS phase when lead lag touches down
                [problem] = addConstraint(problem, tempZHip -...
                    varS.var.len0*...
                    cos(varS.var.(['thetaD' num2str(tempi+1)])), 0, 0);

            end

            % Constraint: Start DS phase when lead lag touches down
            [problem] = addConstraint(problem,...
                varS.var.len0 -...
                sqrt((tempZHip*...
                tan(varS.var.(['thetaD' num2str(tempi)])))^2 +...
                tempZHip^2), 0, 0);

            % Constraint: SS foot position must be at DS lead foot position
            [problem] = addConstraint(problem,...
                varS.var.(['footSStance' num2str(tempi)]) -...
                tempXHip -...
                tempZHip.*tan(varS.var.(['thetaD' num2str(tempi)])),...
                0, 0);

            % Check if using VPP or BSLIP model
            switch varS.params.model

                case 'VPP'

                    % Call into VPP DS dynamics
                    [varS, problem, tempStateEnd] =...
                        nlpVPPDS(varS, dynS, problem, tempi, tempStateK);

                case 'BSLIP'

                    % Call into BSLIP DS dynamics
                    [varS, problem, tempStateEnd] =...
                        nlpBSLIPDS(varS, dynS, problem, tempi, tempStateK);

            end

            %%%%%%%%%%%%%%%%%%%% Single Support Phase %%%%%%%%%%%%%%%%%%%%%
            
            % Update time tracking variable
            varS.params.timeTrack = varS.params.timeTrack +...
                varS.var.(['tfD' num2str(tempi)]);

            % Create state variables
            [varS, problem] = varCreation(varS, problem, 'state',...
                varS.params.shootCount,...
                length(varS.init.state(1,:)),...
                varS.init.state(varS.params.shootCount,:)',...
                varS.lb.state, varS.ub.state);

            % Grab state variables for easier coding
            tempStateK = varS.var.(['state'...
                num2str(varS.params.shootCount)]);

            % Constraint: State propogation must be continuous
            [problem] = addConstraint(problem,...
                tempStateEnd - tempStateK,...
                zeros(length(varS.init.state(1,:)),1),...
                zeros(length(varS.init.state(1,:)),1));

            % Create leg stiffness ROC variables
            [varS, problem] = varCreation(varS, problem, 'kSLagDot',...
                varS.params.inputCountS, 1,...
                varS.init.kLagDot(varS.params.inputCountS),...
                varS.lb.kDot, varS.ub.kDot);

            [varS, problem] = varCreation(varS, problem, 'kSLeadDot',...
                varS.params.inputCountS, 1,...
                varS.init.kLeadDot(varS.params.inputCountS),...
                varS.lb.kDot, varS.ub.kDot);

            % Check if using constant or varying spring stiffness
            if ~strcmp(varS.params.springType, 'Constant')

                % Constraint: Leg stiffness change must maintain bounds
                [problem] = addConstraint(problem,...
                    varS.var.(['kSLagDot'...
                    num2str(varS.params.inputCountS)]) -...
                    varS.var.(['kDLagDot'...
                    num2str(varS.params.inputCountD-1)]),...
                    -varS.params.stiffVaryMax, varS.params.stiffVaryMax);

                [problem] = addConstraint(problem,...
                    varS.var.(['kSLeadDot'...
                    num2str(varS.params.inputCountS)]) -...
                    varS.var.(['kDLeadDot'...
                    num2str(varS.params.inputCountD-1)]),...
                    -varS.params.stiffVaryMax, varS.params.stiffVaryMax);

            end

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate hip position
                tempXHip = tempStateK(1) -...
                    varS.params.rH*sin(tempStateK(3));
                tempZHip = tempStateK(2) -...
                    varS.params.rH*cos(tempStateK(3));

            elseif strcmp(varS.params.model, 'BSLIP')

                % Calculate hip position
                tempXHip = tempStateK(1);
                tempZHip = tempStateK(2);

            end

            % Constraint: DS phase ends when lag leg reaches nominal length
            [problem] = addConstraint(problem,...
                varS.var.len0 -...
                sqrt((tempXHip -...
                varS.var.(['footDLag' num2str(tempi)])).^2 +...
                tempZHip.^2), 0, 0);

            % Constraint: Vertical CoM velocity must be positive
            [problem] = addConstraint(problem,...
                tempStateK((length(tempStateK)/2 + 1)), 0.005, 0.5);

            % Check if using VPP or BSLIP template model
            switch varS.params.model

                case 'VPP'

                    % Call into SS VPP dynamics
                    [varS, problem, tempStateEnd] =...
                        nlpVPPSS(varS, problem, dynS, tempi, tempStateK);

                case 'BSLIP'

                    % Call into SS BSLIP dynamics
                    [varS, problem, tempStateEnd] =...
                        nlpBSLIPSS(varS, problem, dynS, tempi, tempStateK);

            end

            % Update time tracking variable
            varS.params.timeTrack = varS.params.timeTrack +...
                varS.var.(['tfS' num2str(tempi)]);

        end

        % Constraint: Time tracking must be close to human gait duration
        [problem] = addConstraint(problem, varS.params.timeTrack,...
            datS.data.timeGait - 0.005,...
            datS.data.timeGait + 0.005);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% NLP Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Call to NLP solver
        [varS, datS] = nlpSolver(varS, datS, dynS, problem);
        
        % Set first run flag to false
        varS.params.run1 = false;
        tempRun1 = false;
        
        % Iterate number of optimization runs for trial
        tempJ = tempJ + 1;

        % Check if objective cost and trial run thresholds have been met
        if (varS.optims.obj > 1e-4) && (tempRuns < 5)

            % Iterate number of trial runs to conduct
            tempRuns = tempRuns + 1;
            varS.params.runs = tempRuns;

            % Reset element iterators for each phase
            tempInputCountD = 1;
            tempInputCountS = 1;
            tempShootCount = 1;

            % Reset number of trial runs
            tempParamsVar{10} = tempRun1;

        end


    end

    % Clean up variable structure for future data storage
    varS = rmfield(varS, 'var');
    varS.params = rmfield(varS.params, 'timeTrack');

    % Call to save data for future analysis
    saveData(varS, datS)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Post Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if temporary variable clear flag is set
    if tempParamsSetup{3}

        % Clear all variables that start with 'temp...'
        clear temp*

    end

end