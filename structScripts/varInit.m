%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% varInit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 09 December 2021
% Last Updated: 04 May 2022

% This function is used to initialize the variable (varS) and NLP (problem)
% structures that are used for template model optimization. Initial bounds,
% parameters, and variables/constraint storage is initialized by calling
% this function.

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   dynS - Structure used for dynamics storage
%   datS - Structure used for subject data storage
%   params - Cell containing different initialization options

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   problem - Structure used to pass necessary information to NLP solver

%% %%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, problem] = varInit(varS, dynS, datS, params)

    % Initialize problem struct
    problem.vars = {};
    problem.varsInit = [];
    problem.varsLB = [];
    problem.varsUB = [];

    problem.cost = 0;
    problem.resVert = 0;

    problem.constraints = {};
    problem.constraintsLB = [];
    problem.constraintsUB = [];

    % Initialize variable struct
    varS.inds.state = [];
    varS.inds.stateFull = [];
    varS.inds.foot = [];
    varS.inds.theta = [];
    varS.inds.u = [];
    varS.inds.kLagDot = [];
    varS.inds.kLeadDot = [];
    varS.inds.time = [];
    varS.inds.gen = [];

    varS.params.model = params{15};
    varS.params.method = params{16};
    varS.params.size = 'Base';
    varS.params.source = params{17};
    varS.params.springType = params{18};
    varS.params.fitType = params{20};

    varS.params.len0 = datS.subj.len0;
    varS.params.m = datS.subj.weight;
    varS.params.g = datS.subj.g;

    % Check if model type is VPP to initialize VPP specific variables
    if strcmp(varS.params.model, 'VPP')

        varS.params.vppType = params{19};

        % Check if VP parameter is set to varying, initialize VP variable
        if strcmp(params{19}, 'Varying')

            varS.inds.rVPP = [];

        end

        varS.params.rH = params{3};

        % Check if VP parameter is constant, initialize VP variable
        if strcmp(varS.params.vppType, 'Constant')

            varS.params.rVPPD = params{4}(1);
            varS.params.rVPPS = params{4}(2);

        end

        varS.params.gamma = params{2};
        varS.params.J = params{5};

    end

    varS.params.xInit = datS.subj.treadmill;
    varS.params.vDiff = params{6};
    varS.params.N = params{8};
    varS.params.M = params{9};
    varS.params.dt = 1/(params{8}*params{9});
    varS.params.steps = params{1};
    varS.params.inputCountD = params{11};
    varS.params.inputCountS = params{12};
    varS.params.shootCount = params{13};
    varS.params.timeTrack = params{14};
    varS.params.Q = params{21};
    varS.params.residualVert = 0;

    % Check if model stiffness is constant, initialize varying properties
    if strcmp(varS.params.springType, 'Constant')

        varS.params.stiffVaryMax = 0;
        varS.params.stiffPhaseMax = 0;

    else

        varS.params.stiffVaryMax = 10;
        varS.params.stiffPhaseMax = 0.05;

    end

    % Check model type for VPP or BSLIP
    switch params{15}

        case 'VPP'

            % Set upper and lower bounds for state variables
            varS.lb.state = [0; datS.subj.len0*cos(pi/6); -15*pi/180;...
                datS.subj.treadmill-0.5; -0.5; -100*pi/180;...
                5; 5];

            varS.ub.state = [2*params{1};...
                1.5*(datS.subj.len0+params{3}); 15*pi/180;...
                datS.subj.treadmill+0.5; 0.5; 100*pi/180; 30; 30];

        case 'BSLIP'

            % Set upper and lower bounds for state variables
            varS.lb.state = [0; datS.subj.len0*cos(pi/6);...
                datS.subj.treadmill-0.5; -0.5; 5; 5];

            varS.ub.state = [2*params{1}; 2*datS.subj.len0;...
                datS.subj.treadmill+0.5; 0.5; 50; 50];

        otherwise

            errMsg = 'Inputted Template Model Not Yet Implemented';
            error(errMsg);

    end

    % Check if there is more than one gait cycle being evaluated
    if length(datS.data.timeDS1)>1

        % NEEDS UPDATING
        % Set upper and lower bounds for double support phases
        varS.lb.tfD = min(min(datS.data.timeDS1),...
            min(datS.data.timeDS2)) - 0.02;

        varS.ub.tfD = max(max(datS.data.timeDS1),...
            max(datS.data.timeDS2)) + 0.02;

        % Set upper and lower bounds for single support phases
        varS.lb.tfS = min(min(datS.data.timeSS1),...
            min(datS.data.timeSS2)) - 0.02;

        varS.ub.tfS = max(max(datS.data.timeSS1),...
            max(datS.data.timeSS2)) + 0.02;

    else

        % Set upper and lower bounds for double support phase
        varS.lb.tfD = min(datS.data.timeDS1, datS.data.timeDS2) - 0.02;

        varS.ub.tfD = max(datS.data.timeDS1, datS.data.timeDS2) + 0.02;

        % Set upper and lower bounds for single support phase
        varS.lb.tfS = min(datS.data.timeSS1, datS.data.timeSS2) - 0.02;

        varS.ub.tfS = max(datS.data.timeSS1, datS.data.timeSS2) + 0.02;

    end

    % Set upper and lower bound for foot positions
    varS.lb.foot = -1;
    varS.ub.foot = 5;

    % Set upper and lower bound for touchdown angle
    varS.lb.theta = 5*pi/180;
    varS.ub.theta = 45*pi/180;

    % Check if model stiffness is constant, set upper and lower ROC bounds
    if strcmp(varS.params.springType, 'Constant')

        varS.lb.kDot = 0;
        varS.ub.kDot = 0;

    else

        varS.lb.kDot = -100;
        varS.ub.kDot = 100;

    end

    % Set upper and lower control bounds
    varS.lb.u = 0;
    varS.ub.u = 0;

    % Set upper and lower epsilon bounds
    varS.lb.eps = -params{7};
    varS.ub.eps = params{7};

    % Set upper and lower leg length bounds
    varS.lb.leg = datS.subj.len0 - 0.2;
    varS.ub.leg = datS.subj.len0 + 0.2;

    % Check if VP is constant or varying, set upper and lower bounds
    if strcmp(params{19}, 'Varying')

        varS.lb.rVPP = 0;
        varS.ub.rVPP = 0.5;

    end

    % Store coefficients of vertical CoM trajectory wrt time
    varS.fit.humPosVert = datS.fit.coefPosVert;
    varS.fit.humPosVertF = datS.fit.coefPosVertF;

    % Store coefficients of horizontal CoM trajectory wrt time
    varS.fit.humPosHor = datS.fit.coefPosHor;
    varS.fit.humPosHorF = datS.fit.coefPosHorF;

    % Store coefficients of vertical CoM velocity wrt time
    varS.fit.humVelVert = datS.fit.coefVelVert;

    % Store coefficients of horizontal CoM velocity wrt time
    varS.fit.humVelHor = datS.fit.coefVelHor;

    % Check if first optimization run, initialize variable seeding
    if params{10}

        % Store number of optimization attempts
        varS.params.runs = 1;

        % Initialize time vector for subject data
        varS.fit.tHum = 0;

        % For each gait cycle
        for i = 1:length(datS.data.timeDS1)

            % Check if initial gait cycle
            if i == 1

                % Generate time vector for subject data
                varS.fit.tHum = [linspace(varS.fit.tHum(end),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i)),...
                    params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    (datS.data.timeSS1(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i)), params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) +...
                    (datS.data.timeDS2(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i)),...
                    params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                    (datS.data.timeSS2(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                    datS.data.timeSS2(i)), params{8})];

            else

                % Append to time vector for subject data
                varS.fit.tHum = [varS.fit.tHum,...
                    linspace(varS.fit.tHum(end),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i)),...
                    params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    (datS.data.timeSS1(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i)), params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) +...
                    (datS.data.timeDS2(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i)),...
                    params{8}),...
                    ...
                    linspace((varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                    (datS.data.timeSS2(i)/params{8})),...
                    (varS.fit.tHum(end) + datS.data.timeDS1(i) +...
                    datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                    datS.data.timeSS2(i)), params{8})];

            end

        end

        % Initialize full time vector for subject data
        varS.fit.tHumFull = 0;

        % Check if using collocation or multishooting
        if strcmp(varS.params.method, 'Collocation')

            % For each gait phase
            for i = 1:length(datS.data.timeFull)

                % For each shooting element
                for j = 1:params{8}

                    % Check if first phase and element
                    if (j == 1) && (i == 1)

                        % Generate full time vector
                        varS.fit.tHumFull = [(varS.fit.tHumFull(end) +...
                            dynS.colPoints(2:end)*...
                            (datS.data.timeFull(i)/params{8}))];

                    else

                        % Append to full time vector
                        varS.fit.tHumFull = [varS.fit.tHumFull,...
                            (varS.fit.tHumFull(end) +...
                            dynS.colPoints(2:end)*...
                            (datS.data.timeFull(i)/params{8}))];

                    end

                end

            end

        elseif strcmp(varS.params.method, 'Multishooting')

            % For each gait cycle
            for i = 1:length(datS.data.timeDS1)

                % Check if first gait cycle
                if i == 1

                    % Generate full time vector
                    varS.fit.tHumFull =...
                        [linspace(varS.fit.tHumFull(end),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i)),...
                        params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) +...
                        (datS.data.timeSS1(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i)), params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) + datS.data.timeSS1(i) +...
                        (datS.data.timeDS2(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i) + datS.data.timeDS2(i)),...
                        params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) + datS.data.timeSS1(i) +...
                        datS.data.timeDS2(i) +...
                        (datS.data.timeSS2(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                        datS.data.timeSS2(i)), params{9}*params{8})];

                else

                    % Append to full time vector
                    varS.fit.tHumFull = [varS.fit.tHumFull,...
                        linspace(varS.fit.tHumFull(end),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i)),...
                        params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) +...
                        (datS.data.timeSS1(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i)), params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) + datS.data.timeSS1(i) +...
                        (datS.data.timeDS2(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i) + datS.data.timeDS2(i)),...
                        params{9}*params{8}),...
                        ...
                        linspace((varS.fit.tHumFull(end) +...
                        datS.data.timeDS1(i) + datS.data.timeSS1(i) +...
                        datS.data.timeDS2(i) +...
                        (datS.data.timeSS2(i)/(params{9}*params{8}))),...
                        (varS.fit.tHumFull(end) + datS.data.timeDS1(i) +...
                        datS.data.timeSS1(i) + datS.data.timeDS2(i) +...
                        datS.data.timeSS2(i)), params{9}*params{8})];

                end


            end

        end

        % Store polyfit evaluation of CoM kinematics
        varS.fit.posCOMSegHor = polyval(varS.fit.humPosHor,varS.fit.tHum);
        varS.fit.posCOMSegVert = polyval(varS.fit.humPosVert,varS.fit.tHum);
        varS.fit.velCOMSegHor = polyval(varS.fit.humVelHor,varS.fit.tHum);
        varS.fit.velCOMSegVert = polyval(varS.fit.humVelVert,varS.fit.tHum);

        % Store Fourier evaluation of CoM kinematics
        varS.fit.posCOMSegHorF = fitFourier(varS.fit.humPosHorF,...
            varS.fit.tHum, 'Hor');
        varS.fit.posCOMSegVertF = fitFourier(varS.fit.humPosVert,...
            varS.fit.tHum, 'Vert');
        varS.fit.velCOMSegHorF = fitFourierDer(varS.fit.humPosHorF,...
            varS.fit.tHum, 'Hor');
        varS.fit.velCOMSegVertF = fitFourierDer(varS.fit.humPosVertF,...
            varS.fit.tHum, 'Vert');

        % Check if using VPP or BSLIP template model
        switch params{15}

            case 'VPP'

                % Initialize state variables seed
                varS.init.state = zeros(2*params{1}*params{8}, 8);
                varS.init.state(:,3) = 10*pi/180;

                % Check if using polyfit or Fourier fit
                if params{20} > 1

                    % Seed CoM velocity
                    varS.init.state(:,4) = varS.fit.velCOMSegHor;
                    varS.init.state(:,5) = varS.fit.velCOMSegVert;

                else

                    % Seed CoM velocity
                    varS.init.state(:,4) = varS.fit.velCOMSegHorF;
                    varS.init.state(:,5) = varS.fit.velCOMSegVertF;

                end

                % Check if VP is set to varying or constant
                if strcmp(params{19}, 'Varying')

                    % Seed VP
                    varS.init.rVPPS = zeros(params{1},1);
                    varS.init.rVPPD = 0.3*ones(params{1},1);

                end

            case 'BSLIP'

                % Initialize state variables seed
                varS.init.state = zeros(2*params{1}*params{8},6);

                % Check if using polyfit or Fourier fit
                if params{20} > 1

                    % Seed CoM velocity
                    varS.init.state(:,3) = varS.fit.velCOMSegHor;
                    varS.init.state(:,4) = varS.fit.velCOMSegVert;

                else

                    % Seed CoM velocity
                    varS.init.state(:,3) = varS.fit.velCOMSegHorF;
                    varS.init.state(:,4) = varS.fit.velCOMSegVertF;

                end


            otherwise

                disp('Non-implemented model chosen')

        end

        % Check if using polyfit or Fourier fit
        if params{20} > 1

            % Seed CoM position
            varS.init.state(:,1) = varS.fit.posCOMSegHor;
            varS.init.state(:,2) = varS.fit.posCOMSegVert;

        else

            % Seed CoM position
            varS.init.state(:,1) = varS.fit.posCOMSegHorF;
            varS.init.state(:,2) = varS.fit.posCOMSegVertF;

        end

        % Seed leg stiffness
        varS.init.state(:,end-1) = 20;
        varS.init.state(:,end) = 20;

        % Seed leg stiffness ROC
        varS.init.kLagDot = zeros(2*params{1}*params{8},1);
        varS.init.kLeadDot = zeros(2*params{1}*params{8},1);

        % Seed nominal leg length
        varS.init.leg = varS.params.len0;

        % Seed gait phase durations
        tempTfD = [datS.data.timeDS1'; datS.data.timeDS2'];
        tempTfS = [datS.data.timeSS1'; datS.data.timeSS2'];
        varS.init.tFD = tempTfD(:);
        varS.init.tFS = tempTfS(:);

        % Seed touchdown angle
        varS.init.theta = 20*pi/180*ones((params{1} + 1), 1);

        % Seed foot position
        varS.init.footD = 0.5*(0:params{1});
        varS.init.footS = (2/3) + (2/3)*(0:params{1}-1);

        % Seed input
        varS.init.u = zeros(params{1}*params{8},1);

    end

    % For each step taken by template model
    for i=1:params{1}

        %%%%%%%%%%%%%%%%%%%%%% Double support phase %%%%%%%%%%%%%%%%%%%%%%%
        
        % Create limit of integration variable
        [problem, varS.var.(['tfD' num2str(i)]),...
            varS.inds.time(end+1,:)] = addVariable(problem,...
            ['tfD' num2str(i)], 1, varS.init.tFD(i),...
            varS.lb.tfD, varS.ub.tfD);
        
        % Create touchdown angle variable
        [problem, varS.var.(['thetaD' num2str(i)]),...
            varS.inds.theta(end+1,:)] = addVariable(problem,...
            ['thetaD' num2str(i)], 1, varS.init.theta(i),...
            varS.lb.theta, varS.ub.theta);
        
        % Create foot location variable
        [problem, varS.var.(['footDLag' num2str(i)]),...
            varS.inds.foot(end+1,:)] = addVariable(problem,...
            ['footDLag' num2str(i)], 1, varS.init.footD(i),...
            varS.lb.foot, varS.ub.foot);

        %%%%%%%%%%%%%%%%%%%%%% Single Support Phase %%%%%%%%%%%%%%%%%%%%%%%
        
        % Create limit of integration variable
        [problem, varS.var.(['tfS' num2str(i)]),...
            varS.inds.time(end+1,:)] = addVariable(problem,...
            ['tfS' num2str(i)], 1, varS.init.tFS(i),...
            varS.lb.tfS, varS.ub.tfS);
        
        % Create foot location variable
        [problem, varS.var.(['footSStance' num2str(i)]),...
            varS.inds.foot(end+1,:)] = addVariable(problem,...
            ['footSStance' num2str(i)], 1, varS.init.footS(i),...
            varS.lb.foot, varS.ub.foot);

        %%%%%%%%%%%%%%%%%%%%%%%%%% VPP Specific %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Check if using VPP or BSLIP template model
        if strcmp(varS.params.model, 'VPP')

            % Check if using varying or constant VP
            if strcmp(varS.params.vppType, 'Varying')

                % Create double support VP variable
                [problem, varS.var.(['rVPPD' num2str(i)]),...
                    varS.inds.rVPP(end+1,:)] = addVariable(problem,...
                    ['rVPPD' num2str(i)], 1,...
                    varS.init.rVPPD(i), varS.lb.rVPP, varS.ub.rVPP);

                % Create single support VP variable
                [problem, varS.var.(['rVPPS' num2str(i)]),...
                    varS.inds.rVPP(end+1,:)] = addVariable(problem,...
                    ['rVPPS' num2str(i)], 1,...
                    varS.init.rVPPS(i), varS.lb.rVPP, varS.ub.rVPP);

            end

        end

    end

    % Create final touchdown angle variable
    [problem, varS.var.(['thetaD' num2str(i+1)]),...
        varS.inds.theta(end+1,:)] = addVariable(problem,...
        ['thetaD' num2str(i+1)], 1,...
        varS.init.theta(end), varS.lb.theta, varS.ub.theta);
    
    % Create final foot position variable
    [problem, varS.var.(['footDLag' num2str(i+1)]),...
        varS.inds.foot(end+1,:)] = addVariable(problem,...
        ['footDLag' num2str(i+1)], 1,...
        varS.init.footD(end), varS.lb.foot, varS.ub.foot);
    
    % Create nominal leg length variable
    [problem, varS.var.len0, varS.inds.leg(1)] = addVariable(problem,...
        'len0', 1, varS.init.leg, varS.lb.leg, varS.ub.leg);

end