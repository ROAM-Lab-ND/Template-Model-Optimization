%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nlpSolver %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 2 September 2021
% Last Updated: 2 September 2021

% This function is used to invoke the solver IPOPT after the CasADi
% dynamics model has already been built for the BSLIP and VPP models. Note
% that this function should be called for the case where the nominal leg
% length is optimized by CasADi

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage
%   dynS - Structure used for dynamic storage
%   problem - Structure used to pass necessary information to NLP solver

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, datS] = nlpSolver(varS, datS, dynS, problem)

    import casadi.*

    % Create prob struct to pass required information to NLP solver
    prob = struct('f', problem.cost, 'x', vertcat(problem.vars{:}),...
        'g', vertcat(problem.constraints{:}));

    % IPOPT Options Setting
    sOpts = struct;
    sOpts.ipopt.max_iter = 500;
    sOpts.ipopt.print_info_string = 'yes';
    sOpts.ipopt.print_level = 5;
    sOpts.ipopt.file_print_level = 4;
    % sOpts.ipopt.linear_solver = 'ma57';
    % sOpts.ipopt.hessian_approximation = 'limited-memory';
    % sOpts.ipopt.tol = 1e-02;
    
    % Store file location where results should be stored
    fileLoc = cd;
    fileLoc = fullfile(fileLoc, '\Data');
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')
        
        % Check if using constant VP and VP is different between DS and SS
        if strcmp(varS.params.vppType, 'Constant') &&...
                (abs(varS.params.rVPPD - varS.params.rVPPS) == 0)
            
            % Store file name for data storage
            fileName = [fileLoc '\s' datS.subj.subjectID '_t'...
                datS.subj.testNum '_i' num2str(datS.inds.start)...
                num2str(datS.inds.end) '_' varS.params.model '_'...
                varS.params.method '_' varS.params.springType 'K_'...
                varS.params.vppType 'SameVPP_Ipopt.txt'];
            
        else
            
            % Store file name for data storage
            fileName = [fileLoc '\s' datS.subj.subjectID '_t'...
                datS.subj.testNum '_i' num2str(datS.inds.start)...
                num2str(datS.inds.end) '_' varS.params.model '_'...
                varS.params.method '_' varS.params.springType 'K_'...
                varS.params.vppType 'DiffVPP_Ipopt.txt'];
            
        end
        
    else
        
        % Store file name for data storage
        fileName = [fileLoc '\s' datS.subj.subjectID '_t' datS.subj.testNum...
            '_i' num2str(datS.inds.start) num2str(datS.inds.end) '_'...
            varS.params.model '_' varS.params.method '_'...
            varS.params.springType 'K_Ipopt.txt'];
        
    end
    
    % Store filename in IPOPT options struct
    sOpts.ipopt.output_file = fileName;

    % Set solver to IPOPT
    solver = nlpsol('solver', 'ipopt', prob, sOpts);

    % Solve using NLP
    sol = solver('x0', problem.varsInit, 'lbx', problem.varsLB,...
        'ubx', problem.varsUB, 'lbg', problem.constraintsLB,...
        'ubg', problem.constraintsUB);

    % Grab resulting solution for all variables
    wOpt = full(sol.x);
    
    % Store resulting objective cost
    varS.optims.obj = full(sol.f);

    % Grab state values from wOpt, reseed initial guess values
    stateOpt = reshape(wOpt(varS.inds.state(:)), [],...
        length(varS.inds.state(1,:)));
    varS.optims.state = stateOpt;
    varS.init.state = stateOpt;
    
    % Check if collocation method was used, grab full state values
    if strcmp(varS.params.method, 'Collocation')
        
        stateFullOpt = reshape(wOpt(varS.inds.stateFull(:)),[],...
            length(varS.inds.stateFull(1,:)));
        varS.optims.stateFull = stateFullOpt;
        
    end

    % Grab foot position values from wOpt
    fOpt = reshape(wOpt(varS.inds.foot(:)), 2*varS.params.steps + 1, 1);
    varS.optims.fOpt = fOpt;
    varS.optims.fOptOrig = fOpt;
    varS.init.footD = fOpt(1:2:end);
    varS.init.footS = fOpt(2:2:end);

    % Grab touchdown angle values from wOpt, reseed initial guess values
    angOpt = reshape(wOpt(varS.inds.theta(:)), varS.params.steps + 1,1);
    varS.optims.angOpt = angOpt;
    varS.init.theta = angOpt;
    varS.optims.angOptDeg = angOpt*180/pi;
    
    % Grab input values from wOpt, reseed initial guess values
    uOpt = reshape(wOpt(varS.inds.u(:)), [], 1);
    varS.optims.uOpt = uOpt;
    varS.init.u = uOpt;

    % Grab leg stiffness values from wOpt
    kOptLag = stateOpt(:, end-1);
    kOptLead = stateOpt(:, end);
    varS.optims.kOptLag = kOptLag;
    varS.optims.kOptLead = kOptLead;
    
    % Grab leg stiffness ROC values from wOpt, reseed initial guess values
    kOptLagDot = reshape(wOpt(varS.inds.kLagDot(:)), [], 1);
    kOptLeadDot = reshape(wOpt(varS.inds.kLeadDot(:)), [] ,1);
    varS.optims.kOptLagDot = kOptLagDot;
    varS.optims.kOptLeadDot = kOptLeadDot;
    varS.init.kLagDot = kOptLagDot;
    varS.init.kLeadDot = kOptLeadDot;
    
    % Check if using VPP template model
    if strcmp(varS.params.model, 'VPP')
        
        % Check if using varying stiffness model
        if strcmp(varS.params.vppType, 'Varying')
        
            % Grab VP values from wOpt
            rVPP = reshape(wOpt(varS.inds.rVPP(:)), [], 1);
            varS.optims.rVPP = rVPP;
            
        else
            
            % Store VP values in format that aligns with varying VP method
            varS.optims.rVPP = repmat([varS.params.rVPPD;...
                varS.params.rVPPS], varS.params.steps, 1);
            
        end
        
        % Store VP values for each phase, reseed initial guess values
        varS.optims.rVPPD = varS.optims.rVPP(1:2:end);
        varS.optims.rVPPS = varS.optims.rVPP(2:2:end);
        varS.init.rVPPD = varS.optims.rVPPD;
        varS.init.rVPPS = varS.optims.rVPPS;
        
    end

    % Grab limit of integration variables from wOpt
    timeOpt = reshape(wOpt(varS.inds.time(:)),[],1);
    varS.optims.timeOpt = timeOpt;
    varS.init.tFD = timeOpt(1:2:end);
    varS.init.tFS = timeOpt(2:2:end);
    
    % Initialize time vector
    timeVec = 0; %[linspace(0, timeOpt(1), vS.params.N)];
    
    % For each gait phase
    for i=1:length(vars.optims.timeOpt)
        
        % Check if on first gait phase
        if i == 1
            
            % Build time vector
            timeVec = [linspace(timeVec(end),...
                (timeVec(end) + timeOpt(i)), varS.params.N)];
            
        else
            
            % Append to time vector
            timeVec = [timeVec, linspace((timeVec(end) +...
                timeOpt(i)/varS.params.N),...
                (timeVec(end) + timeOpt(i)), varS.params.N)];
            
        end
        
    end
    
    % Store time vector
    varS.optims.timeElem = timeVec;
    
    % Initialize time vector for full state
    varS.optims.timeFull = 0;
    
    % Check if collocation method was used
    if strcmp(varS.params.method, 'Collocation')
       
        % For each gait phase
        for i = 1:length(varS.optims.timeOpt)

            % For each shooting element
            for j = 1:varS.params.N
                
                % Check if first gait phase and first shooting element
                if (i == 1) && (j == 1)
                
                    % Build time vector
                    varS.optims.timeFull = [(varS.optims.timeFull(end) +...
                        dynS.colPoints(2:end)*...
                        (varS.optims.timeOpt(i)/varS.params.N))];
                    
                else
                    
                    % Append to time vector
                    varS.optims.timeFull = [varS.optims.timeFull,...
                        (varS.optims.timeFull(end) +...
                        dynS.colPoints(2:end)*...
                        (varS.optims.timeOpt(i)/varS.params.N))];
                    
                end

            end

        end
        
    elseif strcmp(varS.params.method, 'Multishooting')
        
        % Initialize time vector
        timeVec = 0;
    
        % For each gait phase
        for i=1:length(timeOpt)
            
            % Check if on first gait phase 
            if i == 1
                
                % Build time vector
                timeVec = [linspace(timeVec(end),...
                    (timeVec(end) + timeOpt(i)),...
                    varS.params.N*varS.params.M)];
                
            else
                
                % Append to time vector
                timeVec = [timeVec, linspace((timeVec(end) +...
                    timeOpt(i)/(varS.params.N*varS.params.M)),...
                    (timeVec(end) + timeOpt(i)), varS.params.N*varS.params.M)];
                
            end

        end
        
        % Store time vector for full state
        varS.optims.timeFull = timeVec;
        
    end
    
    % Store normalized time vectors
    varS.optims.timeNorm = varS.optims.timeElem/max(varS.optims.timeElem);
    varS.optims.timeNormFull = varS.optims.timeFull/max(varS.optims.timeFull);

    % Grab nominal leg length variable from wOpt
    len0Opt = reshape(wOpt(varS.inds.leg(:)),[],1);
    varS.optims.len0Opt = len0Opt;
    
    %% %%%%%%%%%%%%%%%%%%%%% Full System Creation %%%%%%%%%%%%%%%%%%%%%% %%

    % Update Optimized Foot Data for Future Use
    varS = footUpdate(varS);

    % Run Post Optimization RK4 Integration if Using Multishooting
    if strcmp(varS.params.method, 'Multishooting')

        varS = intRK4(varS);

    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% Error Analysis %%%%%%%%%%%%%%%%%%%%%%%%% %%
    
    % Check if objective cost is below  threshold or maximum runs achieved
    if (varS.optims.obj < 1e-04) || (varS.params.runs >= 4)
        
        % Calculate GRF Data From Optimized COM Data
        varS.params.size = 'Base';
        [varS, datS] = calcGRF(varS, datS);
        
        varS.params.size = 'Full';
        [varS, datS] = calcGRF(varS, datS);

        % Calculate Quantified Measures for Goodness of Fit Analysis
        varS = quantOptim(varS, datS);
        
    end
    
end