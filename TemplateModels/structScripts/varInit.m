%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% varInit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 9 December 2021
% Last Updated: 9 December 2021

% This function is used to initialize the variable (vS) and NLP (problem)
% structures that are used for template model optimization. Initial bounds, 
% parameters, and variables/constraint storage is initialized by calling 
% this function.

%% %%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, problem] = varInit(vS, dynS, dataS, params)
    
    % Initialize problem structure for storing variables, constraints, and cost
    problem.vars = {};
    problem.varsInit = [];
    problem.varsLB = [];
    problem.varsUB = [];

    problem.cost = 0;
    problem.resVert = 0;

    problem.constraints = {};
    problem.constraintsLB = [];
    problem.constraintsUB = [];

    % Initialize variable structure for storing variables that can be
    % incorporated into constraints
    vS.inds.state = [];
    vS.inds.stateFull = [];
    vS.inds.foot = [];
    vS.inds.theta = [];
    vS.inds.u = [];
    vS.inds.kLagDot = [];
    vS.inds.kLeadDot = [];
    vS.inds.time = [];
    vS.inds.gen = [];
    
    if strcmp(params{19}, 'Varying')
        
        vS.inds.rVPP = [];
        
    end
    
    vS.params.model = params{15};
    vS.params.method = params{16};
    vS.params.size = 'Base';
    vS.params.source = params{17};
    vS.params.springType = params{18};
    vS.params.fitType = params{20};
    
    vS.params.len0 = dataS.subj.len0;
    vS.params.m = dataS.subj.weight;
    vS.params.g = dataS.subj.g;
    
    if strcmp(vS.params.model, 'VPP')
        
        vS.params.vppType = params{19};
        
        vS.params.rH = params{3};
        
        if strcmp(vS.params.vppType, 'Constant')
            
            vS.params.rVPPD = params{4}(1);
            vS.params.rVPPS = params{4}(2);
        
        end
        
        vS.params.gamma = params{2};
        vS.params.J = params{5};
        
    end
    
    vS.params.xInit = dataS.subj.treadmill;
    vS.params.vDiff = params{6};
    vS.params.N = params{8};
    vS.params.M = params{9};
    vS.params.dt = 1/(params{8}*params{9});
    vS.params.steps = params{1};
    vS.params.inputCountD = params{11};
    vS.params.inputCountS = params{12};
    vS.params.shootCount = params{13};
    vS.params.timeTrack = params{14};
    vS.params.Q = params{21};
    vS.params.residualVert = 0;
    
    if strcmp(vS.params.springType, 'Constant')
        
        vS.params.stiffVaryMax = 0;
        vS.params.stiffPhaseMax = 0;
        
    else
        
        vS.params.stiffVaryMax = 10;
        vS.params.stiffPhaseMax = 0.05;
        
    end
    
    switch params{15}
        
        case 'VPP'
            
            % Set Bounds on Variables
            vS.lb.state = [0; dataS.subj.len0*cos(pi/6); -5*pi/180;...
                dataS.subj.treadmill-0.5; -0.5; -100*pi/180;...
                5; 5]; % State Lower Bounds
            vS.ub.state = [2*params{1}; 1.5*(dataS.subj.len0+params{3});...
                10*pi/180; dataS.subj.treadmill+0.5; 0.5; 100*pi/180;...
                50; 50]; % State Upper Bounds
            
            % Set Bounds on Variables
            vS.lb.stateDS = [0; dataS.subj.len0*cos(pi/6); -5*pi/180;...
                dataS.subj.treadmill-0.5; -0.5; -100*pi/180;...
                5; 5]; % State Lower Bounds
            vS.ub.stateDS = [2*params{1}; 1.5*(dataS.subj.len0+params{3});...
                10*pi/180; dataS.subj.treadmill+0.5; 0.5; 100*pi/180;...
                50; 50]; % State Upper Bounds
            
        case 'BSLIP'
            
            % Set Bounds on Variables
            vS.lb.state = [0; dataS.subj.len0*cos(pi/6);...
                dataS.subj.treadmill-0.5; -0.5; 5; 5]; % State Lower Bounds
            vS.ub.state = [2*params{1}; 2*dataS.subj.len0;...
                dataS.subj.treadmill+0.5; 0.5; 50; 50]; % State Upper Bounds
            
            % Set Bounds on Variables
            vS.lb.stateDS = [0; dataS.subj.len0*cos(pi/6);...
                dataS.subj.treadmill-0.5; -0.5; 5; 5]; % State Lower Bounds
            vS.ub.stateDS = [2*params{1}; 2*dataS.subj.len0;...
                dataS.subj.treadmill+0.5; 0.5; 50; 50]; % State Upper Bounds
            
        otherwise
            
            errMsg = 'Inputted Template Model Not Yet Implemented';
            error(errMsg);
            
    end
    
    % Limit of Integration Lower Bound
    if length(dataS.data.timeDS1)>1
        
        vS.lb.tfD = min(min(dataS.data.timeDS1),...
            min(dataS.data.timeDS2)) - 0.02; 
        % Limit of Integration Upper Bound Double Support
        vS.ub.tfD = max(max(dataS.data.timeDS1),...
            max(dataS.data.timeDS2)) + 0.02; 
        % Limit of Integration Lower Bound
        vS.lb.tfS = min(min(dataS.data.timeSS1),...
            min(dataS.data.timeSS2)) - 0.02; 
        % Limit of Integration Upper Bound Single Support
        vS.ub.tfS = max(max(dataS.data.timeSS1),...
            max(dataS.data.timeSS2)) + 0.02;
        
%         vS.lb.tfD = min(dataS.data.timeDS1(dataS.inds.start),...
%             dataS.data.timeDS2(dataS.inds.start)) - 0.02; 
%         % Limit of Integration Upper Bound Double Support
%         vS.ub.tfD = max(dataS.data.timeDS1(dataS.inds.start),...
%             dataS.data.timeDS2(dataS.inds.start)) + 0.02; 
%         % Limit of Integration Lower Bound
%         vS.lb.tfS = min(dataS.data.timeSS1(dataS.inds.start),...
%             dataS.data.timeSS2(dataS.inds.start)) - 0.02; 
%         % Limit of Integration Upper Bound Single Support
%         vS.ub.tfS = max(dataS.data.timeSS1(dataS.inds.start),...
%             dataS.data.timeSS2(dataS.inds.start)) + 0.02;
        
    else
        
        vS.lb.tfD = min(dataS.data.timeDS1, dataS.data.timeDS2) - 0.02; 
        % Limit of Integration Upper Bound Double Support
        vS.ub.tfD = max(dataS.data.timeDS1, dataS.data.timeDS2) + 0.02; 
        % Limit of Integration Lower Bound
        vS.lb.tfS = min(dataS.data.timeSS1, dataS.data.timeSS2) - 0.02;
        % Limit of Integration Upper Bound Single Support
        vS.ub.tfS = max(dataS.data.timeSS1, dataS.data.timeSS2) + 0.02;
        
    end
    
    vS.lb.foot = -1; % Foot X-Position Lower Bound
    vS.ub.foot = 5; % Foot X-Position Upper Bound
    vS.lb.theta = 5*pi/180; % Touchdown Angle Lower Bound
    vS.ub.theta = 45*pi/180; % Touchdown Angle Upper Bound
    
    if strcmp(vS.params.springType, 'Constant')
        
        vS.lb.kDotDS = 0; % Leg Stiffness Differential Lower Bound
        vS.ub.kDotDS = 0; % Leg Stiffness Differential Upper Bound
        vS.lb.kDotSS = 0; % Leg Stiffness Differential Lower Bound
        vS.ub.kDotSS = 0; % Leg Stiffness Differential Upper Bound
        
    else
        
        vS.lb.kDotDS = -500; % Leg Stiffness Differential Lower Bound
        vS.ub.kDotDS = 500; % Leg Stiffness Differential Upper Bound
        vS.lb.kDotSS = -100; % Leg Stiffness Differential Lower Bound
        vS.ub.kDotSS = 100; % Leg Stiffness Differential Upper Bound
        
    end
    
    vS.lb.u = 0; % Control Input Lower Bound
    vS.ub.u = 0; % Control Input Upper Bound
    vS.lb.eps = -params{7};
    vS.ub.eps = params{7};
    vS.lb.leg = dataS.subj.len0 - 0.4;
    vS.ub.leg = dataS.subj.len0 + 0.4;
    
    if strcmp(params{19}, 'Varying')
        
        vS.lb.rVPP = 0;
        vS.ub.rVPP = 0.5;
        
    end
    
    % Fit Coefficients of Vertical COM Trajectory wrt Time
    vS.fit.humPosVert = dataS.fit.coefPosVert;
    vS.fit.humPosVertF = dataS.fit.coefPosVertF;

    % Fit Coefficients of Horizontal COM Trajectory wrt Time
    vS.fit.humPosHor = dataS.fit.coefPosHor;
    vS.fit.humPosHorF = dataS.fit.coefPosHorF;
        
    % Polyfit Coefficients of Vertical COM Velocity wrt Time    
    vS.fit.humVelVert = dataS.fit.coefVelVert;
    
    % Polyfit Coefficients of Horizontal COM Velocity wrt Time
    vS.fit.humVelHor = dataS.fit.coefVelHor;
    
    % Initialize Initial Guesses for Variables
    if params{10}
        
        vS.params.runs = 1;
       
        vS.fit.tHum = 0;
        
        for i = 1:length(dataS.data.timeDS1)
            
            if i == 1
                
                vS.fit.tHum = [linspace(vS.fit.tHum(end),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   (dataS.data.timeSS1(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + (dataS.data.timeDS2(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                   (dataS.data.timeSS2(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                   dataS.data.timeSS2(i)), params{8})];
               
            else
                
                vS.fit.tHum = [vS.fit.tHum, linspace(vS.fit.tHum(end),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   (dataS.data.timeSS1(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + (dataS.data.timeDS2(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i)), params{8}),...
                   ...
                   linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                   (dataS.data.timeSS2(i)/params{8})),...
                   (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
                   dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                   dataS.data.timeSS2(i)), params{8})];          
                
            end
           
%             vS.fit.tHum = [vS.fit.tHum,...
%                 linspace((vS.fit.tHum(end) + (dataS.data.timeDS1(i)/params{8})),...
%                (vS.fit.tHum(end) + dataS.data.timeDS1(i)), params{8}),...
%                ...
%                linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                (dataS.data.timeSS1(i)/params{8})),...
%                (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                dataS.data.timeSS1(i)), params{8}),...
%                ...
%                linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                dataS.data.timeSS1(i) + (dataS.data.timeDS2(i)/params{8})),...
%                (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                dataS.data.timeSS1(i) + dataS.data.timeDS2(i)), params{8}),...
%                ...
%                linspace((vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
%                (dataS.data.timeSS2(i)/params{8})),...
%                (vS.fit.tHum(end) + dataS.data.timeDS1(i) +...
%                dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
%                dataS.data.timeSS2(i)), params{8})];
           
        end
               
%            vS.fit.tHum = [0, linspace((dataS.data.timeDS1/params{8}),...
%                dataS.data.timeDS1, params{8}),...
%                linspace(dataS.data.timeDS1 + (dataS.data.timeSS1/params{8}),...
%                dataS.data.timeDS1 + dataS.data.timeSS1, params{8}),...
%                linspace(dataS.data.timeDS1 + dataS.data.timeSS1 +...
%                (dataS.data.timeDS2/params{8}), dataS.data.timeDS1 +...
%                dataS.data.timeSS1 + dataS.data.timeDS2, params{8}),...
%                linspace(dataS.data.timeDS1 + dataS.data.timeSS1 +...
%                dataS.data.timeDS2 + (dataS.data.timeSS2/params{8}),...
%                dataS.data.timeGait, params{8})];

       vS.fit.tHumFull = 0;
       
       if strcmp(vS.params.method, 'Collocation')
       
           for i = 1:length(dataS.data.timeFull)

               for j = 1:params{8}

                   if (j == 1) && (i == 1)
                       
                       vS.fit.tHumFull = [(vS.fit.tHumFull(end) +...
                           dynS.colPoints(2:end)*...
                           (dataS.data.timeFull(i)/params{8}))];
                       
                   else
                       
                       vS.fit.tHumFull = [vS.fit.tHumFull,...
                           (vS.fit.tHumFull(end) +...
                           dynS.colPoints(2:end)*...
                           (dataS.data.timeFull(i)/params{8}))];
                       
                   end

               end

           end
           
       elseif strcmp(vS.params.method, 'Multishooting')
           
           for i = 1:length(dataS.data.timeDS1)
               
               if i == 1
           
                    vS.fit.tHumFull = [linspace(vS.fit.tHumFull(end),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i)),...
                       params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) +...
                       (dataS.data.timeSS1(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i)), params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) + dataS.data.timeSS1(i) +...
                       (dataS.data.timeDS2(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i) + dataS.data.timeDS2(i)),...
                       params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) + dataS.data.timeSS1(i) +...
                       dataS.data.timeDS2(i) +...
                       (dataS.data.timeSS2(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                       dataS.data.timeSS2(i)), params{9}*params{8})];
                   
               else
                   
                   
                    vS.fit.tHumFull = [vS.fit.tHumFull, linspace(vS.fit.tHumFull(end),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i)),...
                       params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) +...
                       (dataS.data.timeSS1(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i)), params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) + dataS.data.timeSS1(i) +...
                       (dataS.data.timeDS2(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i) + dataS.data.timeDS2(i)),...
                       params{9}*params{8}),...
                       ...
                       linspace((vS.fit.tHumFull(end) +...
                       dataS.data.timeDS1(i) + dataS.data.timeSS1(i) +...
                       dataS.data.timeDS2(i) +...
                       (dataS.data.timeSS2(i)/(params{9}*params{8}))),...
                       (vS.fit.tHumFull(end) + dataS.data.timeDS1(i) +...
                       dataS.data.timeSS1(i) + dataS.data.timeDS2(i) +...
                       dataS.data.timeSS2(i)), params{9}*params{8})];
                   
               end
                   

            end
           
%            vS.fit.tHumFull = [0, linspace((dataS.data.timeDS1/(params{9}*params{8})), dataS.data.timeDS1,...
%                params{9}*params{8}),...
%                linspace(dataS.data.timeDS1 +...
%                (dataS.data.timeSS1/(params{9}*params{8})),...
%                dataS.data.timeDS1 + dataS.data.timeSS1, params{9}*params{8}),...
%                linspace(dataS.data.timeDS1 + dataS.data.timeSS1 +...
%                (dataS.data.timeDS2/(params{9}*params{8})), dataS.data.timeDS1 +...
%                dataS.data.timeSS1 + dataS.data.timeDS2, params{9}*params{8}),...
%                linspace(dataS.data.timeDS1 + dataS.data.timeSS1 +...
%                dataS.data.timeDS2 + (dataS.data.timeSS2/(params{9}*params{8})),...
%                dataS.data.timeGait, params{9}*params{8})];
           
       end
       
       vS.fit.posCOMSegHor = polyval(vS.fit.humPosHor,vS.fit.tHum);
       vS.fit.posCOMSegVert = polyval(vS.fit.humPosVert,vS.fit.tHum);
       vS.fit.velCOMSegHor = polyval(vS.fit.humVelHor,vS.fit.tHum);
       vS.fit.velCOMSegVert = polyval(vS.fit.humVelVert,vS.fit.tHum);
       
       vS.fit.posCOMSegHorF = fitFourier(vS.fit.humPosHorF,...
           vS.fit.tHum, 'Hor');
       vS.fit.posCOMSegVertF = fitFourier(vS.fit.humPosVertF,...
           vS.fit.tHum, 'Vert');
       vS.fit.velCOMSegHorF = fitFourierDer(vS.fit.humPosHorF,...
           vS.fit.tHum, 'Hor');
       vS.fit.velCOMSegVertF = fitFourierDer(vS.fit.humPosVertF,...
           vS.fit.tHum, 'Vert');
       
       switch params{15}
           
           case 'VPP'
       
               vS.init.state = zeros(2*params{1}*params{8},8);
               vS.init.state(:,3) = 10*pi/180;
               
               if params{20} > 1
                   
                   vS.init.state(:,4) = vS.fit.velCOMSegHor;
                   vS.init.state(:,5) = vS.fit.velCOMSegVert;
                   
               else
                   
                   vS.init.state(:,4) = vS.fit.velCOMSegHorF;
                   vS.init.state(:,5) = vS.fit.velCOMSegVertF;
                   
               end
               
               if strcmp(params{19}, 'Varying')
                   
                   vS.init.rVPPS = zeros(params{1},1);
                   vS.init.rVPPD = 0.3*ones(params{1},1);
                   
               end
               
           case 'BSLIP'
               
               vS.init.state = zeros(2*params{1}*params{8},6);
               
               if params{20} > 1
                   
                   vS.init.state(:,3) = vS.fit.velCOMSegHor;
                   vS.init.state(:,4) = vS.fit.velCOMSegVert;
                   
               else
                   
                   vS.init.state(:,3) = vS.fit.velCOMSegHorF;
                   vS.init.state(:,4) = vS.fit.velCOMSegVertF;
                   
               end
               
               
           otherwise
               
               disp('Non-implemented model chosen')
               
       end
       
       if params{20} > 1
           
           vS.init.state(:,1) = vS.fit.posCOMSegHor;
           vS.init.state(:,2) = vS.fit.posCOMSegVert;
           
       else
           
           vS.init.state(:,1) = vS.fit.posCOMSegHorF;
           vS.init.state(:,2) = vS.fit.posCOMSegVertF;
           
       end
       
       vS.init.state(:,end-1) = 20;
       vS.init.state(:,end) = 20;
       
       vS.init.kLagDot = zeros(2*params{1}*params{8},1);
       vS.init.kLeadDot = zeros(2*params{1}*params{8},1);
       
       vS.init.leg = vS.params.len0;
                  
       tempTfD = [dataS.data.timeDS1'; dataS.data.timeDS2'];
       tempTfS = [dataS.data.timeSS1'; dataS.data.timeSS2'];
       vS.init.tFD = tempTfD(:);
       vS.init.tFS = tempTfS(:);
       
       vS.init.theta = 20*pi/180*ones((params{1} + 1), 1);
       
       vS.init.footD = 0.5*(0:params{1});
       vS.init.footS = (2/3) + (2/3)*(0:params{1}-1);
       
       vS.init.u = zeros(params{1}*params{8},1);       
        
    end
    
    for i=1:params{1}

        % Double Support Phase
        % Add limit of integration variable for double support phase
        [problem, vS.var.(['tfD' num2str(i)]), vS.inds.time(end+1,:)] =...
            addVariable(problem, ['tfD' num2str(i)], 1, vS.init.tFD(i),...
            vS.lb.tfD, vS.ub.tfD);
        % Add touchdown angle variable for double support phase
        [problem, vS.var.(['thetaD' num2str(i)]),...
            vS.inds.theta(end+1,:)] = addVariable(problem,...
            ['thetaD' num2str(i)], 1, vS.init.theta(i),...
            vS.lb.theta, vS.ub.theta);
        % Add foot location variable for double support phase
        [problem, vS.var.(['footDLag' num2str(i)]),...
            vS.inds.foot(end+1,:)] = addVariable(problem,...
            ['footDLag' num2str(i)], 1, vS.init.footD(i),...
            vS.lb.foot, vS.ub.foot);

        % Single Support Phase
        % Add limit of integration variable for single support phase
        [problem, vS.var.(['tfS' num2str(i)]), vS.inds.time(end+1,:)] =...
            addVariable(problem, ['tfS' num2str(i)], 1, vS.init.tFS(i),...
            vS.lb.tfS, vS.ub.tfS);
        % Add foot location variable for single support phase
        [problem, vS.var.(['footSStance' num2str(i)]),...
            vS.inds.foot(end+1,:)] = addVariable(problem,...
            ['footSStance' num2str(i)], 1, vS.init.footS(i),...
            vS.lb.foot, vS.ub.foot);
        
        % Create Variables for Varying VPP Method
        if strcmp(vS.params.model, 'VPP')
            
            if strcmp(vS.params.vppType, 'Varying')
            
                [problem, vS.var.(['rVPPD' num2str(i)]),...
                    vS.inds.rVPP(end+1,:)] = addVariable(problem,...
                    ['rVPPD' num2str(i)], 1,...
                    vS.init.rVPPD(i), vS.lb.rVPP, vS.ub.rVPP);

                [problem, vS.var.(['rVPPS' num2str(i)]),...
                    vS.inds.rVPP(end+1,:)] = addVariable(problem,...
                    ['rVPPS' num2str(i)], 1,...
                    vS.init.rVPPS(i), vS.lb.rVPP, vS.ub.rVPP);
                
            end
            
        end        

    end

    % Add final touchdown angle variable
    [problem, vS.var.(['thetaD' num2str(i+1)]), vS.inds.theta(end+1,:)]...
        = addVariable(problem, ['thetaD' num2str(i+1)], 1,...
        vS.init.theta(end), vS.lb.theta, vS.ub.theta);
    [problem, vS.var.(['footDLag' num2str(i+1)]), vS.inds.foot(end+1,:)]...
        = addVariable(problem, ['footDLag' num2str(i+1)], 1,...
        vS.init.footD(end), vS.lb.foot, vS.ub.foot);
    [problem, vS.var.len0, vS.inds.leg(1)] = addVariable(problem,...
        'len0', 1, vS.init.leg, vS.lb.leg, vS.ub.leg);
    
end