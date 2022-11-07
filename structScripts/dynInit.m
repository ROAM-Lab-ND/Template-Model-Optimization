%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dynInit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 09 December 2021
% Last Updated: 12 April 2022

% This function is used to initialize the dynamics (dynS) structure that is
% used for template model optimization.

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   dynS - Structure used for dynamics storage
%   problem - Structure used to pass necessary information to NLP solver

% OUTPUTS:
%   dynS - Structure used for dynamics storage
%   problem - Structure used to pass necessary information to NLP solver

%% %%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%% %%

function [dynS, problem] = dynInit(varS, dynS, problem)

    % Import CasADi toolbox
    import casadi.*

    % Create CoM kinematic variables
    pCOM = SX.sym('pCOM',2); % COM position variables
    vCOM = SX.sym('vCOM',2); % COM velocity variables

    % Create leg stiffness variables
    kLag = SX.sym('kLag'); % Back leg spring stiffness
    kLead = SX.sym('kLead'); % Front/single support leg spring stiffness

    % Create time variable
    tHum = SX.sym('tHum'); % Limit of integration variable

    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model,'VPP')

        % Create trunk kinematic variables
        angTrunk = SX.sym('angTrunk',1); % Trunk angular position
        rotTrunk = SX.sym('rotTrunk',1); % Trunk angular velocity

        % Check if using varying or constant VP
        if strcmp(varS.params.vppType, 'Varying')

            % Create VP variables
            rVPPD = SX.sym('rVPPD'); % VPP distance from COM along trunk
            rVPPS = SX.sym('rVPPS');

        end

        % Stack state variables
        state = [pCOM; angTrunk; vCOM; rotTrunk; kLag; kLead];

    elseif strcmp(varS.params.model,'BSLIP')

        % Stack state variables
        state = [pCOM; vCOM; kLag; kLead];

    end

    %%%%%%%%%%%%%%%%%%%%%% Double Support Symbolics %%%%%%%%%%%%%%%%%%%%%%%

    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model,'VPP')

        % Create state variables
        stateD = SX.sym('stateD', 8);
        
    elseif strcmp(varS.params.model,'BSLIP')

        % Create state variables
        stateD = SX.sym('stateD', 6);

    end

    % Create foot position variable
    footDLag = SX.sym('footDLag');
    
    % Create touchdown angle variable
    thetaD = SX.sym('thetaD');
    
    % Create leg stiffness ROC variables
    kDLagDot = SX.sym('kDLagDot');
    kDLeadDot = SX.sym('kDLeadDot');
    
    % Create limit of integration variable (phase duration)
    tfD = SX.sym('tfD');
    
    % Create nominal leg length variable
    legD = SX.sym('legD');

    % Check if using VPP or BSLIP model and if using varying or constant VP
    if isfield(varS.params, 'vppType') &&...
            strcmp(varS.params.vppType, 'Varying')

        % Stack non-state variables
        pD = [stateD; thetaD; kDLagDot; kDLeadDot; footDLag;...
            tfD; tHum; legD; rVPPD];

    else

        % Stack non-state variables
        pD = [stateD; thetaD; kDLagDot; kDLeadDot; footDLag;...
            tfD; tHum; legD];

    end
    
    % Create input variable
    u = SX.sym('u');

    %%%%%%%%%%%%%%%%%%%%%% Single Support Symbolics %%%%%%%%%%%%%%%%%%%%%%%
    
    % Create foot position variable
    footSStance = SX.sym('footSStance');
    
    % Create leg stiffness ROC variables
    kSLagDot = SX.sym('kSLagDot');
    kSLeadDot = SX.sym('kSLeadDot');
    
    % Create limit of integration variable (phase duration)
    tfS = SX.sym('tfS');
    
    % Create nominal leg length variable
    legS = SX.sym('legS');

    % Check if using VPP or BSLIP model and if using varying or constant VP
    if isfield(varS.params, 'vppType') &&...
            strcmp(varS.params.vppType, 'Varying')

        % Stack non-state variables
        pS = [kSLagDot; kSLeadDot; footSStance; tfS; tHum; legS; rVPPS];

    else

        % Stack non-state variables
        pS = [kSLagDot; kSLeadDot; footSStance; tfS; tHum; legS];

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if using polyfit or Fourier fit
    if varS.params.fitType == 1

        % Store subject CoM position
        humDatVert = fitFourier(varS.fit.humPosVertF, tHum, 'Vert');
        humDatHor = fitFourier(varS.fit.humPosHorF, tHum, 'Hor');

    else

        % Store subject CoM position
        humDatVert = polyval(varS.fit.humPosVert,tHum);
        humDatHor = polyval(varS.fit.humPosHor,tHum);

    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    %%%%%%%%%%%%%%%%%%%%%%% Double Support Dynamics %%%%%%%%%%%%%%%%%%%%%%%

    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model,'VPP')

        % Calculate initial hip position
        xHip0D = stateD(1) - varS.params.rH*sin(stateD(3));
        zHip0D = stateD(2) - varS.params.rH*cos(stateD(3));
        hip0D = [xHip0D; zHip0D];

        % Calculate current hip position
        xHipD = state(1) - varS.params.rH*sin(state(3));
        zHipD = state(2) - varS.params.rH*cos(state(3));
        hipD = [xHipD; zHipD];

    elseif strcmp(varS.params.model,'BSLIP')

        % Calculate initial hip position
        xHip0D = stateD(1);
        zHip0D = stateD(2);
        hip0D = [xHip0D; zHip0D];

        % Calculate current hip position
        xHipD = state(1);
        zHipD = state(2);
        hipD = [xHipD; zHipD];

    end

    % Calculate stride length
    stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag;

    % Calculate current leg lengths
    lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2);
    lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2);

    % Calculate spring force along legs
    forceDSpringLag = 1000*state(end-1).*(legD - lenDLag);
    forceDSpringLead = 1000*state(end).*(legD - lenDLead);
    
    % Calculate leg orientation terms
    sinDLegLag = (hipD(1) - footDLag)./lenDLag;
    cosDLegLag = hipD(2)./lenDLag;
    sinDLegLead = -(stride - (hipD(1) - footDLag))./lenDLead;
    cosDLegLead = hipD(2)./lenDLead;
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')

        % Calculate leg-trunk orientation terms
        sinDPsiLag = sinDLegLag*cos(state(3)) - cosDLegLag*sin(state(3));
        cosDPsiLag = cosDLegLag*cos(state(3)) + sinDLegLag*sin(state(3));
        sinDPsiLead = sinDLegLead*cos(state(3)) - cosDLegLead*sin(state(3));
        cosDPsiLead = cosDLegLead*cos(state(3)) + sinDLegLead*sin(state(3));

        % Check if using varying or constant VP
        if strcmp(varS.params.vppType, 'Varying')

            % Calculate hip-leg to VP orientation terms
            tanBetaLag = (varS.params.rH + rVPPD)*sinDPsiLag/...
                (lenDLag + (varS.params.rH + rVPPD)*cosDPsiLag);
            tanBetaLead = (varS.params.rH + rVPPD)*sinDPsiLead/...
                (lenDLead + (varS.params.rH + rVPPD)*cosDPsiLead);

        else

            % Calculate hip-leg to VP orientation terms
            tanBetaLag = (varS.params.rH + varS.params.rVPPD)*sinDPsiLag/...
                (lenDLag + (varS.params.rH + varS.params.rVPPD)*cosDPsiLag);
            tanBetaLead = (varS.params.rH + varS.params.rVPPD)*sinDPsiLead/...
                (lenDLead + (varS.params.rH + varS.params.rVPPD)*cosDPsiLead);

        end

        % Calculate required torque to direct GRF to VP
        tauLead = forceDSpringLead*lenDLead*tanBetaLead;
        tauLag = forceDSpringLag*lenDLag*tanBetaLag;

        % Calculate overall resulting GRF terms
        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u) -...
            cosDLegLag.*tauLag./lenDLag;
        forceXLead = sinDLegLead.*forceDSpringLead -...
            cosDLegLead.*tauLead./lenDLead;

        forceZLag = cosDLegLag.*(forceDSpringLag + 1000*u) +...
            sinDLegLag.*tauLag./lenDLag;
        forceZLead = cosDLegLead.*forceDSpringLead +...
            sinDLegLead.*tauLead./lenDLead;

        % Build state dynamics
        dynS.stateDotD = tfD.*[state(4); state(5); state(6);...
            (1./varS.params.m).*(forceXLag + forceXLead);...
            (1./varS.params.m).*(forceZLead + forceZLag)-...
            varS.params.g;...
            (1./varS.params.J).*(tauLag + tauLead -...
            (forceXLag + forceXLead)*varS.params.rH*cos(state(3)) +...
            (forceZLag + forceZLead)*varS.params.rH*sin(state(3)));...
            kDLagDot; kDLeadDot];

    elseif strcmp(varS.params.model,'BSLIP')

        % Calculate overall resulting GRF terms
        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u);
        forceXLead = sinDLegLead.*forceDSpringLead;

        forceZLag = cosDLegLag.*(forceDSpringLag + 1000*u);
        forceZLead = cosDLegLead.*forceDSpringLead;

        % Build state dynamics
        dynS.stateDotD = tfD.*[state(3); state(4);...
            (1./varS.params.m).*(forceXLag + forceXLead);...
            (1/varS.params.m).*(forceZLag + forceZLead)-varS.params.g;...
            kDLagDot; kDLeadDot];

    end

    % Calculate objective cost
    dynS.costD = varS.params.Q(1,1)*(state(1) - humDatHor)^2 +...
        varS.params.Q(2,2)*(state(2) - humDatVert)^2;
    dynS.resVertD = abs(varS.params.Q(2,2)*(state(2) - humDatVert));
    %     dynStruct.costD = 0;

    % Generate CasADi function object
    dynS.fD = Function('fD', {state, pD, u},...
        {dynS.stateDotD, dynS.costD, dynS.resVertD});

    %%%%%%%%%%%%%%%%%%%%%%% Single Support Dynamics %%%%%%%%%%%%%%%%%%%%%%%

    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model,'VPP')

        % Calculate current hip position
        xHipS = state(1) - varS.params.rH*sin(state(3));
        zHipS = state(2) - varS.params.rH*cos(state(3));
        hipS = [xHipS; zHipS];

    elseif strcmp(varS.params.model,'BSLIP')

        % Calculate current hip position
        xHipS = state(1);
        zHipS = state(2);
        hipS = [xHipS; zHipS];

    end

    % Calculate current leg length
    lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2);

    % Calculate spring force along leg
    forceSSpring = 1000*state(end).*(legS - lenStance);

    % Calculate leg orientation terms
    sinSLeg = (hipS(1) - footSStance)./lenStance; % X-Axis Weight
    cosSLeg = hipS(2)./lenStance; % Y-Axis Weight

    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model,'VPP')

        % Calculate leg-trunk orientation terms
        sinSPsi = sinSLeg*cos(state(3)) - cosSLeg*sin(state(3));
        cosSPsi = cosSLeg*cos(state(3)) + sinSLeg*sin(state(3));

        % Check if using varying or constant VP
        if strcmp(varS.params.vppType, 'Varying')

            % Calculate hip-leg to VP orientation term
            tanBetaS = (varS.params.rH + rVPPS)*sinSPsi/...
                (lenStance + (varS.params.rH + rVPPS)*cosSPsi);

        else

            % Calculate hip-leg to VP orientation term
            tanBetaS = (varS.params.rH + varS.params.rVPPS)*sinSPsi/...
                (lenStance + (varS.params.rH + varS.params.rVPPS)*cosSPsi);

        end

        % Calculate required torque to direct GRF to VP
        tauS = forceSSpring*lenStance*tanBetaS;

        % Calculate overall resulting GRF terms
        forceX = sinSLeg.*forceSSpring - cosSLeg.*tauS./lenStance;
        forceZ = cosSLeg.*forceSSpring + sinSLeg.*tauS./lenStance;

        % Build state dynamics
        dynS.stateDotS = tfS.*[state(4);state(5);state(6);...
            (1./varS.params.m).*forceX;...
            (1./varS.params.m).*forceZ - varS.params.g;
            (1./varS.params.J).*(tauS -...
            forceX*varS.params.rH*cos(state(3)) +...
            forceZ*varS.params.rH*sin(state(3)));
            kSLagDot; kSLeadDot];

    elseif strcmp(varS.params.model,'BSLIP')

        % Calculate overall resulting GRF terms
        forceX = sinSLeg.*forceSSpring;
        forceZ = cosSLeg.*forceSSpring;

        % Build state dynamics
        dynS.stateDotS = tfS*[state(3); state(4);...
            (1./varS.params.m).*(forceX);...
            (1./varS.params.m).*forceZ-varS.params.g;...
            kSLagDot; kSLeadDot];

    end

    % Calculate objective cost
    dynS.costS = varS.params.Q(1,1)*(state(1) - humDatHor)^2 +...
        varS.params.Q(2,2)*(state(2) - humDatVert)^2;
    dynS.resVertS = abs(varS.params.Q(2,2)*(state(2) - humDatVert));
    %     dynStruct.costS = 0;

    % Generate CasADi function object
    dynS.fS = Function('fS', {state, pS},...
        {dynS.stateDotS, dynS.costS, dynS.resVertS});

    %% %%%%%%%%%%%%%%%%%%%%%% Multishooting (RK4) %%%%%%%%%%%%%%%%%%%%%% %%

    % Check if using multishooting or collocation
    if strcmp(varS.params.method, 'Multishooting')

        % Create state variables
        STATE0 = MX.sym('STATE0',length(state));
        
        % Grab time step
        DT = varS.params.dt;

        %%%%%%%%%%%%%%%%%%%%%%% Double Support RK4 %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Create state variable
        STATED = STATE0;
        
        % Create non-state variables
        PD = MX.sym('PD',length(pD));
        
        % Create input variable
        U = MX.sym('U');
        
        % Create cost variables
        QD = 0;
        RD = 0;

        % For each finite elements
        for i=1:varS.params.M

            % Call CasADi function object for double support dynamics
            [k1D, k1qD, k1rD] = dynS.fD(STATED,PD,U);
            [k2D, k2qD, k2rD] = dynS.fD(STATED + (DT*k1D)/2,PD,U);
            [k3D, k3qD, k3rD] = dynS.fD(STATED + (DT*k2D)/2,PD,U);
            [k4D, k4qD, k4rD] = dynS.fD(STATED + (DT*k3D),PD,U);

            % Propogate state and cost based on RK4 method
            STATED = STATED + (DT/6)*(k1D + 2*k2D + 2*k3D + k4D);
            QD = QD + (DT/6)*(k1qD + 2*k2qD + 2*k3qD + k4qD);
            RD = RD + (DT/6)*(k1rD + 2*k2rD + 2*k3rD + k4rD);

        end

        % Generate CasADi function object
        dynS.FD = Function('FD', {STATE0, PD, U}, {STATED, QD},...
            {'state0', 'p', 'u'}, {'stated', 'qd', 'rd'});

        %%%%%%%%%%%%%%%%%%%%%%% Single Support RK4 %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Create state variables
        STATES = STATE0;
        
        % Create non-state variables
        PS = MX.sym('PS',length(pS));
        
        % Create cost variables
        QS = 0;
        RS = 0;

        % For each finite element
        for i=1:varS.params.M

            % Call CasADi function object for single support dynamics
            [k1S, k1qS, k1rS] = dynS.fS(STATES,PS);
            [k2S, k2qS, k2rS] = dynS.fS(STATES + (DT*k1S)/2,PS);
            [k3S, k3qS, k3rS] = dynS.fS(STATES + (DT*k2S)/2,PS);
            [k4S, k4qS, k4rS] = dynS.fS(STATES + (DT*k3S),PS);

            % Propogate state and cost based on RK4 method
            STATES = STATES + (DT/6)*(k1S + 2*k2S + 2*k3S + k4S);
            QS = QS + (DT/6)*(k1qS + 2*k2qS + 2*k3qS + k4qS);
            RS = RS + (DT/6)*(k1rS + 2*k2rS + 2*k3rS + k4rS);

        end

        % Generate CasADi function object
        dynS.FS = Function('FS', {STATE0, PS}, {STATES, QS},...
            {'state0', 'p'}, {'states', 'qs', 'rs'});

    end

end