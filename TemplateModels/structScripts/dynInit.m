function [dynStruct, problem] = dynInit(vS, dynStruct, problem)
    
    import casadi.*

    % State Symbolics
    pCOM = SX.sym('pCOM',2); % COM position variables
    vCOM = SX.sym('vCOM',2); % COM velocity variables

    kLag = SX.sym('kLag'); % Back leg spring stiffness
    kLead = SX.sym('kLead'); % Front/single support leg spring stiffness
    
    tHum = SX.sym('tHum'); % Limit of integration variable
    
    if strcmp(vS.params.model,'VPP')
        
        angTrunk = SX.sym('angTrunk',1); % Trunk angular position
        rotTrunk = SX.sym('rotTrunk',1); % Trunk angular velocity
        
        if strcmp(vS.params.vppType, 'Varying')
            
            rVPPD = SX.sym('rVPPD'); % VPP distance from COM along trunk
            rVPPS = SX.sym('rVPPS');
            
        end

        state = [pCOM; angTrunk; vCOM; rotTrunk; kLag; kLead]; % COM full state
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        state = [pCOM; vCOM; kLag; kLead]; % COM full state
        
    end    

    % Double Support Symbolics
    
    if strcmp(vS.params.model,'VPP')
        
        stateD = SX.sym('stateD', 8); % Double support state capture
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        stateD = SX.sym('stateD', 6);
        
    end
    
    footDLag = SX.sym('footDLag'); % Double support lag foot position
    thetaD = SX.sym('thetaD'); % Double support touchdown angle
    kDLagDot = SX.sym('kDLagDot'); % Double support leg stiffness derivative
    kDLeadDot = SX.sym('kDLeadDot');
    tfD = SX.sym('tfD'); % Double support limit of integration
    legD = SX.sym('legD');
    
    if isfield(vS.params, 'vppType') && strcmp(vS.params.vppType, 'Varying')
            
        pD = [stateD; thetaD; kDLagDot; kDLeadDot; footDLag;...
            tfD; tHum; legD; rVPPD];
        
    else
        
        pD = [stateD; thetaD; kDLagDot; kDLeadDot; footDLag;...
            tfD; tHum; legD];
        
    end

    % Single Stance Symbolics
    footSStance = SX.sym('footSStance'); % Single support stance foot position
    kSLagDot = SX.sym('kSLagDot'); % Single support leg stiffness
    kSLeadDot = SX.sym('kSLeadDot');
    tfS = SX.sym('tfS'); % Single support limit of integration
    legS = SX.sym('legS');
    
    if isfield(vS.params, 'vppType') && strcmp(vS.params.vppType, 'Varying')

        pS = [kSLagDot; kSLeadDot; footSStance; tfS; tHum; legS; rVPPS];
        
    else
        
        pS = [kSLagDot; kSLeadDot; footSStance; tfS; tHum; legS];
        
    end

    % Control Input Symbolic
    u = SX.sym('u'); % Force input on lag leg during double stance
    
    % Human Polyfit
    if vS.params.fitType == 1
        
        humDatVert = fitFourier(vS.fit.humPosVertF, tHum, 'Vert');
        humDatHor = fitFourier(vS.fit.humPosHorF, tHum, 'Hor');
        
    else
        
        humDatVert = polyval(vS.fit.humPosVert,tHum);
        humDatHor = polyval(vS.fit.humPosHor,tHum);
        
    end
    

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    % ##################### Double Support Dynamics ##################### %
    
    if strcmp(vS.params.model,'VPP')
        
        xHip0D = stateD(1) - vS.params.rH*sin(stateD(3));
        zHip0D = stateD(2) - vS.params.rH*cos(stateD(3));
        hip0D = [xHip0D; zHip0D];
    
        xHipD = state(1) - vS.params.rH*sin(state(3));
        zHipD = state(2) - vS.params.rH*cos(state(3));
        hipD = [xHipD; zHipD];
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        xHip0D = stateD(1);
        zHip0D = stateD(2);
        hip0D = [xHip0D; zHip0D];
    
        xHipD = state(1);
        zHipD = state(2);
        hipD = [xHipD; zHipD];
        
    end
    
    stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag; % Stride Length
    
    lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2); % Lag Leg Length
    lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2); % Lead Leg Length

    forceDSpringLag = 1000*state(end-1).*(legD - lenDLag); % Lag Leg Spring Force
    forceDSpringLead = 1000*state(end).*(legD - lenDLead); % Lead Leg Spring Force
    
    sinDLegLag = (hipD(1) - footDLag)./lenDLag; % X-Axis Weight Lag Leg
    cosDLegLag = hipD(2)./lenDLag; % Y-Axis Weight Lag Leg
    sinDLegLead = -(stride - (hipD(1) - footDLag))./lenDLead; % X-Axis Weight Lead Leg
    cosDLegLead = hipD(2)./lenDLead; % Y-Axis Weight Lead Leg
    
    if strcmp(vS.params.model, 'VPP')
    
        sinDPsiLag = sinDLegLag*cos(state(3)) - cosDLegLag*sin(state(3));
        cosDPsiLag = cosDLegLag*cos(state(3)) + sinDLegLag*sin(state(3));
        sinDPsiLead = sinDLegLead*cos(state(3)) - cosDLegLead*sin(state(3));
        cosDPsiLead = cosDLegLead*cos(state(3)) + sinDLegLead*sin(state(3));
        
        if strcmp(vS.params.vppType, 'Varying')
            
            tanBetaLag = (vS.params.rH + rVPPD)*sinDPsiLag/...
                (lenDLag + (vS.params.rH + rVPPD)*cosDPsiLag);    
            tanBetaLead = (vS.params.rH + rVPPD)*sinDPsiLead/...
                (lenDLead + (vS.params.rH + rVPPD)*cosDPsiLead);
            
        else
 
            tanBetaLag = (vS.params.rH + vS.params.rVPPD)*sinDPsiLag/...
                (lenDLag + (vS.params.rH + vS.params.rVPPD)*cosDPsiLag);    
            tanBetaLead = (vS.params.rH + vS.params.rVPPD)*sinDPsiLead/...
                (lenDLead + (vS.params.rH + vS.params.rVPPD)*cosDPsiLead);
            
        end
    
        tauLead = forceDSpringLead*lenDLead*tanBetaLead;
        tauLag = forceDSpringLag*lenDLag*tanBetaLag;

        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u) -...
            cosDLegLag.*tauLag./lenDLag;
        forceXLead = sinDLegLead.*forceDSpringLead -...
            cosDLegLead.*tauLead./lenDLead;

        forceZLag = cosDLegLag.*(forceDSpringLag + 1000*u) +...
            sinDLegLag.*tauLag./lenDLag;    
        forceZLead = cosDLegLead.*forceDSpringLead +...
            sinDLegLead.*tauLead./lenDLead;
    
        % State Dynamics
        dynStruct.stateDotD = tfD.*[state(4); state(5); state(6);...
            (1./vS.params.m).*(forceXLag + forceXLead);...
            (1./vS.params.m).*(forceZLead + forceZLag)-...
            vS.params.g;...
            (1./vS.params.J).*(tauLag + tauLead -...
            (forceXLag + forceXLead)*vS.params.rH*cos(state(3)) +...
            (forceZLag + forceZLead)*vS.params.rH*sin(state(3)));...
            kDLagDot; kDLeadDot];
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u);
        forceXLead = sinDLegLead.*forceDSpringLead;

        forceZLag = cosDLegLag.*(forceDSpringLag + 1000*u);    
        forceZLead = cosDLegLead.*forceDSpringLead;
        
        % State Dynamics
        dynStruct.stateDotD = tfD.*[state(3); state(4);...
            (1./vS.params.m).*(forceXLag + forceXLead);...
            (1/vS.params.m).*(forceZLag + forceZLead)-vS.params.g;...
            kDLagDot; kDLeadDot];
        
    end

    % Cost Function
    dynStruct.costD = vS.params.Q(1,1)*(state(1) - humDatHor)^2 +...
        vS.params.Q(2,2)*(state(2) - humDatVert)^2;
    dynStruct.resVertD = abs(vS.params.Q(2,2)*(state(2) - humDatVert));
%     dynStruct.costD = 0;

    % CasADi Function Object
    dynStruct.fD = Function('fD', {state, pD, u},...
        {dynStruct.stateDotD, dynStruct.costD, dynStruct.resVertD});

    % ##################### Single Support Dynamics ##################### %

    if strcmp(vS.params.model,'VPP')
        
        xHipS = state(1) - vS.params.rH*sin(state(3));
        zHipS = state(2) - vS.params.rH*cos(state(3));
        hipS = [xHipS; zHipS];
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        xHipS = state(1);
        zHipS = state(2);
        hipS = [xHipS; zHipS];
        
    end
    
    lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2); % Stance Leg Length
    
    forceSSpring = 1000*state(end).*(legS - lenStance); % Spring Force

    sinSLeg = (hipS(1) - footSStance)./lenStance; % X-Axis Weight
    cosSLeg = hipS(2)./lenStance; % Y-Axis Weight
    
    if strcmp(vS.params.model,'VPP')
    
        sinSPsi = sinSLeg*cos(state(3)) - cosSLeg*sin(state(3));
        cosSPsi = cosSLeg*cos(state(3)) + sinSLeg*sin(state(3));
        
        if strcmp(vS.params.vppType, 'Varying')

            tanBetaS = (vS.params.rH + rVPPS)*sinSPsi/...
                (lenStance + (vS.params.rH + rVPPS)*cosSPsi);
            
        else
            
            tanBetaS = (vS.params.rH + vS.params.rVPPS)*sinSPsi/...
                (lenStance + (vS.params.rH + vS.params.rVPPS)*cosSPsi);
            
        end

        tauS = forceSSpring*lenStance*tanBetaS;

        forceX = sinSLeg.*forceSSpring - cosSLeg.*tauS./lenStance;
        forceZ = cosSLeg.*forceSSpring + sinSLeg.*tauS./lenStance;


        % State Dynamics
        dynStruct.stateDotS = tfS.*[state(4);state(5);state(6);...
            (1./vS.params.m).*forceX;...
            (1./vS.params.m).*forceZ - vS.params.g;
            (1./vS.params.J).*(tauS -...
            forceX*vS.params.rH*cos(state(3)) +...
            forceZ*vS.params.rH*sin(state(3)));
            kSLagDot; kSLeadDot];
        
    elseif strcmp(vS.params.model,'BSLIP')
        
        forceX = sinSLeg.*forceSSpring;
        forceZ = cosSLeg.*forceSSpring;
        
        % State Dynamics
        dynStruct.stateDotS = tfS*[state(3); state(4);...
            (1./vS.params.m).*(forceX);...
            (1./vS.params.m).*forceZ-vS.params.g;...
            kSLagDot; kSLeadDot];
        
    end

    % Cost Function
    dynStruct.costS = vS.params.Q(1,1)*(state(1) - humDatHor)^2 +...
        vS.params.Q(2,2)*(state(2) - humDatVert)^2;
    dynStruct.resVertS = abs(vS.params.Q(2,2)*(state(2) - humDatVert));
%     dynStruct.costS = 0;

    % CadADi Function Object
    dynStruct.fS = Function('fS', {state, pS},...
        {dynStruct.stateDotS, dynStruct.costS, dynStruct.resVertS});
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Runga-Kutta 4 %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    
    if strcmp(vS.params.method, 'Multishooting')

        % RK4 Initial Symbolics
        STATE0 = MX.sym('STATE0',length(state));
        PD = MX.sym('PD',length(pD));
        U = MX.sym('U');
        DT = vS.params.dt;

        % Double Support RK4
        STATED = STATE0;
        QD = 0;
        RD = 0;

        for i=1:vS.params.M

            % RK4 Step Calculations from CasADi Function Object

            [k1D, k1qD, k1rD] = dynStruct.fD(STATED,PD,U);
            [k2D, k2qD, k2rD] = dynStruct.fD(STATED + (DT*k1D)/2,PD,U);
            [k3D, k3qD, k3rD] = dynStruct.fD(STATED + (DT*k2D)/2,PD,U);
            [k4D, k4qD, k4rD] = dynStruct.fD(STATED + (DT*k3D),PD,U);

            STATED = STATED + (DT/6)*(k1D + 2*k2D + 2*k3D + k4D);
            QD = QD + (DT/6)*(k1qD + 2*k2qD + 2*k3qD + k4qD);
            RD = RD + (DT/6)*(k1rD + 2*k2rD + 2*k3rD + k4rD);

        end

        % RK4 Full RK4 Integration Function Object
        dynStruct.FD = Function('FD', {STATE0, PD, U}, {STATED, QD},...
            {'state0', 'p', 'u'}, {'stated', 'qd', 'rd'});

        % Single Stance RK4
        PS = MX.sym('PS',length(pS));
        STATES = STATE0;
        QS = 0;
        RS = 0;

        for i=1:vS.params.M

            % RK4 Step Calculations from CasADi Function Object

            [k1S, k1qS, k1rS] = dynStruct.fS(STATES,PS);
            [k2S, k2qS, k2rS] = dynStruct.fS(STATES + (DT*k1S)/2,PS);
            [k3S, k3qS, k3rS] = dynStruct.fS(STATES + (DT*k2S)/2,PS);
            [k4S, k4qS, k4rS] = dynStruct.fS(STATES + (DT*k3S),PS);

            STATES = STATES + (DT/6)*(k1S + 2*k2S + 2*k3S + k4S);
            QS = QS + (DT/6)*(k1qS + 2*k2qS + 2*k3qS + k4qS);
            RS = RS + (DT/6)*(k1rS + 2*k2rS + 2*k3rS + k4rS);

        end

        % RK4 Full RK4 Integration Function Object
        dynStruct.FS = Function('FS', {STATE0, PS}, {STATES, QS},...
            {'state0', 'p'}, {'states', 'qs', 'rs'});
        
    end

end