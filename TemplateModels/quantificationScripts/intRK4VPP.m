%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotVPP %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 13 July 2021
% Last Updated: 13 July 2021

% This function is used to plot the dynamics of an optimized BSLIP model
% using the dynamics

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function stateDot = intRK4VPP(vS, state, params, u, phase)

    if phase == 1
        
        % Double Stance Dynamics
        stateD = params(1:8);
        thetaD = params(9);
        kDLagDot = params(10);
        kDLeadDot = params(11);
        footDLag = params(12);
        tfD = params(13);
        rVPP = params(end);

        % ##################### Double Support Dynamics ##################### %
    
        xHip0D = stateD(1) - vS.params.rH*sin(stateD(3));
        yHip0D = stateD(2) - vS.params.rH*cos(stateD(3));
        hip0D = [xHip0D;yHip0D];

        xHipD = state(1) - vS.params.rH*sin(state(3));
        yHipD = state(2) - vS.params.rH*cos(state(3));
        hipD = [xHipD;yHipD];

        stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag; % Stride Length

        lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2); % Lag Leg Length
        lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2); % Lead Leg Length

        forceDSpringLag = 1000*state(7).*(vS.optims.len0Opt - lenDLag); % Lag Leg Spring Force
        forceDSpringLead = 1000*state(8).*(vS.optims.len0Opt - lenDLead); % Lead Leg Spring Force

        sinDLegLag = (hipD(1) - footDLag)./lenDLag; % X-Axis Weight Lag Leg
        cosDLegLag = hipD(2)./lenDLag; % Y-Axis Weight Lag Leg
        sinDLegLead = -(stride - (hipD(1) - footDLag))./lenDLead; % X-Axis Weight Lead Leg
        cosDLegLead = hipD(2)./lenDLead; % Y-Axis Weight Lead Leg

        sinDPsiLag = sinDLegLag*cos(state(3)) - cosDLegLag*sin(state(3));
        cosDPsiLag = cosDLegLag*cos(state(3)) + sinDLegLag*sin(state(3));
        sinDPsiLead = sinDLegLead*cos(state(3)) - cosDLegLead*sin(state(3));
        cosDPsiLead = cosDLegLead*cos(state(3)) + sinDLegLead*sin(state(3));

        tanBetaLag = (vS.params.rH + rVPP)*sinDPsiLag/...
            (lenDLag + (vS.params.rH + rVPP)*cosDPsiLag);    
        tanBetaLead = (vS.params.rH + rVPP)*sinDPsiLead/...
            (lenDLead + (vS.params.rH + rVPP)*cosDPsiLead);

        tauLead = forceDSpringLead*lenDLead*tanBetaLead;
        tauLag = forceDSpringLag*lenDLag*tanBetaLag;

        forceXLag = sinDLegLag.*forceDSpringLag - cosDLegLag.*tauLag./lenDLag;
        forceXLead = sinDLegLead.*forceDSpringLead - cosDLegLead.*tauLead./lenDLead;

        forceYLag = cosDLegLag.*forceDSpringLag + sinDLegLag.*tauLag./lenDLag;    
        forceYLead = cosDLegLead.*forceDSpringLead + sinDLegLead.*tauLead./lenDLead;

        % State Dynamics
        stateDot = tfD.*[state(4); state(5); state(6);...
            (1./vS.params.m).*(forceXLag + forceXLead + 1000*u*sinDLegLag);...
            (1./vS.params.m).*(forceYLead + forceYLag + 1000*u*cosDLegLag)-...
            vS.params.g;...
            (1./vS.params.J).*(tauLag + tauLead -...
            (forceXLag + forceXLead)*vS.params.rH*cos(state(3)) +...
            (forceYLag + forceYLead)*vS.params.rH*sin(state(3)));...
            kDLagDot; kDLeadDot];
        
    else
        
        % Single Support Dynamics
        kSLagDot = params(1);
        kSLeadDot = params(2);
        footSStance = params(3);
        tfS = params(4);
        rVPP = params(end);
        
        % ##################### Single Support Dynamics ##################### %

        xHipS = state(1) - vS.params.rH*sin(state(3));
        yHipS = state(2) - vS.params.rH*cos(state(3));
        hipS = [xHipS;yHipS];

        lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2); % Stance Leg Length

        forceSSpring = 1000*state(8).*(vS.optims.len0Opt - lenStance); % Spring Force

        sinSLeg = (hipS(1) - footSStance)./lenStance; % X-Axis Weight
        cosSLeg = hipS(2)./lenStance; % Y-Axis Weight

        sinSPsi = sinSLeg*cos(state(3)) - cosSLeg*sin(state(3));
        cosSPsi = cosSLeg*cos(state(3)) + sinSLeg*sin(state(3));

        tanBetaS = (vS.params.rH + rVPP)*sinSPsi/...
            (lenStance + (vS.params.rH + rVPP)*cosSPsi);

        tauS = forceSSpring*lenStance*tanBetaS;

        forceX = sinSLeg.*forceSSpring - cosSLeg.*tauS./lenStance;
        forceY = cosSLeg.*forceSSpring + sinSLeg.*tauS./lenStance;

        % State Dynamics
        stateDot = tfS*[state(4);state(5);state(6);...
            (1./vS.params.m).*forceX;...
            (1./vS.params.m).*forceY - vS.params.g;
            (1./vS.params.J).*(tauS -...
            forceX*vS.params.rH*cos(state(3)) +...
            forceY*vS.params.rH*sin(state(3)));...
            kSLagDot;kSLeadDot];
    
    end    

end