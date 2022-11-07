%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotVPP %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 13 July 2021
% Last Updated: 13 July 2021

% This function is used to plot the dynamics of an optimized BSLIP model
% using the dynamics

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   state - State values at current time step
%   params - Non-state values needed for propogation at current time step
%   u - Input control value at current time step
%   phase - Phase indicator

% OUTPUTS:
%   stateDot - State ROC values based on model dynamics

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function stateDot = intRK4VPP(varS, state, params, u, phase)

    % Check if in double support or single support
    if phase == 1
        
        % Parse out pertinate non-current-state parameters 
        stateD = params(1:8);
        thetaD = params(9);
        kDLagDot = params(10);
        kDLeadDot = params(11);
        footDLag = params(12);
        tfD = params(13);
        rVPP = params(end);

        %%%%%%%%%%%%%%%%%%%%% Double Support Dynamics %%%%%%%%%%%%%%%%%%%%%
    
        % Calculate initial hip position
        xHip0D = stateD(1) - varS.params.rH*sin(stateD(3));
        zHip0D = stateD(2) - varS.params.rH*cos(stateD(3));
        hip0D = [xHip0D; zHip0D];

        % Calculate current hip position
        xHipD = state(1) - varS.params.rH*sin(state(3));
        zHipD = state(2) - varS.params.rH*cos(state(3));
        hipD = [xHipD; zHipD];

        % Calculate stride length
        stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag;

        % Calculate current leg lengths
        lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2);
        lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2);

        % Calculate spring force along legs
        forceDSpringLag = 1000*state(7).*(varS.optims.len0Opt - lenDLag);
        forceDSpringLead = 1000*state(8).*(varS.optims.len0Opt - lenDLead);

        % Calculate leg orientation terms
        sinDLegLag = (hipD(1) - footDLag)./lenDLag;
        cosDLegLag = hipD(2)./lenDLag;
        sinDLegLead = -(stride - (hipD(1) - footDLag))./lenDLead;
        cosDLegLead = hipD(2)./lenDLead;

        % Calculate leg-trunk orientation terms
        sinDPsiLag = sinDLegLag*cos(state(3)) - cosDLegLag*sin(state(3));
        cosDPsiLag = cosDLegLag*cos(state(3)) + sinDLegLag*sin(state(3));
        sinDPsiLead = sinDLegLead*cos(state(3)) -...
            cosDLegLead*sin(state(3));
        cosDPsiLead = cosDLegLead*cos(state(3)) +...
            sinDLegLead*sin(state(3));

        % Calculate hip-leg to VP orientation terms
        tanBetaLag = (varS.params.rH + rVPP)*sinDPsiLag/...
            (lenDLag + (varS.params.rH + rVPP)*cosDPsiLag);    
        tanBetaLead = (varS.params.rH + rVPP)*sinDPsiLead/...
            (lenDLead + (varS.params.rH + rVPP)*cosDPsiLead);

        % Calculate hip torque terms
        tauLead = forceDSpringLead*lenDLead*tanBetaLead;
        tauLag = forceDSpringLag*lenDLag*tanBetaLag;

        % Calculate overall resulting force terms
        forceXLag = sinDLegLag.*forceDSpringLag -...
            cosDLegLag.*tauLag./lenDLag;
        forceXLead = sinDLegLead.*forceDSpringLead -...
            cosDLegLead.*tauLead./lenDLead;

        forceZLag = cosDLegLag.*forceDSpringLag +...
            sinDLegLag.*tauLag./lenDLag;    
        forceZLead = cosDLegLead.*forceDSpringLead +...
            sinDLegLead.*tauLead./lenDLead;

        % Calculate state dynamics
        stateDot = tfD.*[state(4); state(5); state(6);...
            (1./varS.params.m).*...
            (forceXLag + forceXLead + 1000*u*sinDLegLag);...
            (1./varS.params.m).*...
            (forceZLead + forceZLag + 1000*u*cosDLegLag)-...
            varS.params.g;...
            (1./varS.params.J).*(tauLag + tauLead -...
            (forceXLag + forceXLead)*varS.params.rH*cos(state(3)) +...
            (forceZLag + forceZLead)*varS.params.rH*sin(state(3)));...
            kDLagDot; kDLeadDot];
        
    else
        
        % Parse out pertinate non-current-state parameters
        kSLagDot = params(1);
        kSLeadDot = params(2);
        footSStance = params(3);
        tfS = params(4);
        rVPP = params(end);
        
        %%%%%%%%%%%%%%%%%%%%% Single Support Dynamics %%%%%%%%%%%%%%%%%%%%%

        % Calculate current hip position
        xHipS = state(1) - varS.params.rH*sin(state(3));
        zHipS = state(2) - varS.params.rH*cos(state(3));
        hipS = [xHipS; zHipS];

        % Calculate current leg length
        lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2);

        % Calculate spring force along leg
        forceSSpring = 1000*state(8).*(varS.optims.len0Opt - lenStance);

        % Calculate leg orientation terms
        sinSLeg = (hipS(1) - footSStance)./lenStance;
        cosSLeg = hipS(2)./lenStance;

        % Calculate leg-trunk orientation terms
        sinSPsi = sinSLeg*cos(state(3)) - cosSLeg*sin(state(3));
        cosSPsi = cosSLeg*cos(state(3)) + sinSLeg*sin(state(3));

        % Calculate hip-leg to VP orientation terms
        tanBetaS = (varS.params.rH + rVPP)*sinSPsi/...
            (lenStance + (varS.params.rH + rVPP)*cosSPsi);

        % Calculate hip torque term
        tauS = forceSSpring*lenStance*tanBetaS;

        % Calculate overall resulting force terms
        forceX = sinSLeg.*forceSSpring - cosSLeg.*tauS./lenStance;
        forceZ = cosSLeg.*forceSSpring + sinSLeg.*tauS./lenStance;

        % Calculate state dynamics
        stateDot = tfS*[state(4);state(5);state(6);...
            (1./varS.params.m).*forceX;...
            (1./varS.params.m).*forceZ - varS.params.g;
            (1./varS.params.J).*(tauS -...
            forceX*varS.params.rH*cos(state(3)) +...
            forceZ*varS.params.rH*sin(state(3)));...
            kSLagDot;kSLeadDot];
    
    end    

end