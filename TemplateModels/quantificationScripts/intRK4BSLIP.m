%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotBSLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 13 July 2021
% Last Updated: 13 July 2021

% This function is used to plot the dynamics of an optimized BSLIP model
% using the dynamics

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function stateDot = intRK4BSLIP(vS, state, params, u, phase)

    if phase == 1
        
        % Double Stance Dynamics
        stateD = params(1:6);
        thetaD = params(7);
        kDLagDot = params(8);
        kDLeadDot = params(9);
        footDLag = params(10);
        tfD = params(11);
        
        % ##################### Double Support Dynamics ##################### %
        
        xHip0D = stateD(1);
        yHip0D = stateD(2);
        hip0D = [xHip0D; yHip0D];

        xHipD = state(1);
        yHipD = state(2);
        hipD = [xHipD; yHipD];

        stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag; % Stride Length

        lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2); % Lag Leg Length
        lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2); % Lead Leg Length

        forceDSpringLag = 1000*state(end-1).*(vS.optims.len0Opt - lenDLag); % Lag Leg Spring Force
        forceDSpringLead = 1000*state(end).*(vS.optims.len0Opt - lenDLead); % Lead Leg Spring Force

        sinDLegLag = (hipD(1) - footDLag)./lenDLag; % X-Axis Weight Lag Leg
        cosDLegLag = hipD(2)./lenDLag; % Y-Axis Weight Lag Leg
        sinDLegLead = (stride - (hipD(1) - footDLag))./lenDLead; % X-Axis Weight Lead Leg
        cosDLegLead = hipD(2)./lenDLead; % Y-Axis Weight Lead Leg

        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u);
        forceXLead = -sinDLegLead.*forceDSpringLead;

        forceYLag = cosDLegLag.*(forceDSpringLag + 1000*u);    
        forceYLead = cosDLegLead.*forceDSpringLead;

        % State Dynamics
        stateDot = tfD.*[state(3); state(4);...
            (1./vS.params.m).*(forceXLag + forceXLead);...
            (1./vS.params.m).*(forceYLag + forceYLead)-vS.params.g;...
            kDLagDot; kDLeadDot];
        
    else
        
        % Single Support Parameters
        kSLagDot = params(1);
        kSLeadDot = params(2);
        footSStance = params(3);
        tfS = params(4);
        
        % ##################### Single Support Dynamics ##################### %

        xHipS = state(1);
        yHipS = state(2);
        hipS = [xHipS; yHipS];

        lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2); % Stance Leg Length

        forceSSpring = 1000*state(end).*(vS.optims.len0Opt - lenStance); % Spring Force

        sinSLeg = (hipS(1) - footSStance)./lenStance; % X-Axis Weight
        cosSLeg = hipS(2)./lenStance; % Y-Axis Weight
        
        forceX = sinSLeg.*forceSSpring;
        forceY = cosSLeg.*forceSSpring;

        % State Dynamics
        stateDot = tfS*[state(3); state(4);...
            (1./vS.params.m).*forceX;...
            (1./vS.params.m).*forceY - vS.params.g;...
            kSLagDot; kSLeadDot];
    
    end    

end