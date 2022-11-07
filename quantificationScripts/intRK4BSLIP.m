%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotBSLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
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

function stateDot = intRK4BSLIP(varS, state, params, u, phase)

    % Check if in double support or single support
    if phase == 1
        
        % Parse out pertinate non-current-state parameters
        stateD = params(1:6);
        thetaD = params(7);
        kDLagDot = params(8);
        kDLeadDot = params(9);
        footDLag = params(10);
        tfD = params(11);
        
        %%%%%%%%%%%%%%%%%%%%% Double Support Dynamics %%%%%%%%%%%%%%%%%%%%%
        
        % Calculate initial hip position
        xHip0D = stateD(1);
        zHip0D = stateD(2);
        hip0D = [xHip0D; zHip0D];

        % Calculate current hip position
        xHipD = state(1);
        zHipD = state(2);
        hipD = [xHipD; zHipD];

        % Calculate stride length
        stride = hip0D(1) + hip0D(2).*tan(thetaD) - footDLag;

        % Calculate current leg lengths
        lenDLag = sqrt((hipD(1) - footDLag).^2 + hipD(2).^2);
        lenDLead = sqrt((stride - (hipD(1) - footDLag)).^2 + hipD(2).^2);

        % Calculate spring force along legs
        forceDSpringLag = 1000*state(end-1).*...
            (varS.optims.len0Opt - lenDLag);
        forceDSpringLead = 1000*state(end).*...
            (varS.optims.len0Opt - lenDLead);

        % Calculate leg orientation terms
        sinDLegLag = (hipD(1) - footDLag)./lenDLag;
        cosDLegLag = hipD(2)./lenDLag;
        sinDLegLead = (stride - (hipD(1) - footDLag))./lenDLead;
        cosDLegLead = hipD(2)./lenDLead;

        % Calculate overall resulting force terms
        forceXLag = sinDLegLag.*(forceDSpringLag + 1000*u);
        forceXLead = -sinDLegLead.*forceDSpringLead;

        forceZLag = cosDLegLag.*(forceDSpringLag + 1000*u);    
        forceZLead = cosDLegLead.*forceDSpringLead;

        % Calculate state dynamics
        stateDot = tfD.*[state(3); state(4);...
            (1./varS.params.m).*(forceXLag + forceXLead);...
            (1./varS.params.m).*(forceZLag + forceZLead)-varS.params.g;...
            kDLagDot; kDLeadDot];
        
    else
        
        % Parse out pertinate non-current-state parameters
        kSLagDot = params(1);
        kSLeadDot = params(2);
        footSStance = params(3);
        tfS = params(4);
        
        %%%%%%%%%%%%%%%%%%%%% Single Support Dynamics %%%%%%%%%%%%%%%%%%%%%

        % Calculate current hip position
        xHipS = state(1);
        zHipS = state(2);
        hipS = [xHipS; zHipS];

        % Calculate leg length
        lenStance = sqrt((hipS(1) - footSStance).^2 + hipS(2).^2);

        % Calculate spring force along leg
        forceSSpring = 1000*state(end).*(varS.optims.len0Opt - lenStance);

        % Calculate leg orientation terms
        sinSLeg = (hipS(1) - footSStance)./lenStance;
        cosSLeg = hipS(2)./lenStance;
        
        % Calculate overall resulting force terms
        forceX = sinSLeg.*forceSSpring;
        forceZ = cosSLeg.*forceSSpring;

        % State Dynamics
        stateDot = tfS*[state(3); state(4);...
            (1./varS.params.m).*forceX;...
            (1./varS.params.m).*forceZ - varS.params.g;...
            kSLagDot; kSLeadDot];
    
    end    

end