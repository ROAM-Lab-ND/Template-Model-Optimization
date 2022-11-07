%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcGRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 10 December 2021
% Last Updated: 10 December 2021

% This function is used to back calculate the GRF determined by the
% optimized BSLIP/VPP model

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used subject data storage

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used subject data storage

%% %%%%%%%%%%%%%%%%%%%%%%% Numerical Integration %%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, datS] = calcGRF(varS, datS)

    % Check if subject GRF data is available
    if isfield(datS.data, 'grf')

        % Grab optimized normalized time vector for easier coding
        timeNorm = varS.optims.timeNorm;

        % Interpolate GRF data at optimized time instances
        datS.interp.grfVert = interp1(datS.data.timeGRFNorm,...
            datS.fit.grfVert, timeNorm, 'spline');
        datS.interp.grfHor = interp1(datS.data.timeGRFNorm,...
            datS.fit.grfHor, timeNorm, 'spline');
        
        % Calculate resulting GRF magnitude and angle
        datS.interp.grfNorm = sqrt(datS.interp.grfHor.^2 +...
            datS.interp.grfVert.^2);
        datS.interp.grfAng = atan2(datS.interp.grfHor,...
            datS.interp.grfVert);
        
        % Fix any angle calculations that go beyond 180 degrees
        indsFlip = find(datS.interp.grfAng >= pi/2);
        datS.interp.grfAng(indsFlip) =...
            datS.interp.grfAng(indsFlip) - pi;

        % Grab optimized normalized full time vector for easier coding
        timeNorm = varS.optims.timeNormFull;

        % Interpolate full GRF data at optimized time instances
        datS.interp.grfVertFull = interp1(datS.data.timeGRFNorm,...
            datS.fit.grfVert, timeNorm, 'spline');
        datS.interp.grfHorFull = interp1(datS.data.timeGRFNorm,...
            datS.fit.grfHor, timeNorm, 'spline');
        
        % Calculate resulting full GRF magnitude and angle
        datS.interp.grfNormFull = sqrt(datS.interp.grfHorFull.^2 +...
            datS.interp.grfVertFull.^2);
        datS.interp.grfAngFull = atan2(datS.interp.grfHorFull,...
            datS.interp.grfVertFull);
        
        % Fix any angle calculations that go beyond 180 degrees
        indsFlip = find(datS.interp.grfAngFull >= pi/2);
        datS.interp.grfAngFull(indsFlip) =...
            datS.interp.grfAngFull(indsFlip) - pi;

    end

    % Check if completing calculations for all time instances
    if ~strcmp(varS.params.size, 'Full')

        % Grab length of state vector
        lenVec = length(varS.optims.state(:,1));
        
        % Grab normalized time vector
        timeNorm = varS.optims.timeNorm;

        % Grab state and input vectors
        stateOpt = varS.optims.state;
        optForce = varS.optims.uOpt;

        % Grab number of shooting elements
        phaseDur = varS.params.N;

    else

        % Grab length of full state vector
        lenVec = length(varS.optims.stateFull(:,1));
        
        % Grab normalized full time vector
        timeNorm = varS.optims.timeNormFull;

        % Grab full state and input vectors
        stateOpt = varS.optims.stateFull;
        optForce = repelem(varS.optims.uOpt, varS.params.M);

        % Grab number of shooting and finite elements
        phaseDur = varS.params.N*varS.params.M;

    end
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')
           
        % Grab VP values for easier coding
        rVPP = repelem(varS.optims.rVPP(:), phaseDur);
        rVPP = [rVPP; varS.optims.rVPP(1)];
        
        % Initialize hip torque values
        tauLead = zeros(1, lenVec);
        tauLag = zeros(1, lenVec);
        
        % Initialize hip-leg to VP orientation terms
        tanBetaLead = zeros(1, lenVec);
        tanBetaLag = zeros(1, lenVec);

    end
    
    % Initialize GRF terms
    grf1x = zeros(1, lenVec);
    grf1y = zeros(1, lenVec);
    grf2x = zeros(1, lenVec);
    grf2y = zeros(1, lenVec);
    
    % Grab nominal leg length
    legRight = varS.optims.len0Opt*ones(1, lenVec);
    legLeft = varS.optims.len0Opt*ones(1, lenVec);
    
    % Initialize cycle iterator
    cycle = 1;

    % For each time instance
    for i=1:lenVec

        % Check if using VPP or BSLIP template model
        if strcmp(varS.params.model, 'VPP')

            % Calculate hip position
            xHip = stateOpt(i,1) - varS.params.rH*sin(stateOpt(i,3));
            zHip = stateOpt(i,2) - varS.params.rH*cos(stateOpt(i,3));

        else

            % Calculate hip position
            xHip = stateOpt(i,1);
            zHip = stateOpt(i,2);

        end

        % Check if currently in double support 1
        if ((i-1) < ((4*(cycle - 1) + 1)*phaseDur))

            % Grab current foot positions
            fPos1 = varS.optims.fOptOrig(4*(cycle-1) + 1);
            fPos2 = varS.optims.fOptOrig(4*(cycle-1) + 2);

            % Calculate current leg length
            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);
            len2 = sqrt((fPos2 - xHip)^2 + zHip^2);

            % Store current leg lengths for easier coding
            legLeft(i) = len1;
            legRight(i) = len2;

            % Calculate input and spring forces along legs
            forceOpt = 0;%1000*optForce(i);
            forceSpringLag = 1000*stateOpt(i,end-1)*...
                (varS.optims.len0Opt - len1);
            forceSpringLead = 1000*stateOpt(i,end)*...
                (varS.optims.len0Opt - len2);

            % Calculate leg orientation terms
            sinLegLag = (xHip - fPos1)/len1;
            cosLegLag = zHip/len1;

            sinLegLead = (xHip - fPos2)/len2;
            cosLegLead = zHip/len2;

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate leg-trunk orientation terms
                sinPsiLag = sinLegLag*cos(stateOpt(i,3)) -...
                    cosLegLag*sin(stateOpt(i,3));
                cosPsiLag = cosLegLag*cos(stateOpt(i,3)) +...
                    sinLegLag*sin(stateOpt(i,3));

                sinPsiLead = sinLegLead*cos(stateOpt(i,3)) -...
                    cosLegLead*sin(stateOpt(i,3));
                cosPsiLead = cosLegLead*cos(stateOpt(i,3)) +...
                    sinLegLead*sin(stateOpt(i,3));

                % Calculate hip-leg to VP orientation terms
                tanBetaLag(i) = (varS.params.rH + rVPP(i))*...
                    sinPsiLag/(len1 + (varS.params.rH +...
                    rVPP(i))*cosPsiLag);
                tanBetaLead(i) = (varS.params.rH + rVPP(i))*...
                    sinPsiLead/(len2 + (varS.params.rH +...
                    rVPP(i))*cosPsiLead);

                % Calculate hip torque terms
                tauLag(i) = forceSpringLag*len1*tanBetaLag(i);
                tauLead(i) = forceSpringLead*len2*tanBetaLead(i);

                % Calculate overall resulting force terms
                forceXLag = sinLegLag*forceSpringLag -...
                    cosLegLag*tauLag(i)/len1;
                forceXLead = sinLegLead*forceSpringLead -...
                    cosLegLead*tauLead(i)/len2;

                forceYLag = cosLegLag*forceSpringLag +...
                    sinLegLag*tauLag(i)/len1;
                forceYLead = cosLegLead*forceSpringLead +...
                    sinLegLead*tauLead(i)/len2;

            else

                % Calculate overall resulting force terms
                forceXLag = sinLegLag*(forceSpringLag + forceOpt);
                forceXLead = sinLegLead*forceSpringLead;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt);
                forceYLead = cosLegLead*forceSpringLead;

            end

            % Store GRF terms for future use
            grf1y(i) = forceYLag;
            grf2y(i) = forceYLead;

            grf1x(i) = forceXLag;
            grf2x(i) = forceXLead;

        % Check if currently in single support 1
        elseif ((i-1) < (4*(cycle - 1) + 2)*phaseDur)

            % Grab current foot position
            fPos1 = varS.optims.fOptOrig(4*(cycle-1) + 2);

            % Calculate current leg length
            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            % Store current leg length
            legRight(i) = len1;

            % Calculate spring force along leg
            forceSpring = 1000*stateOpt(i,end)*...
                (varS.optims.len0Opt - len1);

            % Calculate leg orientation terms
            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate leg-trunk orientation terms
                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                % Calculate hip-leg to VP orientation terms
                tanBetaLead(i) = (varS.params.rH + rVPP(i))*...
                    sinPsi/(len1 + (varS.params.rH +...
                    rVPP(i))*cosPsi);

                % Calculate hip torque terms
                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end
            
            % Store GRF terms
            grf1y(i) = 0;
            grf2y(i) = forceY;

            grf1x(i) = 0;
            grf2x(i) = forceX;

        % Check if currently in double support 2
        elseif ((i-1) < (4*(cycle - 1) + 3)*phaseDur)

            % Grab current foot positions
            fPos1 = varS.optims.fOptOrig(4*(cycle-1) + 4);
            fPos2 = varS.optims.fOptOrig(4*(cycle-1) + 3);

            % Calculate current leg lengths
            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);
            len2 = sqrt((fPos2 - xHip)^2 + zHip^2);

            % Store current leg lengths
            legLeft(i) = len1;
            legRight(i) = len2;

            % Calculate input and spring forces along legs
            forceOpt = 0;%1000*optForce(i-phaseDur);
            forceSpringLag = 1000*stateOpt(i,end-1)*...
                (varS.optims.len0Opt - len2);
            forceSpringLead = 1000*stateOpt(i,end)*...
                (varS.optims.len0Opt - len1);

            % Calculate leg orientation terms
            sinLegLag = (xHip - fPos2)/len2;
            cosLegLag = zHip/len2;

            sinLegLead = (xHip - fPos1)/len1;
            cosLegLead = zHip/len1;

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate leg-trunk orientation terms
                sinPsiLag = sinLegLag*cos(stateOpt(i,3)) -...
                    cosLegLag*sin(stateOpt(i,3));
                cosPsiLag = cosLegLag*cos(stateOpt(i,3)) +...
                    sinLegLag*sin(stateOpt(i,3));

                sinPsiLead = sinLegLead*cos(stateOpt(i,3)) -...
                    cosLegLead*sin(stateOpt(i,3));
                cosPsiLead = cosLegLead*cos(stateOpt(i,3)) +...
                    sinLegLead*sin(stateOpt(i,3));

                % Calculate hip-leg to VP orientation terms
                tanBetaLag(i) = (varS.params.rH + rVPP(i))*...
                    sinPsiLag/(len2 + (varS.params.rH +...
                    rVPP(i))*cosPsiLag);
                tanBetaLead(i) = (varS.params.rH + rVPP(i))*...
                    sinPsiLead/(len1 + (varS.params.rH +...
                    rVPP(i))*cosPsiLead);

                % Calculate hip torque terms
                tauLag(i) = forceSpringLag*len2*tanBetaLag(i);
                tauLead(i) = forceSpringLead*len1*tanBetaLead(i);

                % Calculate overall resulting force terms
                forceXLag = sinLegLag*(forceSpringLag + forceOpt) -...
                    cosLegLag*tauLag(i)/len2;
                forceXLead = sinLegLead*forceSpringLead -...
                    cosLegLead*tauLead(i)/len1;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt) +...
                    sinLegLag*tauLag(i)/len2;
                forceYLead = cosLegLead*forceSpringLead +...
                    sinLegLead*tauLead(i)/len1;

            else

                % Calculate overall resulting force terms
                forceXLag = sinLegLag*(forceSpringLag + forceOpt);
                forceXLead = sinLegLead*forceSpringLead;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt);
                forceYLead = cosLegLead*forceSpringLead;

            end

            % Store GRF terms
            grf1y(i) = forceYLead;
            grf2y(i) = forceYLag;

            grf1x(i) = forceXLead;
            grf2x(i) = forceXLag;

        % Check if currently in single support 2
        elseif ((i-1) < (4*(cycle - 1) + 4)*phaseDur)

            % Grab current foot position
            fPos1 = varS.optims.fOptOrig(4*(cycle-1) + 5);

            % Calculate current leg length
            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            % Store current leg length
            legLeft(i) = len1;

            % Calculate spring force along leg
            forceSpring = 1000*stateOpt(i,end)*...
                (varS.optims.len0Opt - len1);

            % Calculate leg orientation terms
            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate leg-trunk orientation terms
                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                % Calculate hip-leg to VP orientation term
                tanBetaLead(i) = (varS.params.rH + rVPP(i))*...
                    sinPsi/(len1 + (varS.params.rH +...
                    rVPP(i))*cosPsi);

                % Calculate hip torque term
                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end

            % Store GRF terms
            grf1y(i) = forceY;
            grf2y(i) = 0;

            grf1x(i) = forceX;
            grf2x(i) = 0;

        else

            % Grab current foot position
            fPos1 = varS.optims.fOptOrig(4*(cycle-1) + 5);

            % Calculate current leg length
            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            % Store current leg length
            legLeft(i) = len1;

            % Calculate spring force along leg
            forceSpring = 1000*stateOpt(i,end)*...
                (varS.optims.len0Opt - len1);

            % Calculate leg orientation terms
            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            % Check if using VPP or BSLIP template model
            if strcmp(varS.params.model, 'VPP')

                % Calculate leg-trunk orientation terms
                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                % Calculate hip-leg to VP orientation term
                tanBetaLead(i) = (varS.params.rH + rVPP(end))*...
                    sinPsi/(len1 + (varS.params.rH +...
                    rVPP(end))*cosPsi);

                % Calculate hip torque term
                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                % Calculate overall resulting force terms
                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end

            % Store GRF terms
            grf1y(i) = forceY;
            grf2y(i) = 0;

            grf1x(i) = forceX;
            grf2x(i) = 0;

        end
        
        % Check if current gait cycle has ended
        if (i/cycle) > 4*phaseDur
            
            % Update gait cycle iterator
            cycle = cycle + 1;
            
        end

    end
    
    % Check if working with base or full time instance vector
    if ~strcmp(varS.params.size, 'Full')

        % Store GRF values
        varS.optims.grfVert = [grf1y', grf2y'];
        varS.optims.grfHor = [grf1x', grf2x'];
        varS.optims.grfNorm = sqrt(varS.optims.grfVert.^2 +...
            varS.optims.grfHor.^2);

        % Store GRF angles
        varS.optims.grfAng = atan2(varS.optims.grfHor, varS.optims.grfVert);
        indsFlip = find(varS.optims.grfAng >= pi/2);
        varS.optims.grfAng(indsFlip) = varS.optims.grfAng(indsFlip) - pi;
        
        % Store optimized leg lengths
        varS.optims.legL = legLeft;
        varS.optims.legR = legRight;
        
        % Check if using VPP or BSLIP template model
        if strcmp(varS.params.model, 'VPP')
            
            % Store hip torque values
            varS.optims.tauLag = tauLag;
            varS.optims.tauLead = tauLead;
            
            % Store hip-leg to VP orientation values
            varS.optims.tanBetaLag = tanBetaLag;
            varS.optims.tanBetaLead = tanBetaLead;
            
        end

    else

        % Store full GRF values
        varS.optims.grfVertFull = [grf1y', grf2y'];
        varS.optims.grfHorFull = [grf1x', grf2x'];
        varS.optims.grfNormFull = sqrt(varS.optims.grfVertFull.^2 +...
            varS.optims.grfHorFull.^2);
        varS.optims.timeNormFull = timeNorm;

        % Store full GRF angles
        varS.optims.grfAngFull = atan2(varS.optims.grfHorFull,...
            varS.optims.grfVertFull);
        indsFlip = find(varS.optims.grfAngFull >= pi/2);
        varS.optims.grfAngFull(indsFlip) =...
            varS.optims.grfAngFull(indsFlip) - pi;
        
        % Store full optimized leg lengths
        varS.optims.legLFull = legLeft;
        varS.optims.legRFull = legRight;
        
        % Check if using VPP or BSLIP template models
        if strcmp(varS.params.model, 'VPP')
            
            % Store hip torque values
            varS.optims.tauLagFull = tauLag;
            varS.optims.tauLeadFull = tauLead;
            
            % Store hip-leg to VP orientation values
            varS.optims.tanBetaLagFull = tanBetaLag;
            varS.optims.tanBetaLeadFull = tanBetaLead;
            
        end

    end

end
