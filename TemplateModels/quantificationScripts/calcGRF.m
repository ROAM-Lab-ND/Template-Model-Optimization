%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcGRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 10 December 2021
% Last Updated: 10 December 2021

% This function is used to back calculate the GRF determined by the
% optimized BSLIP/VPP model

%% %%%%%%%%%%%%%%%%%%%%%%% Numerical Integration %%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, dS] = calcGRF(vS, dS)

    if isfield(dS.data,'grf')

        timeNorm = vS.optims.timeNorm;

        dS.interp.grfVert = interp1(dS.data.timeGRFNorm,...
            dS.fit.grfVert, timeNorm,'spline');
        dS.interp.grfHor = interp1(dS.data.timeGRFNorm,...
            dS.fit.grfHor, timeNorm,'spline');
        
        dS.interp.grfNorm = sqrt(dS.interp.grfHor.^2 +...
            dS.interp.grfVert.^2);
        dS.interp.grfAng = atan2(dS.interp.grfHor,...
            dS.interp.grfVert);
        
        indsFlip = find(dS.interp.grfAng >= pi/2);
        dS.interp.grfAng(indsFlip) =...
            dS.interp.grfAng(indsFlip) - pi;

        timeNorm = vS.optims.timeNormFull;

        dS.interp.grfVertFull = interp1(dS.data.timeGRFNorm,...
            dS.fit.grfVert, timeNorm,'spline');
        dS.interp.grfHorFull = interp1(dS.data.timeGRFNorm,...
            dS.fit.grfHor, timeNorm,'spline');
        
        dS.interp.grfNormFull = sqrt(dS.interp.grfHorFull.^2 +...
            dS.interp.grfVertFull.^2);
        dS.interp.grfAngFull = atan2(dS.interp.grfHorFull,...
            dS.interp.grfVertFull);
        
        indsFlip = find(dS.interp.grfAngFull >= pi/2);
        dS.interp.grfAngFull(indsFlip) =...
            dS.interp.grfAngFull(indsFlip) - pi;

    end

    if ~strcmp(vS.params.size,'Full')

        lenVec = length(vS.optims.state(:,1));
        timeNorm = vS.optims.timeNorm;

        stateOpt = vS.optims.state;
        optForce = vS.optims.uOpt;

        phaseDur = vS.params.N;

    else

        lenVec = length(vS.optims.stateFull(:,1));
        timeNorm = vS.optims.timeNormFull;

        stateOpt = vS.optims.stateFull;
        optForce = repelem(vS.optims.uOpt, vS.params.M);

        phaseDur = vS.params.N*vS.params.M;

    end
    
    if strcmp(vS.params.model,'VPP')
           
        rVPP = repelem(vS.optims.rVPP(:), phaseDur);
        rVPP = [rVPP; vS.optims.rVPP(1)];
        
        tauLead = zeros(1, lenVec);
        tauLag = zeros(1, lenVec);
        
        tanBetaLead = zeros(1, lenVec);
        tanBetaLag = zeros(1, lenVec);

    end
    
    % GRF
    grf1x = zeros(1, lenVec);
    grf1y = zeros(1, lenVec);
    grf2x = zeros(1, lenVec);
    grf2y = zeros(1, lenVec);
    
    % Legs
    legRight = vS.optims.len0Opt*ones(1, lenVec);
    legLeft = vS.optims.len0Opt*ones(1, lenVec);
    
    cycle = 1;

    for i=1:lenVec

        if strcmp(vS.params.model,'VPP')

            xHip = stateOpt(i,1) - vS.params.rH*sin(stateOpt(i,3));
            zHip = stateOpt(i,2) - vS.params.rH*cos(stateOpt(i,3));

        else

            xHip = stateOpt(i,1);
            zHip = stateOpt(i,2);

        end

        if ((i-1) < ((4*(cycle - 1) + 1)*phaseDur))

            fPos1 = vS.optims.fOptOrig(4*(cycle-1) + 1);
            fPos2 = vS.optims.fOptOrig(4*(cycle-1) + 2);

            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);
            len2 = sqrt((fPos2 - xHip)^2 + zHip^2);

            legLeft(i) = len1;
            legRight(i) = len2;

            forceOpt = 0;%1000*optForce(i);
            forceSpringLag = 1000*stateOpt(i,end-1)*...
                (vS.optims.len0Opt - len1);
            forceSpringLead = 1000*stateOpt(i,end)*...
                (vS.optims.len0Opt - len2);

            sinLegLag = (xHip - fPos1)/len1;
            cosLegLag = zHip/len1;

            sinLegLead = (xHip - fPos2)/len2;
            cosLegLead = zHip/len2;

            if strcmp(vS.params.model,'VPP')

                sinPsiLag = sinLegLag*cos(stateOpt(i,3)) -...
                    cosLegLag*sin(stateOpt(i,3));
                cosPsiLag = cosLegLag*cos(stateOpt(i,3)) +...
                    sinLegLag*sin(stateOpt(i,3));

                sinPsiLead = sinLegLead*cos(stateOpt(i,3)) -...
                    cosLegLead*sin(stateOpt(i,3));
                cosPsiLead = cosLegLead*cos(stateOpt(i,3)) +...
                    sinLegLead*sin(stateOpt(i,3));

                tanBetaLag(i) = (vS.params.rH + rVPP(i))*...
                    sinPsiLag/(len1 + (vS.params.rH +...
                    rVPP(i))*cosPsiLag);
                tanBetaLead(i) = (vS.params.rH + rVPP(i))*...
                    sinPsiLead/(len2 + (vS.params.rH +...
                    rVPP(i))*cosPsiLead);

                tauLag(i) = forceSpringLag*len1*tanBetaLag(i);
                tauLead(i) = forceSpringLead*len2*tanBetaLead(i);

                forceXLag = sinLegLag*forceSpringLag -...
                    cosLegLag*tauLag(i)/len1;
                forceXLead = sinLegLead*forceSpringLead -...
                    cosLegLead*tauLead(i)/len2;

                forceYLag = cosLegLag*forceSpringLag +...
                    sinLegLag*tauLag(i)/len1;
                forceYLead = cosLegLead*forceSpringLead +...
                    sinLegLead*tauLead(i)/len2;

            else

                forceXLag = sinLegLag*(forceSpringLag + forceOpt);
                forceXLead = sinLegLead*forceSpringLead;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt);
                forceYLead = cosLegLead*forceSpringLead;

            end

            grf1y(i) = forceYLag;
            grf2y(i) = forceYLead;

            grf1x(i) = forceXLag;
            grf2x(i) = forceXLead;

        elseif ((i-1) < (4*(cycle - 1) + 2)*phaseDur)

            fPos1 = vS.optims.fOptOrig(4*(cycle-1) + 2);

            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            legRight(i) = len1;

            forceSpring = 1000*stateOpt(i,end)*...
                (vS.optims.len0Opt - len1);

            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            if strcmp(vS.params.model,'VPP')

                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                tanBetaLead(i) = (vS.params.rH + rVPP(i))*...
                    sinPsi/(len1 + (vS.params.rH +...
                    rVPP(i))*cosPsi);

                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end


            grf1y(i) = 0;
            grf2y(i) = forceY;

            grf1x(i) = 0;
            grf2x(i) = forceX;

        elseif ((i-1) < (4*(cycle - 1) + 3)*phaseDur)

            fPos1 = vS.optims.fOptOrig(4*(cycle-1) + 4);
            fPos2 = vS.optims.fOptOrig(4*(cycle-1) + 3);

            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);
            len2 = sqrt((fPos2 - xHip)^2 + zHip^2);

            legLeft(i) = len1;
            legRight(i) = len2;

            forceOpt = 0;%1000*optForce(i-phaseDur);
            forceSpringLag = 1000*stateOpt(i,end-1)*...
                (vS.optims.len0Opt - len2);
            forceSpringLead = 1000*stateOpt(i,end)*...
                (vS.optims.len0Opt - len1);

            sinLegLag = (xHip - fPos2)/len2;
            cosLegLag = zHip/len2;

            sinLegLead = (xHip - fPos1)/len1;
            cosLegLead = zHip/len1;

            if strcmp(vS.params.model,'VPP')

                sinPsiLag = sinLegLag*cos(stateOpt(i,3)) -...
                    cosLegLag*sin(stateOpt(i,3));
                cosPsiLag = cosLegLag*cos(stateOpt(i,3)) +...
                    sinLegLag*sin(stateOpt(i,3));

                sinPsiLead = sinLegLead*cos(stateOpt(i,3)) -...
                    cosLegLead*sin(stateOpt(i,3));
                cosPsiLead = cosLegLead*cos(stateOpt(i,3)) +...
                    sinLegLead*sin(stateOpt(i,3));

                tanBetaLag(i) = (vS.params.rH + rVPP(i))*...
                    sinPsiLag/(len2 + (vS.params.rH +...
                    rVPP(i))*cosPsiLag);
                tanBetaLead(i) = (vS.params.rH + rVPP(i))*...
                    sinPsiLead/(len1 + (vS.params.rH +...
                    rVPP(i))*cosPsiLead);

                tauLag(i) = forceSpringLag*len2*tanBetaLag(i);
                tauLead(i) = forceSpringLead*len1*tanBetaLead(i);

                forceXLag = sinLegLag*(forceSpringLag + forceOpt) -...
                    cosLegLag*tauLag(i)/len2;
                forceXLead = sinLegLead*forceSpringLead -...
                    cosLegLead*tauLead(i)/len1;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt) +...
                    sinLegLag*tauLag(i)/len2;
                forceYLead = cosLegLead*forceSpringLead +...
                    sinLegLead*tauLead(i)/len1;

            else

                forceXLag = sinLegLag*(forceSpringLag + forceOpt);
                forceXLead = sinLegLead*forceSpringLead;

                forceYLag = cosLegLag*(forceSpringLag + forceOpt);
                forceYLead = cosLegLead*forceSpringLead;

            end

            grf1y(i) = forceYLead;
            grf2y(i) = forceYLag;

            grf1x(i) = forceXLead;
            grf2x(i) = forceXLag;

        elseif ((i-1) < (4*(cycle - 1) + 4)*phaseDur)

            fPos1 = vS.optims.fOptOrig(4*(cycle-1) + 5);

            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            legLeft(i) = len1;

            forceSpring = 1000*stateOpt(i,end)*...
                (vS.optims.len0Opt - len1);

            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            if strcmp(vS.params.model,'VPP')

                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                tanBetaLead(i) = (vS.params.rH + rVPP(i))*...
                    sinPsi/(len1 + (vS.params.rH +...
                    rVPP(i))*cosPsi);

                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end


            grf1y(i) = forceY;
            grf2y(i) = 0;

            grf1x(i) = forceX;
            grf2x(i) = 0;

        else

            fPos1 = vS.optims.fOptOrig(4*(cycle-1) + 5);

            len1 = sqrt((fPos1 - xHip)^2 + zHip^2);

            legLeft(i) = len1;

            forceSpring = 1000*stateOpt(i,end)*...
                (vS.optims.len0Opt - len1);

            sinLeg = (xHip - fPos1)/len1;
            cosLeg = zHip/len1;

            if strcmp(vS.params.model,'VPP')

                sinPsi = sinLeg*cos(stateOpt(i,3)) -...
                    cosLeg*sin(stateOpt(i,3));
                cosPsi = cosLeg*cos(stateOpt(i,3)) +...
                    sinLeg*sin(stateOpt(i,3));

                tanBetaLead(i) = (vS.params.rH + rVPP(end))*...
                    sinPsi/(len1 + (vS.params.rH +...
                    rVPP(end))*cosPsi);

                tauLead(i) = forceSpring*len1*tanBetaLead(i);

                forceX = sinLeg*forceSpring - cosLeg*tauLead(i)/len1;
                forceY = cosLeg*forceSpring + sinLeg*tauLead(i)/len1;

            else

                forceX = sinLeg*forceSpring;
                forceY = cosLeg*forceSpring;

            end


            grf1y(i) = forceY;
            grf2y(i) = 0;

            grf1x(i) = forceX;
            grf2x(i) = 0;

        end
        
        if (i/cycle) > 4*phaseDur
            
            cycle = cycle + 1;
            
        end

    end
    
    if ~strcmp(vS.params.size,'Full')

        vS.optims.grfVert = [grf1y', grf2y'];
        vS.optims.grfHor = [grf1x', grf2x'];
        vS.optims.grfNorm = sqrt(vS.optims.grfVert.^2 +...
            vS.optims.grfHor.^2);

        vS.optims.grfAng = atan2(vS.optims.grfHor, vS.optims.grfVert);
        indsFlip = find(vS.optims.grfAng >= pi/2);
        vS.optims.grfAng(indsFlip) = vS.optims.grfAng(indsFlip) - pi;
        
        vS.optims.legL = legLeft;
        vS.optims.legR = legRight;
        
        if strcmp(vS.params.model,'VPP')
            
            vS.optims.tauLag = tauLag;
            vS.optims.tauLead = tauLead;
            
            vS.optims.tanBetaLag = tanBetaLag;
            vS.optims.tanBetaLead = tanBetaLead;
            
        end

    else

        vS.optims.grfVertFull = [grf1y', grf2y'];
        vS.optims.grfHorFull = [grf1x', grf2x'];
        vS.optims.grfNormFull = sqrt(vS.optims.grfVertFull.^2 +...
            vS.optims.grfHorFull.^2);
        vS.optims.timeNormFull = timeNorm;

        vS.optims.grfAngFull = atan2(vS.optims.grfHorFull,...
            vS.optims.grfVertFull);
        indsFlip = find(vS.optims.grfAngFull >= pi/2);
        vS.optims.grfAngFull(indsFlip) =...
            vS.optims.grfAngFull(indsFlip) - pi;
        
        vS.optims.legLFull = legLeft;
        vS.optims.legRFull = legRight;
        
        if strcmp(vS.params.model,'VPP')
            
            vS.optims.tauLagFull = tauLag;
            vS.optims.tauLeadFull = tauLead;
            
            vS.optims.tanBetaLagFull = tanBetaLag;
            vS.optims.tanBetaLeadFull = tanBetaLead;
            
        end

    end

end
