%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intRK4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 2 September 2021
% Last Updated: 2 September 2021

% This function is used to numerically integrate the optimized trajectory
% of the BSLIP and VPP models

%% %%%%%%%%%%%%%%%%%%%%%%% Numerical Integration %%%%%%%%%%%%%%%%%%%%%%% %%

function vS = intRK4(vS)
    
    if isfield(vS.optims,'stateFull')
    
            vS.optims = rmfield(vS.optims,'stateFull');
    
    end
    
    vecLen = length(vS.optims.state(1,:));
    lenVec = length(vS.optims.state(:,1));
    phaseDur = vS.params.N;
    
    vS.optims.stateFull(1,1:vecLen) = vS.optims.state(1,:);
    
    cycle = 1;

    for i = 1:lenVec

        if (i < ((4*(cycle - 1) + 1)*phaseDur))

            pD = [vS.optims.stateFull(1,:), vS.optims.angOpt(1), vS.optims.kOptLagDot(i),...
                vS.optims.kOptLeadDot(i), vS.optims.fOpt(1), vS.optims.timeOpt(1)];
            phase = 1;

            for j=1:vS.params.M

                if strcmp(vS.params.model, 'VPP')
                    
                    k1D = intRK4VPP(vS, vS.optims.stateFull(end,:)',...
                        [pD, vS.optims.rVPPD(1)], vS.optims.uOpt(i), phase);
                    k2D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1D)/2),...
                        [pD, vS.optims.rVPPD(1)], vS.optims.uOpt(i), phase);
                    k3D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2D)/2),...
                        [pD, vS.optims.rVPPD(1)], vS.optims.uOpt(i), phase);
                    k4D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3D)),...
                        [pD, vS.optims.rVPPD(1)], vS.optims.uOpt(i), phase);
                    
                else
                    
                    k1D = intRK4BSLIP(vS, vS.optims.stateFull(end,:)',...
                        pD, vS.optims.uOpt(i), phase);
                    k2D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1D)/2),...
                        pD, vS.optims.uOpt(i), phase);
                    k3D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2D)/2),...
                        pD, vS.optims.uOpt(i), phase);
                    k4D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3D)),...
                        pD, vS.optims.uOpt(i), phase);
                
                end

                vS.optims.stateFull(end+1,:) =...
                    vS.optims.stateFull(end,:) +...
                    (vS.params.dt/6)*(k1D' + 2*k2D' + 2*k3D' + k4D');

            end

        elseif (i < ((4*(cycle - 1) + 2)*phaseDur))

            pS = [vS.optims.kOptLagDot(i); vS.optims.kOptLeadDot(i); vS.optims.fOpt(2); vS.optims.timeOpt(2)];
            phase = 2;

            for j=1:vS.params.M

                if strcmp(vS.params.model, 'VPP')
                    
                    k1S = intRK4VPP(vS, vS.optims.stateFull(end,:)',...
                        [pS; vS.optims.rVPPS(1)], [], phase);
                    k2S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1S)/2),...
                        [pS; vS.optims.rVPPS(1)], [], phase);
                    k3S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2S)/2),...
                        [pS; vS.optims.rVPPS(1)], [], phase);
                    k4S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3S)),...
                        [pS; vS.optims.rVPPS(1)], [], phase);
                    
                else
                    
                    k1S = intRK4BSLIP(vS, vS.optims.stateFull(end,:)',...
                        pS, [], phase);
                    k2S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1S)/2),...
                        pS, [], phase);
                    k3S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2S)/2),...
                        pS, [], phase);
                    k4S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3S)),...
                        pS, [], phase);
                    
                end

                vS.optims.stateFull(end+1,:) =...
                    vS.optims.stateFull(end,:) +...
                    (vS.params.dt/6)*(k1S' + 2*k2S' + 2*k3S' + k4S');

            end

        elseif (i < ((4*(cycle - 1) + 3)*phaseDur))

            pD = [vS.optims.stateFull((2*vS.params.N*vS.params.M + 1),:), vS.optims.angOpt(2),...
                vS.optims.kOptLagDot(i), vS.optims.kOptLeadDot(i),...
                vS.optims.fOpt(2), vS.optims.timeOpt(3)];
            phase = 1;

            for j=1:vS.params.M

                if strcmp(vS.params.model, 'VPP')
                    
                    k1D = intRK4VPP(vS, vS.optims.stateFull(end,:)',...
                        [pD, vS.optims.rVPPD(2)], vS.optims.uOpt(i-vS.params.N), phase);
                    k2D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1D)/2),...
                        [pD, vS.optims.rVPPD(2)], vS.optims.uOpt(i-vS.params.N), phase);
                    k3D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2D)/2),...
                        [pD, vS.optims.rVPPD(2)], vS.optims.uOpt(i-vS.params.N), phase);
                    k4D = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3D)),...
                        [pD, vS.optims.rVPPD(2)], vS.optims.uOpt(i-vS.params.N), phase);
                    
                else
                    
                    k1D = intRK4BSLIP(vS, vS.optims.stateFull(end,:)',...
                        pD, vS.optims.uOpt(i-vS.params.N), phase);
                    k2D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1D)/2),...
                        pD, vS.optims.uOpt(i-vS.params.N), phase);
                    k3D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2D)/2),...
                        pD, vS.optims.uOpt(i-vS.params.N), phase);
                    k4D = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3D)),...
                        pD, vS.optims.uOpt(i-vS.params.N), phase);
                    
                end

                vS.optims.stateFull(end+1,:) = vS.optims.stateFull(end,:) + (vS.params.dt/6)*(k1D' + 2*k2D' + 2*k3D' + k4D');

            end

        elseif (i < ((4*(cycle - 1) + 4)*phaseDur))

            pS = [vS.optims.kOptLagDot(i); vS.optims.kOptLeadDot(i); vS.optims.fOpt(3); vS.optims.timeOpt(4)];
            phase = 2;

            for j=1:vS.params.M

                if strcmp(vS.params.model, 'VPP')
                    
                    k1S = intRK4VPP(vS, vS.optims.stateFull(end,:)',...
                        [pS; vS.optims.rVPPS(2)], [], phase);
                    k2S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1S)/2),...
                        [pS; vS.optims.rVPPS(2)], [], phase);
                    k3S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2S)/2),...
                        [pS; vS.optims.rVPPS(2)], [], phase);
                    k4S = intRK4VPP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3S)),...
                        [pS; vS.optims.rVPPS(2)], [], phase);
                    
                else
                    
                    k1S = intRK4BSLIP(vS, vS.optims.stateFull(end,:)',...
                        pS, [], phase);
                    k2S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k1S)/2),...
                        pS, [], phase);
                    k3S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k2S)/2),...
                        pS, [], phase);
                    k4S = intRK4BSLIP(vS, (vS.optims.stateFull(end,:)' + (vS.params.dt*k3S)),...
                        pS, [], phase);
                    
                end

                vS.optims.stateFull(end+1,:) = vS.optims.stateFull(end,:) + (vS.params.dt/6)*(k1S' + 2*k2S' + 2*k3S' + k4S');

            end
            
        end
        
        if (i/cycle) > 4*phaseDur
            
            cycle = cycle + 1;
            
        end

    end

end
