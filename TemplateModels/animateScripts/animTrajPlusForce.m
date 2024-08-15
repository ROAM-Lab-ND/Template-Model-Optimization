%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% animTraj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 2 September 2021
% Last Updated: 2 September 2021

% This function is used to generate the animated plot of the BSLIP or VPP
% model after the trajectory has been optimized to human data

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function vS = animTrajPlusForce(vS, saveVid)
    
    if strcmp(vS.params.model,'VPP')
       
        fileName = ['VPPCasADi' datestr(now,1) '.avi'];
        maxVPP = max(vS.optims.rVPP);        
        
    else
        
        fileName = ['BSLIPCasADi' datestr(now,1) '.avi'];
        maxVPP = 0;
        
    end
    
    if ~strcmp(vS.params.size,'Full')

        lenVec = length(vS.optims.state(:,1));
        timeVec = vS.optims.timeElem;

        optState = vS.optims.state;

        phaseDur = vS.params.N;
        
        scaleGRF = 1/max(max(vS.optims.grfNorm));
        optForceX = vS.optims.grfHor;
        optForceY = vS.optims.grfVert;

    else

        lenVec = length(vS.optims.stateFull(:,1));
        timeVec = vS.optims.timeFull;

        optState = vS.optims.stateFull;

        phaseDur = vS.params.N*vS.params.M;
        
        scaleGRF = 1/max(max(vS.optims.grfNormFull));
        optForceX = vS.optims.grfHorFull;
        optForceY = vS.optims.grfVertFull;

    end
    
    if strcmp(vS.params.model,'VPP')
           
        rVPP = repelem(vS.optims.rVPP(:), phaseDur);
        rVPP = [rVPP; vS.optims.rVPP(1)];

    end
    
    if saveVid
        
        % Open video file for recording
        writerObj = VideoWriter(fileName);
        writerObj.FrameRate = 10;
        writerObj.Quality = 100;
        open(writerObj);
        
    end

    % Initialize animate figure
    figTraj = figure;
    
    axesTraj = axes(figTraj);
    
    title(axesTraj,'COM Trajectory')
    xlabel(axesTraj,'Fore/Aft COM Position [m]')
    ylabel(axesTraj,'Vertical COM Position [m]')
    
    hold on

    axis(axesTraj, [vS.optims.fOpt(1) vS.optims.state(end,1)+0.6...
        0 max(vS.optims.state(:,2))+maxVPP+0.1])
    
    if saveVid

                frame = getframe(gcf);
                writeVideo(writerObj,frame);

    end
    
    humDatVert = polyval(vS.fit.humPosVert,timeVec);
    humDatHor = polyval(vS.fit.humPosHor,timeVec);
    
    animHum = animatedline(humDatHor, humDatVert, 'LineStyle', '-.',...
        'Color', 'm', 'LineWidth', 3);
    animModel = animatedline(vS.optims.state(1,1),vS.optims.state(1,2),...
        'Color', 'k', 'LineWidth',3);
    
    plot(axesTraj, NaN, NaN, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    plot(axesTraj, NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    plot(axesTraj, NaN, NaN, 'b', 'LineWidth', 2);
    plot(axesTraj, NaN, NaN, 'r', 'LineWidth', 2);
    
    legend('Human COM Trajectory', [vS.params.model ' COM Trajectory'],...
        'Touchdown Event', 'Liftoff Event',...
        'Right Leg','Left Leg', 'Location', 'NorthEast');
    
    for i = 1:lenVec
        
       addpoints(animModel, optState(i,1), optState(i,2));
       
       if strcmp(vS.params.model,'VPP')

            xHip = optState(i,1) - vS.params.rH*sin(optState(i,3));
            yHip = optState(i,2) - vS.params.rH*cos(optState(i,3));

            xTrunk = xHip + (vS.params.rH + maxVPP + 0.1)*sin(optState(i,3));
            yTrunk = yHip + (vS.params.rH + maxVPP + 0.1)*cos(optState(i,3));

            plotVPP = [(optState(i,1) + rVPP(i)*sin(optState(i,3)));...
                (optState(i,2) + rVPP(i)*cos(optState(i,3)))];

       else

            xHip = optState(i,1);
            yHip = optState(i,2);

       end
       
       switch i
           
           case {1, (4*phaseDur + 1)}
               
               styleLag = 'r';
               styleLead = 'b--';
               
           case {(phaseDur + 1)}
               
               styleLag = 'r:';
               styleLead = 'b';
               
           case {(2*phaseDur + 1)}
               
               styleLag = 'b';
               styleLead = 'r--';
               
           case {(3*phaseDur + 1)}
               
               styleLag = 'b:';
               styleLead = 'r';
               
       end
       
       switch i
           
           case {1, (phaseDur + 1)}
               
               footLag = vS.optims.fOpt(1);
               footLead = vS.optims.fOpt(2);

           case {(2*phaseDur + 1), (3*phaseDur + 1)}
               
               footLag = vS.optims.fOpt(2);
               footLead = vS.optims.fOpt(3);
             
           case {(4*phaseDur + 1)}
               
               footLag = vS.optims.fOpt(3);
               footLead = optState(i,1) +...
                   optState(i,2)*tan(vS.optims.angOpt(end));            
               
       end
       
       if ((i > 2*phaseDur) && (i <= 4*phaseDur))
           
           forceLag = [optForceX(i,2), optForceY(i,2)];
           forceLead = [optForceX(i,1), optForceY(i,1)];
           
       else
           
           forceLag = [optForceX(i,1), optForceY(i,1)];
           forceLead = [optForceX(i,2), optForceY(i,2)];
           
       end
           
           
       
       if ((i > phaseDur) && (i <=  2*phaseDur)) ||...
               ((i > 3*phaseDur) && (i <= 4*phaseDur))
           
           plotPosLead = plot([footLead, xHip], [0, yHip],...
               styleLead, 'LineWidth', 2, 'HandleVisibility', 'off');
           plotForceLead = quiver(footLead, 0,...
               forceLead(1), forceLead(2), scaleGRF,...
               'Color', 'k', 'LineWidth', 2);
               
       else
           
           plotPosLag =  plot([footLag, xHip], [0, yHip],...
               styleLag, 'LineWidth', 2, 'HandleVisibility', 'off');
           plotPosLead = plot([footLead, xHip], [0, yHip],...
               styleLead, 'LineWidth', 2, 'HandleVisibility', 'off');
           plotForceLag = quiver(footLag, 0,...
               forceLag(1), forceLag(2), scaleGRF,...
               'Color', 'k', 'LineWidth', 2);
           plotForceLead = quiver(footLead, 0,...
               forceLead(1), forceLead(2), scaleGRF,...
               'Color', 'k', 'LineWidth', 2);
           
       end
       
       if strcmp(vS.params.model,'VPP')
           
           plotAngleVPP = plot([xHip,xTrunk],[yHip,yTrunk],...
               'g', 'LineWidth', 2, 'HandleVisibility', 'off');
           plotPosVPP = plot(plotVPP(1), plotVPP(2),...
               'r.', 'MarkerSize', 15, 'HandleVisibility', 'off');

       end
       
       drawnow
    
       if saveVid
           
           frame = getframe(gcf);
           writeVideo(writerObj,frame);
           
       end
       
       pause(0.02);
       
       if (i < (4*phaseDur + 1))
       
           if exist('plotPosLag', 'var') == 1

               delete(plotPosLag)

           end

           if exist('plotPosLead', 'var') == 1

               delete(plotPosLead)

           end

           if exist('plotAngleVPP', 'var') == 1

               delete(plotAngleVPP)

           end

           if exist('plotPosVPP', 'var') == 1

               delete(plotPosVPP)

           end

           if exist('plotForceLag', 'var') == 1

               delete(plotForceLag)

           end

           if exist('plotForceLead', 'var') == 1

               delete(plotForceLead)

           end
           
       end
        
    end
    
    hold off
    
    if saveVid

        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        
        close(writerObj);
        
    end
    
end