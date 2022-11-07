%% %%%%%%%%%%%%%%%%%%%%%%%%%%% animVelTime %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 17 November 2021
% Last Updated: 17 November 2021

% This function is used to compare the optimized velocity for the BSLIP
% or VPP model with the parameterized human velocity

%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
% NOTE: Mocap uses derivative of the polyfit to the human trajectory for
% the parameterized human velocity
%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function vS = animVelTime(vS, saveVid)
    
    if strcmp(vS.params.model,'VPP')
       
        fileNameJoint = ['VPPTemporalVel' datestr(now,1) '.avi'];     
        
    else
        
        fileNameJoint = ['BSLIPTemporalVel' datestr(now,1) '.avi'];
        
    end
    
    if ~strcmp(vS.params.size, 'Full')
        
        phaseDur = vS.params.N;
        timeOpt = vS.optims.timeElem;
        timeHum = vS.fit.tHum;
        
        velXOpt = vS.optims.state(:,(length(vS.optims.state(1,:))/2));
        velYOpt = vS.optims.state(:,(length(vS.optims.state(1,:))/2 + 1));
        
    else
        
        phaseDur = vS.params.N*vS.params.M;
        timeOpt = vS.optims.timeFull;
        timeHum = vS.fit.tHumFull;
        
        velXOpt = vS.optims.stateFull(:,...
            (length(vS.optims.stateFull(1,:))/2));
        velYOpt = vS.optims.stateFull(:,...
            (length(vS.optims.stateFull(1,:))/2 + 1));
        
    end
    
    if saveVid
        
        % Open video file for recording
        writerObjJoint = VideoWriter(fileNameJoint);
        writerObjJoint.FrameRate = 10;
        writerObjJoint.Quality = 100;
        open(writerObjJoint);
        
    end
    
    lenVec = length(velXOpt);
        
    humDatVert = polyval(vS.fit.humVelVert, timeHum);
    humDatHor = polyval(vS.fit.humVelHor, timeHum);

    % Initialize animate figure
    jointFig = figure;
    
    axesX = subplot(2,1,1);
    
    title(axesX,'COM Fore/Aft Velocity')
    xlabel(axesX,'Time [s]')
    ylabel(axesX,'Fore/Aft COM Velocity [m/s]')
    
    axis(axesX, [-0.02...
        (max(timeOpt(end), timeHum(end)) + 0.02)...
        (min(min(velXOpt), min(humDatHor)) - 0.02)...
        (max(max(velXOpt), max(humDatHor)) + 0.02)]);
    
    animModelXOpt = animatedline(axesX, timeOpt(1),...
        velXOpt(1), 'Color', 'r', 'Linewidth', 1);
    
    animModelXHum = animatedline(axesX, timeHum(1),...
        humDatHor(1), 'Color', 'k', 'Linewidth', 1);
    
    animContactXOpt = animatedline(axesX, timeOpt(1),...
        velXOpt(1), 'Color', 'g', 'Marker', '*',...
        'MarkerSize', 15, 'LineStyle', 'none');
    
    animContactXHum = animatedline(axesX, timeHum(1),...
        humDatHor(1), 'Color', 'b', 'Marker', '+',...
        'MarkerSize', 15, 'LineStyle', 'none');
    
    legend(axesX, 'Optimized Trajectory', 'Experimental Trajectory',...
        'Optimized Contact Instant', 'Experimental Contact Instant',...
        'Location', 'southeast')

    axesY = subplot(2,1,2);

    title(axesY,'COM Vertical Trajectory')
    xlabel(axesY,'Time [s]')
    ylabel(axesY,'Vertical COM Position [m]')
    
    axis(axesY, [-0.02...
        (max(timeOpt(end), timeHum(end)) + 0.02)...
        (min(min(velYOpt), min(humDatVert)) - 0.02)...
        (max(max(velYOpt), max(humDatVert)) + 0.02)])
    
    animModelYOpt = animatedline(axesY, timeOpt(1),...
        velYOpt(1), 'Color', 'r', 'Linewidth', 1);
    
    animModelYHum = animatedline(axesY, timeHum(1),...
        humDatVert(1), 'Color', 'k', 'Linewidth', 1);
    
    animContactYOpt = animatedline(axesY, timeOpt(1),...
        velYOpt(1), 'Color', 'g', 'Marker', '*',...
        'MarkerSize', 15, 'LineStyle', 'none');
    
    animContactYHum = animatedline(axesY, timeHum(1),...
        humDatVert(1), 'Color', 'b', 'Marker', '+',...
        'MarkerSize', 15, 'LineStyle', 'none');
    
    legend(axesY, 'Optimized Trajectory', 'Experimental Trajectory',...
        'Optimized Contact Instant', 'Experimental Contact Instant',...
        'Location', 'southeast')
    
    if saveVid

            frameJoint = getframe(jointFig);
            writeVideo(writerObjJoint, frameJoint);

    end
    
    for i = 1:lenVec
        
        addpoints(animModelXOpt, timeOpt(i), velXOpt(i));
        addpoints(animModelXHum, timeHum(i), humDatHor(i));

        addpoints(animModelYOpt, timeOpt(i), velYOpt(i));
        addpoints(animModelYHum, timeHum(i), humDatVert(i));
        
        if (mod(i, phaseDur) == 1)
                
            addpoints(animContactXOpt, timeOpt(i), velXOpt(i));
            addpoints(animContactXHum, timeHum(i), humDatHor(i));

            addpoints(animContactYOpt, timeOpt(i), velYOpt(i));
            addpoints(animContactYHum, timeHum(i), humDatVert(i));

        end

        drawnow

        if saveVid

            frameJoint = getframe(jointFig);
            writeVideo(writerObjJoint, frameJoint);

        end

        pause(0.05);

    end
    
    hold off
    
    if saveVid
        
        close(writerObjJoint);
        
    end
    
end