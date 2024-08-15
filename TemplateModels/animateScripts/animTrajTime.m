%% %%%%%%%%%%%%%%%%%%%%%%%%%%% animTrajTime %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 11 November 2021
% Last Updated: 13 January 2022

% This function is used to compare the optimized COM trajectory for the BSLIP
% or VPP model with the parameterized human COM trajectory

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function vS = animTrajTime(vS, saveVid)
    
    if strcmp(vS.params.model,'VPP')
       
        fileNameJoint = ['VPPTemporalCOM' vS.params.springType...
            datestr(now,1) '.avi'];        
        
    else
        
        fileNameJoint = ['BSLIPTemporalCOM' vS.params.springType...
            datestr(now,1) '.avi'];
        
    end
    
    if saveVid
        
        % Open video file for recording
        writerObjJoint = VideoWriter(fileNameJoint);
        writerObjJoint.FrameRate = 10;
        writerObjJoint.Quality = 100;
        open(writerObjJoint);
        
    end
    
    if ~strcmp(vS.params.size,'Full')
        
        timeOpt = vS.optims.timeElem;
        timeHum = vS.fit.tHum;
        stateOpt = vS.optims.state;
        phaseDur = vS.params.N;
        
    else
        
        timeOpt = vS.optims.timeFull;
        timeHum = vS.fit.tHumFull;
        stateOpt = vS.optims.stateFull;
        phaseDur = vS.params.N*vS.params.M;
        
    end
    
    lenVec = length(stateOpt(:,1));
    
    if vS.params.fitType == 1
    
        humDatVert = fitFourier(vS.fit.humPosVertF, timeHum, 'Vert');
        humDatHor = fitFourier(vS.fit.humPosHorF, timeHum, 'Hor');
        
    else
        
        humDatVert = polyval(vS.fit.humPosVert, timeHum);
        humDatHor = polyval(vS.fit.humPosHor, timeHum);
        
    end

    % Initialize animate figure
    jointFig = figure('Position', [10 10 900 600]);
    
    axesX = subplot(2,1,1);
    
    title(axesX,'COM Fore/Aft Trajectory')
    xlabel(axesX,'Time [s]')
    ylabel(axesX,'Fore/Aft COM Position [m]')
    
    axis(axesX, [-0.02 (max(max(timeOpt), max(timeHum)) + 0.02)...
        -0.02 (max(stateOpt(end,1), humDatHor(end)) + 0.02)]);
    
    animModelXOpt = animatedline(axesX, timeOpt(1),...
        stateOpt(1,1), 'Color', 'r', 'LineStyle', '-',...
        'LineWidth', 2);
    
    animModelXHum = animatedline(axesX, timeHum(1),...
        humDatHor(1), 'Color', 'k', 'LineStyle', '--',...
        'LineWidth', 2);
    
    animContactXOpt = animatedline(axesX, timeOpt(1),...
        stateOpt(1,1), 'Color', 'g', 'Marker', '*',...
        'MarkerSize', 15, 'LineStyle', 'None');
    
    animContactXHum = animatedline(axesX, timeHum(1),...
        humDatHor(1), 'Color', 'b', 'Marker', '+',...
        'MarkerSize', 15, 'LineStyle', 'None');
    
    legend(axesX, 'Optimized Trajectory', 'Experimental Trajectory',...
        'Optimized Contact Instant', 'Experimental Contact Instant',...
        'Location', 'southeast')

    axesY = subplot(2,1,2);

    title(axesY, ['COM Vertical Trajectory ['...
        vS.params.model ' ' vS.params.springType ' ' 'Spring]'])
    xlabel(axesY,'Time [s]')
    ylabel(axesY,'Vertical COM Position [m]')
    
    axis(axesY, [-0.02 (max(max(timeOpt), max(timeHum)) + 0.02)...
        (min(min(stateOpt(:,2)), min(humDatVert)) - 0.02)...
        (max(max(stateOpt(:,2)), max(humDatVert)) + 0.02)])

%     axis(axesY, [-0.02 (max(max(timeOpt), max(timeHum)) + 0.02)...
%             0.98 1.08])
    yticks(axesY, [0.98 1 1.02 1.04 1.06 1.08])
    set(gca, 'Fontsize', 14)
    
    animModelYOpt = animatedline(axesY, timeOpt(1),...
        stateOpt(1,2), 'Color', 'r', 'LineStyle', '-',...
        'Linewidth', 2);
    
    animModelYHum = animatedline(axesY, timeHum(1),...
        humDatVert(1), 'Color', 'k', 'LineStyle', '--',...
        'Linewidth', 2);
    
    animContactYHum = animatedline(axesY, timeHum(1),...
        humDatVert(1), 'Color', [0.15 0.15 0.15], 'Marker', '.',...
        'MarkerSize', 25, 'LineStyle', 'none');
    
    animContactYOpt = animatedline(axesY, timeOpt(1),...
        stateOpt(1,2), 'Color', [0.6 0 0], 'Marker', '*',...
        'MarkerSize', 25, 'LineStyle', 'none');
    
    legend(axesY, [animModelYOpt, animContactYOpt,...
            animModelYHum, animContactYHum],...
            {'COM (Opt)', 'Phase Switch (Opt)', 'COM (Exp)',...
            'Phase Switch (Exp)'}, 'Location', 'south', 'NumColumns', 2);
    
    if saveVid

        frameJoint = getframe(jointFig);
        writeVideo(writerObjJoint, frameJoint);

    end
    
    for i = 1:lenVec
        
        addpoints(animModelXHum, timeHum(i), humDatHor(i));        
        addpoints(animModelXOpt, timeOpt(i), stateOpt(i,1));
        
        addpoints(animModelYHum, timeHum(i), humDatVert(i));
        addpoints(animModelYOpt, timeOpt(i), stateOpt(i,2));
        
        if (mod(i, phaseDur) == 1)
            
            addpoints(animContactXHum, timeHum(i), humDatHor(i));
            addpoints(animContactXOpt, timeOpt(i), stateOpt(i,1));

            addpoints(animContactYHum, timeHum(i), humDatVert(i));
            addpoints(animContactYOpt, timeOpt(i), stateOpt(i,2));           
            
        end
        
        drawnow
            
        if saveVid
            
            frameJoint = getframe(jointFig);
            writeVideo(writerObjJoint, frameJoint);
            
        end
        
%         pause(0.01);
        
    end
    
    if saveVid
        
        close(writerObjJoint);
        
    end
    
end