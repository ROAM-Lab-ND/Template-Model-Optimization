%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% animForce %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 September 2021
% Last Updated: 9 September 2021

% This function is used to generate the animated plot of the BSLIP or VPP
% ground reaction forces after the trajectory has been optimized to human
% data

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [vS, dS] = animForce(vS, dS, saveVid)
    
    if ~strcmp(vS.params.size,'Full')
        
        lenVec = length(vS.optims.state(:,1));
        timeNorm = vS.optims.timeElem/max(vS.optims.timeElem);
        
        humDatVert = dS.interp.grfVert;
        humDatHor = dS.interp.grfHor;
        
        optDatVert = vS.optims.grfVert;
        optDatHor = vS.optims.grfHor;

    else
        
        lenVec = length(vS.optims.stateFull(:,1));
        timeNorm = vS.optims.timeFull/max(vS.optims.timeFull);
        
        humDatVert = dS.interp.grfVertFull;
        humDatHor = dS.interp.grfHorFull;
        
        optDatVert = vS.optims.grfVertFull;
        optDatHor = vS.optims.grfHorFull;
        
    end
    
    if strcmp(vS.params.model,'VPP')
       
        fileName = ['VPPCasADiGRF' vS.params.springType...
            datestr(now,1) '.avi'];
        
        
    else
        
        fileName = ['BSLIPCasADiGRF' vS.params.springType...
            datestr(now,1) '.avi'];
        
    end
    
    if saveVid
        
        % Open video file for recording
        writerObj = VideoWriter(fileName);
        writerObj.FrameRate = 10;
        writerObj.Quality = 100;
        open(writerObj);
        
    end

    bodyWeight = vS.params.m*vS.params.g*ones(1,lenVec);
    
        % Initialize animate figure
    figX = figure;%('Position', [10 10 900 300]);
    
    axesX = axes(figX);
    
    title(axesX,'Fore/Aft GRF')
    xlabel(axesX,'Time [Normalized]')
    ylabel(axesX,'Force [N]')
    
    axis(axesX, [0 1,...
        (min(min(min(humDatHor),min(optDatHor))) - 100),...
        (max(max(max(humDatHor),max(optDatHor))) + 100)]);
    
    animForceHorOptL = animatedline(axesX, timeNorm(1),...
        optDatHor(1,1), 'Color', 'r', 'LineStyle', '-', 'Linewidth', 2);
    
    animForceHorOptR = animatedline(axesX, vS.optims.timeElem(1),...
        optDatHor(1,2), 'Color', 'b', 'LineStyle', '-', 'Linewidth', 2);
    
    if isfield(dS.data,'grf')
    
        animForceHorHumL = animatedline(axesX, timeNorm(1), humDatHor(1,1),...
            'Color', 'm', 'LineStyle', '--', 'Linewidth', 2);

        animForceHorHumR = animatedline(axesX, timeNorm(1), humDatHor(1,2),...
            'Color', 'g', 'LineStyle', '--', 'Linewidth', 2);

        legend(axesX, 'Left Leg (Opt)', 'Right Leg (Opt)',...
            'Left Leg (Subject)', 'Right Leg (Subject)',...
            'Location', 'northeast')
        
    else
        
        legend(axesX, 'Left Leg (Opt)', 'Right Leg (Opt)',...
            'Location', 'northeast')
        
    end
    
    figY = figure;%('Position', [10 10 500 300]);

    axesY = axes(figY);

    title(axesY, ['Vertical GRF ['...
        vS.params.model ' ' vS.params.springType ' ' 'Spring]'])
    xlabel(axesY,'Time [Normalized]')
    ylabel(axesY,'Force [N]')
    
    axis(axesY, [0 1,...
        (min(min(min(humDatVert),min(optDatVert))) - 50),...
        (max(max(max(humDatVert),max(optDatVert))) + 500)]);
    
%     axis(axesY, [0 1,...
%         (min(min(min(humDatVert),min(optDatVert))) - 50),...
%         1200]);
    set(gca, 'Fontsize', 14)
    
    animForceBW = animatedline(axesY, timeNorm, bodyWeight,...
        'Color', 'k', 'LineStyle', '-.', 'Linewidth', 1);
    
    animForceVertOptL = animatedline(axesY, timeNorm(1),...
        optDatVert(1,1), 'Color', 'r', 'LineStyle', '-', 'Linewidth', 2);
    
    animForceVertOptR = animatedline(axesY, vS.optims.timeElem(1),...
        optDatVert(1,2), 'Color', 'r', 'LineStyle', '--', 'Linewidth', 2);
    
    if isfield(dS.data,'grf')
    
        animForceVertHumL = animatedline(axesY, timeNorm(1), humDatVert(1,1),...
            'Color', 'k', 'LineStyle', '--', 'Linewidth', 2);

        animForceVertHumR = animatedline(axesY, timeNorm(1), humDatVert(1,2),...
            'Color', 'k', 'LineStyle', '-', 'Linewidth', 2);

        legend(axesY, [animForceVertOptR, animForceVertOptL, animForceBW,...
            animForceVertHumL, animForceVertHumR],...
            {'Right Leg (Opt)', 'Left Leg (Opt)', 'Bodyweight',...
            'Right Leg (Exp)', 'Left Leg (Exp)'}, 'Location', 'northeast',...
            'NumColumns', 2);
        
    else
        
        legend(axesY, 'Left Leg (Opt)', 'Right Leg (Opt)', 'Body Weight',...
            'Location', 'northeast')
        
    end    

    if saveVid

        frame = getframe(gcf);
        writeVideo(writerObj,frame);

    end
    
    for i=1:lenVec
        
        addpoints(animForceHorOptL, timeNorm(i), optDatHor(i,1));
        addpoints(animForceHorOptR, timeNorm(i), optDatHor(i,2));
        
        addpoints(animForceHorHumL, timeNorm(i), humDatHor(i,1));
        addpoints(animForceHorHumR, timeNorm(i), humDatHor(i,2));
        
        addpoints(animForceVertOptL, timeNorm(i), optDatVert(i,1));
        addpoints(animForceVertOptR, timeNorm(i), optDatVert(i,2));
        
        addpoints(animForceVertHumL, timeNorm(i), humDatVert(i,1));
        addpoints(animForceVertHumR, timeNorm(i), humDatVert(i,2));
        
        drawnow;
        
        if saveVid

            frame = getframe(gcf);
            writeVideo(writerObj,frame);

        end
        
    end
    
%     hold on
%     plot(axesY, timeNorm, bodyWeight, '--k')
    
    if saveVid
        
        close(writerObj);
        
    end
    
end