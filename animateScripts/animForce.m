%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% animForce %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 September 2021
% Last Updated: 9 September 2021

% This function is used to generate the animated plot of the BSLIP or VPP
% ground reaction forces after the trajectory has been optimized to human
% data

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage
%   saveVid - Flag to save resulting animation

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [varS, datS] = animForce(varS, datS, saveVid)
    
    % Check if animating full or base state values
    if ~strcmp(varS.params.size, 'Full')
        
        % Calculate length of state vector
        lenVec = length(varS.optims.state(:,1));
        
        % Calculate normalized time
        timeNorm = varS.optims.timeElem/max(varS.optims.timeElem);
        
        % Grab vertical and horizontal subject GRFs for easier coding
        humDatVert = datS.interp.grfVert;
        humDatHor = datS.interp.grfHor;
        
        % Grab vertical and horizontal optimized GRFs for easier coding
        optDatVert = varS.optims.grfVert;
        optDatHor = varS.optims.grfHor;

    else
        
        % Calculate length of state vector
        lenVec = length(varS.optims.stateFull(:,1));
        
        % Calculate normalized time
        timeNorm = varS.optims.timeFull/max(varS.optims.timeFull);
        
        % Grab vertical and horizontal subject GRFs for easier coding
        humDatVert = datS.interp.grfVertFull;
        humDatHor = datS.interp.grfHorFull;
        
        % Grab vertical and horizontal optimized GRFs for easier coding
        optDatVert = varS.optims.grfVertFull;
        optDatHor = varS.optims.grfHorFull;
        
    end
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')
       
        % Build filename to store animation in
        fileName = ['VPPCasADiGRF' varS.params.springType...
            datestr(now,1) '.avi'];        
        
    else
        
        % Build filename to store animation in
        fileName = ['BSLIPCasADiGRF' varS.params.springType...
            datestr(now,1) '.avi'];
        
    end
    
    % Check if save video flag is set
    if saveVid
        
        % Create video object
        writerObj = VideoWriter(fileName);
        
        % Configure video file options
        writerObj.FrameRate = 10;
        writerObj.Quality = 100;
        
        % Open video object
        open(writerObj);
        
    end

    % 
    bodyWeight = varS.params.m*varS.params.g*ones(1,lenVec);
    
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
    
    animForceHorOptR = animatedline(axesX, varS.optims.timeElem(1),...
        optDatHor(1,2), 'Color', 'b', 'LineStyle', '-', 'Linewidth', 2);
    
    if isfield(datS.data,'grf')
    
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
        varS.params.model ' ' varS.params.springType ' ' 'Spring]'])
    xlabel(axesY,'Time [Normalized]')
    ylabel(axesY,'Force [N]')
    
%     axis(axesY, [0 1,...
%         (min(min(min(humDatVert),min(optDatVert))) - 50),...
%         (max(max(max(humDatVert),max(optDatVert))) + 500)]);
    
    axis(axesY, [0 1,...
        (min(min(min(humDatVert),min(optDatVert))) - 50),...
        1200]);
    set(gca, 'Fontsize', 14)
    
    animForceBW = animatedline(axesY, timeNorm, bodyWeight,...
        'Color', 'k', 'LineStyle', '-.', 'Linewidth', 1);
    
    animForceVertOptL = animatedline(axesY, timeNorm(1),...
        optDatVert(1,1), 'Color', 'r', 'LineStyle', '-', 'Linewidth', 2);
    
    animForceVertOptR = animatedline(axesY, varS.optims.timeElem(1),...
        optDatVert(1,2), 'Color', 'r', 'LineStyle', '--', 'Linewidth', 2);
    
    if isfield(datS.data,'grf')
    
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