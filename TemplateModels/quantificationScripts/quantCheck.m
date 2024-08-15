%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% quantCheck %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 19 May 2022
% Last Updated: 19 May 2022

% This function is used to calculate quantification checks for flagging
% whether the optimization process returned a satisfactory solution

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function vS = quantCheck(vS,dS)

    if vS.params.fitType == 1
        
        humDatHorSeg = vS.fit.posCOMSegHorF;
        humDatVertSeg = vS.fit.posCOMSegVertF;
        
    else
        
        humDatHorSeg = vS.fit.posCOMSegHor;
        humDatVertSeg = vS.fit.posCOMSegVert;
        
    end

    % Absolute Error for Center of Mass Trajectory (Segment Aligned)
    vS.quant.errorCOMSegHor = abs(humDatHorSeg -...
        vS.optims.state(:,1)');
    vS.quant.errorCOMSegVert = abs(humDatVertSeg -...
        vS.optims.state(:,2)');
    vS.quant.errorCOMSeg = vecnorm([humDatHorSeg;...
        humDatVertSeg] - vS.optims.state(:,1:2)',1);
    
    vS.quant.errorCOMSegMax = [max(abs(vS.quant.errorCOMSegHor)),...
        max(abs(vS.quant.errorCOMSegVert)), max(abs(vS.quant.errorCOMSeg))];
    vS.quant.errorCOMSegMaxND = vS.quant.errorCOMSegMax./vS.params.len0;
    
    % RMSE for Center of Mass Trajectory (Segment Aligned)
    vS.quant.rmseCOMSegHor = sqrt(mean(vS.quant.errorCOMSegHor.^2));
    vS.quant.rmseCOMSegVert = sqrt(mean(vS.quant.errorCOMSegVert.^2));
    vS.quant.rmseCOMSeg = sqrt(mean(vS.quant.errorCOMSeg.^2));
    
    vS.quant.rmseCOMSegHorND = vS.quant.rmseCOMSegHor./vS.params.len0;
    vS.quant.rmseCOMSegVertND = vS.quant.rmseCOMSegVert./vS.params.len0;
    vS.quant.rmseCOMSegND = vS.quant.rmseCOMSeg./vS.params.len0;
    
    % Absolute Error for Phase Durations
    vS.quant.errorPhase = abs(vS.optims.timeOpt - dS.data.timeFull);
    
    % RMSE Error for Phase Durations
    vS.quant.rmsePhase = sqrt(mean(vS.quant.errorPhase.^2));
    
    % Determine Time Aligned Human Data
    
    vS.fit.posCOMTimeHor = polyval(vS.fit.humPosHor, vS.optims.timeElem);
    vS.fit.posCOMTimeVert = polyval(vS.fit.humPosVert, vS.optims.timeElem);
    
    vS.fit.posCOMTimeHorF = fitFourier(vS.fit.humPosHorF,...
        vS.optims.timeElem, 'Hor');
    vS.fit.posCOMTimeVertF = fitFourier(vS.fit.humPosVertF,...
        vS.optims.timeElem, 'Vert');
    
    if vS.params.fitType == 1
        
        humDatHorTime = vS.fit.posCOMTimeHorF;
        humDatVertTime = vS.fit.posCOMTimeVertF;
        
    else
        
        humDatHorTime = vS.fit.posCOMTimeHor;
        humDatVertTime = vS.fit.posCOMTimeVert;
        
    end
    
    % Absolute Error for Center of Mass Trajectory (Time Instance Aligned)
    vS.quant.errorCOMTimeHor = abs(humDatHorTime -...
        vS.optims.state(:,1)');
    vS.quant.errorCOMTimeVert = abs(humDatVertTime -...
        vS.optims.state(:,2)');
    vS.quant.errorCOMTime = vecnorm([humDatHorTime;...
        humDatVertTime] - vS.optims.state(:,1:2)',1);
    
    vS.quant.errorCOMTimeMax = [max(abs(vS.quant.errorCOMTimeHor)),...
        max(abs(vS.quant.errorCOMTimeVert)), max(abs(vS.quant.errorCOMTime))];
    vS.quant.errorCOMTimeMaxND = vS.quant.errorCOMTimeMax./vS.params.len0;
    
    % RMSE for Center of Mass Trajectory (Time Instance Aligned)
    vS.quant.rmseCOMTimeHor = sqrt(mean(vS.quant.errorCOMTimeHor.^2));
    vS.quant.rmseCOMTimeVert = sqrt(mean(vS.quant.errorCOMTimeVert.^2));
    vS.quant.rmseCOMTime = sqrt(mean(vS.quant.errorCOMTime.^2));
    
    vS.quant.rmseCOMTimeHorND = vS.quant.rmseCOMTimeHor./vS.params.len0;
    vS.quant.rmseCOMTimeVertND = vS.quant.rmseCOMTimeVert./vS.params.len0;
    vS.quant.rmseCOMTimeND = vS.quant.rmseCOMTime./vS.params.len0;
    
    if isfield(dS.data,'grf')
        
        dS.interp.grfVert(find(dS.interp.grfVert(:,1)<0),1) = 0;
        dS.interp.grfVert(find(dS.interp.grfVert(:,2)<0),2) = 0;

        dS.interp.grfVertFull(find(dS.interp.grfVertFull(:,1)<0),1) = 0;
        dS.interp.grfVertFull(find(dS.interp.grfVertFull(:,2)<0),2) = 0;

        % Absolute Error for GRF Profiles
        vS.quant.errorGRFHor = abs(flip(dS.interp.grfHor,2) -...
            vS.optims.grfHor);
        vS.quant.errorGRFVert = abs(flip(dS.interp.grfVert,2) -...
            vS.optims.grfVert);

        vS.quant.errorGRFHorFull = abs(flip(dS.interp.grfHorFull,2) -...
            vS.optims.grfHorFull);
        vS.quant.errorGRFVertFull = abs(flip(dS.interp.grfVertFull,2) -...
            vS.optims.grfVertFull);

        vS.quant.errorGRFMax = [max(abs(vS.quant.errorGRFHor)),...
            max(abs(vS.quant.errorGRFHorFull));...
            max(abs(vS.quant.errorGRFVert)),...
            max(abs(vS.quant.errorGRFVertFull))];
        vS.quant.errorGRFMaxND = vS.quant.errorGRFMax./...
            (vS.params.m*vS.params.g);


        % RMSE for GRF Profile
        vS.quant.rmseGRFHor = sqrt(mean(vS.quant.errorGRFHor.^2,1));
        vS.quant.rmseGRFVert = sqrt(mean(vS.quant.errorGRFVert.^2,1));
        
        vS.quant.rmseGRFHorND = vS.quant.rmseGRFHor./...
            (vS.params.m*vS.params.g);
        vS.quant.rmseGRFVertND = vS.quant.rmseGRFVert./...
            (vS.params.m*vS.params.g);

        vS.quant.rmseGRFHorFull = sqrt(mean(vS.quant.errorGRFHorFull.^2,1));
        vS.quant.rmseGRFVertFull = sqrt(mean(vS.quant.errorGRFVertFull.^2,1));
        
        vS.quant.rmseGRFHorNDFull = vS.quant.rmseGRFHorFull./...
            (vS.params.m*vS.params.g);
        vS.quant.rmseGRFVertNDFull = vS.quant.rmseGRFVertFull./...
            (vS.params.m*vS.params.g);
        
    end
    
    saveData(vS, dS)

end