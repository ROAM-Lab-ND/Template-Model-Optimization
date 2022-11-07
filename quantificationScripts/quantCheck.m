%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% quantCheck %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 19 May 2022
% Last Updated: 19 May 2022

% This function is used to calculate quantification checks for flagging
% whether the optimization process returned a satisfactory solution

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function varS = quantCheck(varS, datS)

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Center of Mass %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if using polyfit or Fourier fit
    if varS.params.fitType == 1
       
        % Grab subject CoM position
        humDatHorSeg = varS.fit.posCOMSegHorF;
        humDatVertSeg = varS.fit.posCOMSegVertF;
        
    else
        
        % Grab subject CoM position
        humDatHorSeg = varS.fit.posCOMSegHor;
        humDatVertSeg = varS.fit.posCOMSegVert;
        
    end
    
    % Calculate segment by segment CoM tracking error
    varS.quant.errorCOMSegHor = abs(humDatHorSeg -...
        varS.optims.state(:,1)');
    
    varS.quant.errorCOMSegVert = abs(humDatVertSeg -...
        varS.optims.state(:,2)');
    
    varS.quant.errorCOMSeg = vecnorm([humDatHorSeg;...
        humDatVertSeg] - varS.optims.state(:,1:2)',1);
    
    % Grab maximum segment by segment CoM tracking error
    varS.quant.errorCOMSegMax = [max(abs(varS.quant.errorCOMSegHor)),...
        max(abs(varS.quant.errorCOMSegVert)),...
        max(abs(varS.quant.errorCOMSeg))];
    varS.quant.errorCOMSegMaxND =...
        varS.quant.errorCOMSegMax./varS.params.len0;
    
    % Calculate RMSE for segment by segment CoM tracking error
    varS.quant.rmseCOMSegHor = sqrt(mean(varS.quant.errorCOMSegHor.^2));
    varS.quant.rmseCOMSegVert = sqrt(mean(varS.quant.errorCOMSegVert.^2));
    varS.quant.rmseCOMSeg = sqrt(mean(varS.quant.errorCOMSeg.^2));
    
    varS.quant.rmseCOMSegHorND =...
        varS.quant.rmseCOMSegHor./varS.params.len0;
    varS.quant.rmseCOMSegVertND =...
        varS.quant.rmseCOMSegVert./varS.params.len0;
    varS.quant.rmseCOMSegND =...
        varS.quant.rmseCOMSeg./varS.params.len0;
    
    % Evaluate polyfit and Fourier fit for subject CoM position   
    varS.fit.posCOMTimeHor = polyval(varS.fit.humPosHor,...
        varS.optims.timeElem);
    varS.fit.posCOMTimeVert = polyval(varS.fit.humPosVert,...
        varS.optims.timeElem);
    
    varS.fit.posCOMTimeHorF = fitFourier(varS.fit.humPosHorF,...
        varS.optims.timeElem, 'Hor');
    varS.fit.posCOMTimeVertF = fitFourier(varS.fit.humPosVertF,...
        varS.optims.timeElem, 'Vert');
    
    % Check if using polyfit or Fourier fit
    if varS.params.fitType == 1
        
        % Grab subject CoM position
        humDatHorTime = varS.fit.posCOMTimeHorF;
        humDatVertTime = varS.fit.posCOMTimeVertF;
        
    else
        
        % Grab subject CoM position
        humDatHorTime = varS.fit.posCOMTimeHor;
        humDatVertTime = varS.fit.posCOMTimeVert;
        
    end
    
    % Calculate temporal CoM tracking error
    varS.quant.errorCOMTimeHor = abs(humDatHorTime -...
        varS.optims.state(:,1)');
    
    varS.quant.errorCOMTimeVert = abs(humDatVertTime -...
        varS.optims.state(:,2)');
    
    varS.quant.errorCOMTime = vecnorm([humDatHorTime;...
        humDatVertTime] - varS.optims.state(:,1:2)',1);
    
    % Calculate maximum temporal CoM tracking error
    varS.quant.errorCOMTimeMax = [max(abs(varS.quant.errorCOMTimeHor)),...
        max(abs(varS.quant.errorCOMTimeVert)),...
        max(abs(varS.quant.errorCOMTime))];
    varS.quant.errorCOMTimeMaxND =...
        varS.quant.errorCOMTimeMax./varS.params.len0;
    
    % Calculate RMSE for temporal CoM tracking error
    varS.quant.rmseCOMTimeHor = sqrt(mean(varS.quant.errorCOMTimeHor.^2));
    varS.quant.rmseCOMTimeVert =...
        sqrt(mean(varS.quant.errorCOMTimeVert.^2));
    varS.quant.rmseCOMTime = sqrt(mean(varS.quant.errorCOMTime.^2));
    
    varS.quant.rmseCOMTimeHorND =...
        varS.quant.rmseCOMTimeHor./varS.params.len0;
    varS.quant.rmseCOMTimeVertND =...
        varS.quant.rmseCOMTimeVert./varS.params.len0;
    varS.quant.rmseCOMTimeND =...
        varS.quant.rmseCOMTime./varS.params.len0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Phase Duration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate phase by phase gait event tracking error
    varS.quant.errorPhase = abs(varS.optims.timeOpt - datS.data.timeFull);
    
    % Calculate RMSE phase by phase gait event tracking error
    varS.quant.rmsePhase = sqrt(mean(varS.quant.errorPhase.^2));
    
    %%%%%%%%%%%%%%%%%%%%%%%% Ground Reaction Force %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if GRF data exists for subject
    if isfield(datS.data,'grf')
        
        % Correct vertical GRF values below 0 to 0
        datS.interp.grfVert(find(datS.interp.grfVert(:,1)<0),1) = 0;
        datS.interp.grfVert(find(datS.interp.grfVert(:,2)<0),2) = 0;

        datS.interp.grfVertFull(...
            find(datS.interp.grfVertFull(:,1)<0),1) = 0;
        datS.interp.grfVertFull(...
            find(datS.interp.grfVertFull(:,2)<0),2) = 0;

        % Calculate GRF tracking error
        varS.quant.errorGRFHor = abs(flip(datS.interp.grfHor,2) -...
            varS.optims.grfHor);
        varS.quant.errorGRFVert = abs(flip(datS.interp.grfVert,2) -...
            varS.optims.grfVert);

        varS.quant.errorGRFHorFull = abs(flip(datS.interp.grfHorFull,2) -...
            varS.optims.grfHorFull);
        varS.quant.errorGRFVertFull = abs(flip(datS.interp.grfVertFull,2) -...
            varS.optims.grfVertFull);

        % Calculate maximum GRF tracking error
        varS.quant.errorGRFMax = [max(abs(varS.quant.errorGRFHor)),...
            max(abs(varS.quant.errorGRFHorFull));...
            max(abs(varS.quant.errorGRFVert)),...
            max(abs(varS.quant.errorGRFVertFull))];
        varS.quant.errorGRFMaxND = varS.quant.errorGRFMax./...
            (varS.params.m*varS.params.g);


        % Calculate RMSE for GRF tracking error
        varS.quant.rmseGRFHor = sqrt(mean(varS.quant.errorGRFHor.^2,1));
        varS.quant.rmseGRFVert = sqrt(mean(varS.quant.errorGRFVert.^2,1));
        
        varS.quant.rmseGRFHorND = varS.quant.rmseGRFHor./...
            (varS.params.m*varS.params.g);
        varS.quant.rmseGRFVertND = varS.quant.rmseGRFVert./...
            (varS.params.m*varS.params.g);

        varS.quant.rmseGRFHorFull =...
            sqrt(mean(varS.quant.errorGRFHorFull.^2,1));
        varS.quant.rmseGRFVertFull =...
            sqrt(mean(varS.quant.errorGRFVertFull.^2,1));
        
        varS.quant.rmseGRFHorNDFull = varS.quant.rmseGRFHorFull./...
            (varS.params.m*varS.params.g);
        varS.quant.rmseGRFVertNDFull = varS.quant.rmseGRFVertFull./...
            (varS.params.m*varS.params.g);
        
    end
    
    % Call to save function to store analysis
    saveData(varS, datS)

end