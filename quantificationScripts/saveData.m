%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% quantOptim %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 8 February 2022
% Last Updated: 8 February 2022

% This function is used to save off data to various files for future use

% INPUTS:
%   varS - Structure used for variable and post-optimization storage
%   datS - Structure used for subject data storage

% OUTPUTS:
%   NONE

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function saveData(varS, datS)

    % Create cell for storing metrics in .csv file
    quant2CSV = {'RMSE COM Segment ND', 'RMSE COM Temporal ND',...
        'RMSE GRF Horizontal ND', 'RMSE GRF Vertical ND', 'RMSE Phase'};
    
    % Check if GRF data exists for subject
    if isfield(datS.data,'grf')
    
        % Create cell of metric values to be stored in .csv file
        quant2CSVVals = {varS.quant.rmseCOMSegND,...
            varS.quant.rmseCOMTimeND, varS.quant.rmseGRFHorNDFull,...
            varS.quant.rmseGRFVertNDFull, varS.quant.rmsePhase};
    
    else
        
        % Create cell of metric values to be stored in .csv file
        quant2CSVVals = {varS.quant.rmseCOMSegND,...
            varS.quant.rmseCOMTimeND, '--', '--', varS.quant.rmsePhase};
    
    end
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')
        
        % Create cell for storing data in .csv file
        optims2CSV = {'Optimal Spring', 'Optimal Phases',...
            'Optimal Leg Length', 'Optimal VPP'};
        
        % Check if using constant or varying spring stiffness
        if strcmp(varS.params.springType, 'Constant')
            
            % Create cell of data values to be stored in .csv file
            optims2CSVVals = {varS.optims.kOptLag(1),...
                varS.optims.timeOpt, varS.optims.len0Opt,...
                varS.optims.rVPP};
            
        else
            
            % Create cell of data values to be stored in .csv file
            optims2CSVVals = {'--', varS.optims.timeOpt,...
                varS.optims.len0Opt, varS.optims.rVPP};
            
        end
        
        % Create header for .csv file
        titleCell = {[varS.params.model '_' varS.params.method '_'...
            varS.params.springType 'Spring_' varS.params.vppType 'VPP'],...
            '-----------------------'};
        
    else
        
        % Create cell for storing data in .csv file 
        optims2CSV = {'Optimal Spring', 'Optimal Phases',...
            'Optimal Leg Length',};
        
        % Check if using constant or varying spring stiffness
        if strcmp(varS.params.springType, 'Constant')
            
            % Create cell of data values to be stored in .csv file
            optims2CSVVals = {varS.optims.kOptLag(1),...
                varS.optims.timeOpt,...
                varS.optims.len0Opt};
            
        else
            
            % Create cell of data values to be stored in .csv file
            optims2CSVVals = {'--', varS.optims.timeOpt,...
                varS.optims.len0Opt};
            
        end
        
        % Create header for .csv file
        titleCell = {[varS.params.model '_' varS.params.method '_'...
            varS.params.springType 'Spring'], '-----------------------'};
        
    end
    
    % Create cell of subject information to be stored in .csv file
    subj2CSV = {'Data Source', 'Subject ID', 'Subject Trial', 'Indices',...
        'Subject Weight', 'Subject Leg Length', 'Subject Treadmill',...
        'Subject Phase Durations'};
    
    % Create cell of subject data to be stored in .csv file
    subj2CSVVals = {varS.params.source, datS.subj.subjectID,...
        datS.subj.testNum, [datS.inds.start, datS.inds.end],...
        datS.subj.weight, datS.subj.len0, datS.subj.treadmill,...
        datS.data.timeFull};
    
    % Combine cell headers and data
    quantCell = [quant2CSV', quant2CSVVals'];
    optimsCell = [optims2CSV', optims2CSVVals'];
    subjCell = [subj2CSV', subj2CSVVals'];
    
    % Create spacer cell for easier reading in .csv file
    spacerCell = {'-----------------------', '-----------------------'};
    
    % Combine all information to be stored in .csv file
    cell2CSV = [titleCell; subjCell; optimsCell; quantCell; spacerCell];
    
    % Grab current directory where main script has been called, go to Data
    fileLoc = cd;
    fileLoc = fullfile(fileLoc, '\Data');
    
    % Write data to .csv file
    writecell(cell2CSV, [fileLoc '\' varS.params.model 'Data.csv'],...
        'WriteMode', 'append')
    
    % Check if using VPP or BSLIP template model
    if strcmp(varS.params.model, 'VPP')
        
        % Check if using constant or varying VP, if VP is same or different
        if strcmp(varS.params.vppType, 'Constant') &&...
                (abs(varS.params.rVPPD - varS.params.rVPPS) == 0)
    
            % Create filename for storing matlab variables
            fileName = [fileLoc '\s' datS.subj.subjectID '_t'...
                datS.subj.testNum '_i' num2str(datS.inds.start)...
                num2str(datS.inds.end) '_' varS.params.model '_'...
                varS.params.method '_' varS.params.springType 'K_'...
                varS.params.vppType 'SameVPP.mat'];
            
        else
            
            % Create filename for storing matlab variables
            fileName = [fileLoc '\s' datS.subj.subjectID '_t'...
                datS.subj.testNum '_i' num2str(datS.inds.start)...
                num2str(datS.inds.end) '_' varS.params.model '_'...
                varS.params.method '_' varS.params.springType 'K_'...
                varS.params.vppType 'DiffVPP.mat'];
            
        end
        
    else
        
        % Create filename for storing matlab variables
        fileName = [fileLoc '\s' datS.subj.subjectID '_t'...
            datS.subj.testNum '_i' num2str(datS.inds.start)...
            num2str(datS.inds.end) '_' varS.params.model '_'...
            varS.params.method '_' varS.params.springType 'K.mat'];
        
    end
    
    % Save variable and subject data structures to .mat file
    save(fileName, 'varS', 'datS')
    
end