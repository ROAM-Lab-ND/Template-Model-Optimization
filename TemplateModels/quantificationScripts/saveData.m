%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% quantOptim %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 8 February 2022
% Last Updated: 8 February 2022

% This function is used to save off data to 

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function saveData(vS,dS)

    quant2CSV = {'RMSE COM Segment ND', 'RMSE COM Temporal ND',...
        'RMSE GRF Horizontal ND', 'RMSE GRF Vertical ND', 'RMSE Phase'};
    
    if isfield(dS.data,'grf')
    
        quant2CSVVals = {vS.quant.rmseCOMSegND, vS.quant.rmseCOMTimeND,...
        vS.quant.rmseGRFHorNDFull, vS.quant.rmseGRFVertNDFull,...
        vS.quant.rmsePhase};
    
    else
        
        quant2CSVVals = {vS.quant.rmseCOMSegND, vS.quant.rmseCOMTimeND,...
        '--', '--', vS.quant.rmsePhase};
    
    end
    
    if strcmp(vS.params.model, 'VPP')
        
        optims2CSV = {'Optimal Spring', 'Optimal Phases',...
            'Optimal Leg Length', 'Optimal VPP'};
        
        if strcmp(vS.params.springType, 'Constant')
            
            optims2CSVVals = {vS.optims.kOptLag(1), vS.optims.timeOpt,...
                vS.optims.len0Opt, vS.optims.rVPP};
            
        else
            
            optims2CSVVals = {'--', vS.optims.timeOpt,...
                vS.optims.len0Opt, vS.optims.rVPP};
            
        end
        
        titleCell = {[vS.params.model '_' vS.params.method '_'...
            vS.params.springType 'Spring_' vS.params.vppType 'VPP'],...
            '-----------------------'};
        
    else
        
        optims2CSV = {'Optimal Spring', 'Optimal Phases',...
            'Optimal Leg Length',};
        
        if strcmp(vS.params.springType, 'Constant')
            
            optims2CSVVals = {vS.optims.kOptLag(1), vS.optims.timeOpt,...
                vS.optims.len0Opt};
            
        else
            
            optims2CSVVals = {'--', vS.optims.timeOpt,...
                vS.optims.len0Opt};
            
        end
        
        titleCell = {[vS.params.model '_' vS.params.method '_'...
            vS.params.springType 'Spring'], '-----------------------'};
        
    end
    
    subj2CSV = {'Data Source', 'Subject ID', 'Subject Trial', 'Indices',...
        'Subject Weight', 'Subject Leg Length', 'Subject Treadmill',...
        'Subject Phase Durations'};
    subj2CSVVals = {vS.params.source, dS.subj.subjectID,...
        dS.subj.testNum, [dS.inds.start, dS.inds.end], dS.subj.weight,...
        dS.subj.len0, dS.subj.treadmill, dS.data.timeFull};
    
    quantCell = [quant2CSV', quant2CSVVals'];
    optimsCell = [optims2CSV', optims2CSVVals'];
    subjCell = [subj2CSV', subj2CSVVals'];
    
    spacerCell = {'-----------------------', '-----------------------'};
    
    cell2CSV = [titleCell; subjCell; optimsCell; quantCell; spacerCell];
    
    %fileLoc = cd;
    %fileLoc = fullfile(fileLoc, '\Data\Subject05\2023-12-21');
    
    fileLoc = 'C:\Users\dkell\Documents\MATLAB\Research\OSL_Modeling\OSL_Templates\HumanData\SubjectData\Subject05\2023-12-21';
    
    writecell(cell2CSV, [fileLoc '\' vS.params.model 'Data.csv'],...
        'WriteMode', 'append')
    
    if strcmp(vS.params.model, 'VPP')
        
        if strcmp(vS.params.vppType, 'Constant') &&...
                (abs(vS.params.rVPPD - vS.params.rVPPS) == 0)
    
            fileName = [fileLoc '\s' dS.subj.subjectID '_t'...
                dS.subj.testNum '_i' num2str(dS.inds.start)...
                num2str(dS.inds.end) '_' vS.params.model '_'...
                vS.params.method '_' vS.params.springType 'K_'...
                vS.params.vppType 'SameVPP.mat'];
            
        else
            
            fileName = [fileLoc '\s' dS.subj.subjectID '_t'...
                dS.subj.testNum '_i' num2str(dS.inds.start)...
                num2str(dS.inds.end) '_' vS.params.model '_'...
                vS.params.method '_' vS.params.springType 'K_'...
                vS.params.vppType 'DiffVPP.mat'];
            
        end
        
    else
        
        fileName = [fileLoc '\s' dS.subj.subjectID '_t' dS.subj.testNum...
            '_i' num2str(dS.inds.start) num2str(dS.inds.end) '_'...
            vS.params.model '_' vS.params.method '_'...
            vS.params.springType 'K.mat'];
        
    end
    
    save(fileName, 'vS', 'dS')
    
    %% YAML Structure Creation
    
    fileName = ['C:\Users\dkell\Desktop\paramsTLC00_PWS.yaml'];
    
    tS = struct;
    tS.posXCoeff = dS.fit.coefPosHorF;
    tS.posZCoeff = dS.fit.coefPosVertF;
    tS.posX = vS.optims.stateFull(:,1)';
    tS.posZ = vS.optims.stateFull(:,2)';
    tS.forceX = vS.optims.grfHorFull(:,2)';
    tS.forceZ = vS.optims.grfVertFull(:,2)';
    tS.gaitCycle = vS.optims.timeNormFull*100;
    tS.gaitDur = tS.gaitCycle(300);
    tS.stiff = vS.optims.stateFull(:,5)';
    tS.len0 = vS.optims.len0Opt;
    
    WriteYaml(fileName, tS)
    
end