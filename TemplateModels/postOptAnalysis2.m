%% %%%%%%%%%%%%%%%%%%%%%%%%%% postOptAnalysis %%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 27 June 2022
% Last Updated: 27 June 2022

% This code uses .mat files from optimized template models (see mainCRC.m)
% to analyze the resulting optimized models 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

close all
clearvars -except qS dataStruct
clc

tempPath = fileparts(which('postOptAnalysis.m'));
addpath(genpath(tempPath));

tempPath = cd;
addpath([tempPath '\Data']);
addpath('C:\Users\dkell\Documents\MATLAB\CasADi');

import casadi.*

%% %%%%%%%%%%%%%%%%%%%%%%%%% Analysis Section %%%%%%%%%%%%%%%%%%%%%%%%%% %%

boolGraph = true;

folderName = uigetdir();
fileNames = dir(fullfile(folderName, '*.mat'));
matFiles = {fileNames.name};

fileNames = dir(fullfile(folderName, '*.txt'));
txtFiles = {fileNames.name};

if ~exist('qS', 'var')
    
    qS = struct;
    
end

count = 0;
subjIDStart = 1;
subjIDEnd = 24;
trialIDStart = 1;
trialIDEnd = 8;
numSubjs = (subjIDEnd - subjIDStart + 1);
numTrials = (trialIDEnd - trialIDStart + 1);

errCOM = NaN(numSubjs, numTrials);
errCOMH = errCOM;
errCOMV = errCOM;
errGRF = errCOM;
errPhase = errCOM;
rangeTorso = errCOM;
optFlag = errCOM;
timeSolve = errCOM;
cellK1 = cell(numSubjs, numTrials);
cellK2 = cell(numSubjs, numTrials);
cellT = cell(numSubjs, numTrials);
cellTau1 = cell(numSubjs, numTrials);
cellTau2 = cell(numSubjs, numTrials);
iter = 0;

for i = ((subjIDStart - 1)*8 + trialIDStart):((subjIDEnd - 1)*8 + trialIDEnd)
    
    load(matFiles{i});
    
    iter = iter + 1;
    indSubj = floor((iter-1)/numTrials) + 1;
    
%     disp(['Treadmill Speed: ' num2str(dS.subj.treadmill)])
    
    if ~mod(i,8)
        
        indTrial = 8;
        
    else
        
        indTrial = mod(i,8);
        
    end
    
    errCOM(indSubj, indTrial) = vS.quant.rmseCOMTimeND;
    errCOMH(indSubj, indTrial) = vS.quant.rmseCOMTimeHorND;
    errCOMV(indSubj, indTrial) = vS.quant.rmseCOMTimeVertND;
    errGRF(indSubj, indTrial) = mean(vS.quant.rmseGRFVertNDFull);
    errPhase(indSubj, indTrial) = vS.quant.rmsePhase/dS.data.timeGait;
    
    rangeTorso(indSubj, indTrial) = range(vS.optims.stateFull(:,3)*180/pi);
    
    springL1 = [];
    springL2 = springL1;

    for j = 1:(vS.params.steps/2)

       springL1 = [springL1; vS.optims.stateFull((4*(j-1)*100+1):(4*(j-1)*100+100), end-1);...
           zeros(100,1); vS.optims.stateFull((4*(j-1)*100+201):(4*(j-1)*100+400),end)];
       springL2 = [springL2; vS.optims.stateFull((4*(j-1)*100+1):(4*(j-1)*100+200), end);...
           vS.optims.stateFull((4*(j-1)*100+201):(4*(j-1)*100+300),end-1); zeros(100,1)];

    end
    
    cellK1(indSubj, indTrial) = {springL1*...
        (dS.subj.len0/dS.subj.bodyWeight)*1000};
    cellK2(indSubj, indTrial) = {springL2*...
        (dS.subj.len0/dS.subj.bodyWeight)*1000};
    
    if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')
    
        cellTau1(indSubj, indTrial) = {vS.optims.tauLagFull/vS.params.m};
        cellTau2(indSubj, indTrial) = {vS.optims.tauLeadFull/vS.params.m};
        
    end
    
    cellT(indSubj, indTrial) = {vS.optims.timeNormFull};
    
    if iter == 1
        
        if contains(matFiles{i}, 'VPP')
            
            modelType = 'VPP';
            
        else
            
            modelType = 'BSLIP';
            
        end
        
        if contains(matFiles{i}, 'VaryingK')
            
            springType = 'varyK';
            
        else
            
            springType = 'constK';
            
        end
        
    end
    
    if litcount(txtFiles{i}, 'Optimal Solution Found.')
        
        %disp(['Subject ' num2str(floor((i-1)/8)+1) ' Trial '...
            %num2str(mod(i,8)) ' Optimally Solved.'])
        count = count + 1;
        
        timeSolve(indSubj, indTrial) = littime(txtFiles{i}, 'Total CPU secs');
        
        if any(dS.data.timeFull < 0)
        
            optFlag(indSubj, indTrial) = 1;
            disp(['Subject ' num2str(indSubj)...
                ' Trial ' num2str(indTrial)...
                ' removed due to incorrect phase calculation.']);
            
        end
        
    else
        
        %disp(['Subject ' num2str(floor((i-1)/8)+1) ' Trial '...
            %num2str(mod(i,8)) ' NOT OPTIMAL.'])
        optFlag(indSubj, indTrial) = 1;
    
    end
    
end

disp(count);

errCOM(optFlag == 1) = NaN;
errCOMAvg = mean(errCOM, 1, 'omitnan');
errCOMSTD = std(errCOM, 0, 1, 'omitnan');
rmTrialInd = ones(size(errCOM));
rmTrialInd(optFlag == 1) = 0;

errCOMMap = ones(size(errCOM));
errCOMMap(isnan(errCOM)) = 0;

%%
for j = 1:length(errCOM(1,:))
    
    for i = 1:length(errCOM(:,1))
        
        if (errCOM(i, j) >= (errCOMAvg(j) + 2*errCOMSTD(j))) %||...
                %(errCOM(i, j) <= (errCOMAvg(j) - 2*errCOMSTD(j)))
            
            errCOM(i, j) = NaN;
            rmTrialInd(i, j) = 0;
            
        end
            
    end
    
end

for j = 1:length(errCOM(1, :))
    
    for i = 1:length(errCOM(:, 1))
        
        tempInd = rmTrialInd(:, j);
        tempInd(i) = 0;
            
        avgCOM = mean(errCOM(tempInd == 1, j), 'omitnan');
        stdCOM = std(errCOM(tempInd == 1, j), 0, 1, 'omitnan');
        
        if (errCOM(i, j) >= (avgCOM + 2*stdCOM)) %||...
               %(errCOM(i, j) <= (avgCOM - 2*stdCOM))
            
            errCOM(i, j) = NaN;
            rmTrialInd(i, j) = 0;
            
        end
            
    end
    
end

errCOMMap((isnan(errCOM) & (errCOMMap ~= 0))) = 0.5;

errCOMH(rmTrialInd == 0) = NaN;
errCOMV(rmTrialInd == 0) = NaN;
errGRF(rmTrialInd == 0) = NaN;
errPhase(rmTrialInd == 0) = NaN;
rangeTorso(rmTrialInd == 0) = NaN;
timeSolve(rmTrialInd == 0) = NaN;
cellK1(rmTrialInd == 0) = {NaN(length(springL1), 1)};
cellK2(rmTrialInd == 0) = {NaN(length(springL1), 1)};
cellT(rmTrialInd == 0) = {NaN(1, 400)};

if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')

    cellTau1(rmTrialInd == 0) = {NaN(1, length(springL1))};
    cellTau2(rmTrialInd == 0) = {NaN(1, length(springL1))};
    
end

%%
stiffAvg1 = cell(1, length(cellK1(1,:)));
stiffAvg2 = stiffAvg1;

stiffSTD1 = stiffAvg1;
stiffSTD2 = stiffAvg1;

if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')

    tauAvg1 = stiffAvg1;
    tauAvg2 = stiffAvg1;

    tauSTD1 = stiffAvg1;
    tauSTD2 = stiffAvg1;
    
end

timeAvg = stiffAvg1;
timeSTD = stiffAvg1;

for j = 1:length(cellK1(1,:))
    
    allK1 = reshape(cell2mat(cellK1(:,j)), 400, []);
    allK2 = reshape(cell2mat(cellK2(:,j)), 400, []);
    
    if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')
    
        allTau1 = reshape(cell2mat(cellTau1(:,j)), [], 400);
        allTau2 = reshape(cell2mat(cellTau2(:,j)), [], 400);
    
        tauSTD1(j) = {std(allTau1, 0, 1, 'omitnan')};
        tauSTD2(j) = {std(allTau2, 0, 1, 'omitnan')};
        
        tauAvg1(j) = {mean(allTau1, 1, 'omitnan')};
        tauAvg2(j) = {mean(allTau2, 1, 'omitnan')};
        
    end
    
    stiffAvg1(j) = {mean(allK1, 2, 'omitnan')};
    stiffAvg2(j) = {mean(allK2, 2, 'omitnan')};
    
    stiffSTD1(j) = {std(allK1, 0, 2, 'omitnan')};
    stiffSTD2(j) = {std(allK2, 0, 2, 'omitnan')};
    
    timeAvg = {mean(cell2mat(cellT(:,j)), 1, 'omitnan')};
    timeSTD = {std(cell2mat(cellT(:,j)), 0, 1, 'omitnan')};
    
end

for i = 1:numel(stiffAvg1)
    
    stiffAvg1{i}(stiffAvg1{i}==0) = NaN;
    stiffAvg2{i}(stiffAvg2{i}==0) = NaN;
    
%     if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')
%     
%         tauAvg1{i}(tauAvg1{i}==0) = NaN;
%         tauAvg2{i}(tauAvg2{i}==0) = NaN;
%         
%     end
    
end

errCOMAvg = mean(errCOM, 1, 'omitnan');
errCOMSTD = std(errCOM, 0, 1, 'omitnan');

errCOMHAvg = mean(errCOMH, 1, 'omitnan');
errCOMHSTD = std(errCOMH, 0, 1, 'omitnan');

errCOMVAvg = mean(errCOMV, 1, 'omitnan');
errCOMVSTD = std(errCOMV, 0, 1, 'omitnan');

errGRFAvg = mean(errGRF, 1, 'omitnan');
errGRFSTD = std(errGRF, 0, 1, 'omitnan');

errPhaseAvg = mean(errPhase, 1, 'omitnan');
errPhaseSTD = std(errPhase, 0, 1, 'omitnan');

rangeTorsoAvg = mean(rangeTorso, 1, 'omitnan');
rangeTorsoSTD = std(rangeTorso, 0, 1, 'omitnan');

timeSolveAvg = mean(timeSolve, 1, 'omitnan');
timeSolveSTD = std(timeSolve, 0, 1, 'omitnan');

qS.(modelType).(springType).errMap = errCOMMap;

qS.(modelType).(springType).errCOM = errCOM;
qS.(modelType).(springType).errCOMAvg = errCOMAvg; %mean(errCOM, 1);
qS.(modelType).(springType).errCOMSTD = errCOMSTD; %std(errCOM, 1);

qS.(modelType).(springType).errCOMHor = errCOMH;
qS.(modelType).(springType).errCOMHorAvg = errCOMHAvg; %mean(errCOMH, 1);
qS.(modelType).(springType).errCOMHorSTD = errCOMHSTD; %std(errCOMH, 1);

qS.(modelType).(springType).errCOMVert = errCOMV;
qS.(modelType).(springType).errCOMVertAvg = errCOMVAvg; %mean(errCOMV, 1);
qS.(modelType).(springType).errCOMVertSTD = errCOMVSTD; %std(errCOMV, 1);

qS.(modelType).(springType).errGRF = errGRF;
qS.(modelType).(springType).errGRFAvg = errGRFAvg; %mean(errGRF, 1);
qS.(modelType).(springType).errGRFSTD = errGRFSTD; %std(errGRF, 1);

qS.(modelType).(springType).errPhase = errPhase;
qS.(modelType).(springType).errPhaseAvg = errPhaseAvg; %mean(errPhase, 1);
qS.(modelType).(springType).errPhaseSTD = errPhaseSTD; %std(errPhase, 1);

qS.(modelType).(springType).K1 = cellK1;
qS.(modelType).(springType).K2 = cellK2;
qS.(modelType).(springType).K1Avg = stiffAvg1;
qS.(modelType).(springType).K2Avg = stiffAvg2;
qS.(modelType).(springType).K1STD = stiffSTD1;
qS.(modelType).(springType).K2STD = stiffSTD2;

qS.(modelType).(springType).timeNormAvg = timeAvg;
qS.(modelType).(springType).timeNormSTD = timeSTD;

qS.(modelType).(springType).torsoAvg = rangeTorsoAvg;
qS.(modelType).(springType).torsoSTD = rangeTorsoSTD;

qS.(modelType).(springType).timeSolve = timeSolve;
qS.(modelType).(springType).timeSolveAvg = timeSolveAvg; %mean(errCOM, 1);
qS.(modelType).(springType).timeSolveSTD = timeSolveSTD; %std(errCOM, 1);

if strcmp(vS.params.model, 'VPP') && isfield(vS.optims,'tauLagFull')
    
    qS.(modelType).(springType).tau1 = cellTau1;
    qS.(modelType).(springType).tau2 = cellTau2;
    qS.(modelType).(springType).tau1Avg = tauAvg1;
    qS.(modelType).(springType).tau2Avg = tauAvg2;
    qS.(modelType).(springType).tau1STD = tauSTD1;
    qS.(modelType).(springType).tau2STD = tauSTD2;
    
end

qS.(modelType).(springType).optFlag = optFlag;

if boolGraph
    
   dataVV.dS.fit.grfHorSD = dataStruct.fit.grfHorSD;
   dataVV.dS.fit.grfVertSD = dataStruct.fit.grfVertSD;
    
   errCOMAll = [qS.BSLIP.constK.errCOMAvg; qS.BSLIP.varyK.errCOMAvg; qS.VPP.constK.errCOMAvg; qS.VPP.varyK.errCOMAvg];
   errGRFAll = [qS.BSLIP.constK.errGRFAvg; qS.BSLIP.varyK.errGRFAvg; qS.VPP.constK.errGRFAvg; qS.VPP.varyK.errGRFAvg];
   errPhaseAll = [qS.BSLIP.constK.errPhaseAvg; qS.BSLIP.varyK.errPhaseAvg; qS.VPP.constK.errPhaseAvg; qS.VPP.varyK.errPhaseAvg];
       
   [qS.stats.avg.pC{1}, qS.stats.avg.hC{1}] = ttest2(errCOMAll(1,:), errCOMAll(2,:));
   [qS.stats.avg.pC{2}, qS.stats.avg.hC{2}] = ttest2(errCOMAll(3,:), errCOMAll(4,:));
   [qS.stats.avg.pC{3}, qS.stats.avg.hC{3}] = ttest2(errCOMAll(1,:), errCOMAll(3,:));
   [qS.stats.avg.pC{4}, qS.stats.avg.hC{4}] = ttest2(errCOMAll(2,:), errCOMAll(4,:));
   
   [qS.stats.avg.pG{1}, qS.stats.avg.hG{1}] = ttest2(errGRFAll(1,:), errGRFAll(2,:));
   [qS.stats.avg.pG{2}, qS.stats.avg.hG{2}] = ttest2(errGRFAll(3,:), errGRFAll(4,:));
   [qS.stats.avg.pG{3}, qS.stats.avg.hG{3}] = ttest2(errGRFAll(1,:), errGRFAll(3,:));
   [qS.stats.avg.pG{4}, qS.stats.avg.hG{4}] = ttest2(errGRFAll(2,:), errGRFAll(4,:));
   
   [qS.stats.avg.pP{1}, qS.stats.avg.hP{1}] = ttest2(errPhaseAll(1,:), errPhaseAll(2,:));
   [qS.stats.avg.pP{2}, qS.stats.avg.hP{2}] = ttest2(errPhaseAll(3,:), errPhaseAll(4,:));
   [qS.stats.avg.pP{3}, qS.stats.avg.hP{3}] = ttest2(errPhaseAll(1,:), errPhaseAll(3,:));
   [qS.stats.avg.pP{4}, qS.stats.avg.hP{4}] = ttest2(errPhaseAll(2,:), errPhaseAll(4,:));
   
   for i = 1:8
       
       errCOMTemp = [qS.BSLIP.constK.errCOM(:,i) qS.BSLIP.varyK.errCOM(:,i)...
           qS.VPP.constK.errCOM(:,i) qS.VPP.varyK.errCOM(:,i)];
       errGRFTemp = [qS.BSLIP.constK.errGRF(:,i) qS.BSLIP.varyK.errGRF(:,i)...
           qS.VPP.constK.errGRF(:,i) qS.VPP.varyK.errGRF(:,i)];
       errPhaseTemp = [qS.BSLIP.constK.errPhase(:,i) qS.BSLIP.varyK.errPhase(:,i)...
           qS.VPP.constK.errPhase(:,i) qS.VPP.varyK.errPhase(:,i)];
       
       [qS.stats.trial.pC{i,1}, qS.stats.trial.hC{i,1}] = signrank(errCOMTemp(:,1), errCOMTemp(:,2));
       [qS.stats.trial.pC{i,2}, qS.stats.trial.hC{i,2}] = signrank(errCOMTemp(:,3), errCOMTemp(:,4));
       [qS.stats.trial.pC{i,3}, qS.stats.trial.hC{i,3}] = signrank(errCOMTemp(:,1), errCOMTemp(:,3));
       [qS.stats.trial.pC{i,4}, qS.stats.trial.hC{i,4}] = signrank(errCOMTemp(:,2), errCOMTemp(:,4));

       [qS.stats.trial.pG{i,1}, qS.stats.trial.hG{i,1}] = signrank(errGRFTemp(:,1), errGRFTemp(:,2));
       [qS.stats.trial.pG{i,2}, qS.stats.trial.hG{i,2}] = signrank(errGRFTemp(:,3), errGRFTemp(:,4));
       [qS.stats.trial.pG{i,3}, qS.stats.trial.hG{i,3}] = signrank(errGRFTemp(:,1), errGRFTemp(:,3));
       [qS.stats.trial.pG{i,4}, qS.stats.trial.hG{i,4}] = signrank(errGRFTemp(:,2), errGRFTemp(:,4));

       [qS.stats.trial.pP{i,1}, qS.stats.trial.hP{i,1}] = signrank(errPhaseTemp(:,1), errPhaseTemp(:,2));
       [qS.stats.trial.pP{i,2}, qS.stats.trial.hP{i,2}] = signrank(errPhaseTemp(:,3), errPhaseTemp(:,4));
       [qS.stats.trial.pP{i,3}, qS.stats.trial.hP{i,3}] = signrank(errPhaseTemp(:,1), errPhaseTemp(:,3));
       [qS.stats.trial.pP{i,4}, qS.stats.trial.hP{i,4}] = signrank(errPhaseTemp(:,2), errPhaseTemp(:,4));
       
       if ismember(i, [2 5 7])
           
           figTitle = {'B-SLIP (C)', 'B-SLIP (V)', 'VPP (C)', 'VPP (V)'};
       
           for j = 1:4
           
               linCG = table(errCOMTemp(:,j), errGRFTemp(:,j),...
                   'VariableNames', {'CoM Tracking Error', 'GRF Tracking Error'});
               mdlCG = fitlm(linCG);

               figure
               plotAdded(mdlCG);
               ax = gca;
               title(figTitle{j}, 'Interpreter', 'latex')
               xlabel('CoM Tracking Error', 'Interpreter', 'latex')
               ylabel('GRF Tracking Error', 'Interpreter', 'latex')
               ax.Legend.String(3) = {'95\% conf. bounds'};
               set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
               set(gcf,'position',[10, 10, 560, 300])
           
           end
           
       end
       %}
       
   end
    
   xTicks = {'40\%', '55\%', '70\%', '85\%', '100\%', '115\%', '130\%', '145\%'};
   %xTicks = {'1', '2', '3', '4', '5', '6', '7', '8'};
   barX = categorical(xTicks);
   barX = reordercats(barX, xTicks);
   
   fieldQS1 = fieldnames(qS);
   numFields1 = numel(fieldQS1)-1;
   fieldQS2 = {};
   barYCOMFull = [];
   barYGRFFull = [];
   barYPhaseFull = [];
   
   barYCOM2Full = [];
   barYGRF2Full = [];
   barYPhase2Full = [];
   dataStrFull = {};
   
   for i = 1:numFields1
       
       numFields2 = numel(fieldnames(qS.(fieldQS1{i})));
       fieldQS2{i} = fieldnames(qS.(fieldQS1{i}));
       
       for j = 1:numFields2
           
           barYCOM = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errCOMAvg;
           barYCOMFull = [barYCOMFull; barYCOM];
           barYCOM2 = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errCOMSTD;
           negBar = zeros(size(barYCOM2));
           barYCOM2Full = [barYCOM2Full; barYCOM2];
           
           barYGRF = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errGRFAvg;
           barYGRFFull = [barYGRFFull; barYGRF];
           barYGRF2 = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errGRFSTD;
           barYGRF2Full = [barYGRF2Full; barYGRF2];
           
           barYPhase = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errPhaseAvg;
           barYPhaseFull = [barYPhaseFull; barYPhase];
           barYPhase2 = qS.(fieldQS1{i}).(fieldQS2{i}{j}).errPhaseSTD;
           barYPhase2Full = [barYPhase2Full; barYPhase2];
           
           dataStr = [fieldQS1{i} ' ' fieldQS2{i}{j}];
           dataStrFull{end+1} = dataStr;
           
           %{
           figure
           hold on
           bar(barX, barYCOM, 0.5, 'FaceColor', [0.6 0.9 0.1])
           errorbar(barX, barYCOM, negBar, barYCOM2,...
               'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
               'LineWidth', 1, 'CapSize', 10)
           xlabel('errPhase Speed (\% of Preferred Walking Speed)')
           ylabel(['COM Error Magnitude '...
               '[$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]'])
           legend(dataStr)
           set(gca, 'FontSize', 14)
           
           figure
           hold on
           bar(barX, barYGRF, 0.5, 'FaceColor', [0.1 0.7 0.9])
           errorbar(barX, barYGRF, negBar, barYGRF2,...
               'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
               'LineWidth', 1, 'CapSize', 10)
           xlabel('Trial Speed (\% of Preferred Walking Speed)')
           ylabel(['GRF Error Magnitude '...
               '[$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]'])
           legend(dataStr)
           set(gca, 'FontSize', 14)
           
           figure
           hold on
           bar(barX, barYPhase, 0.5, 'FaceColor', [1.0 0.7 0.0])
           errorbar(barX, barYPhase, negBar, barYPhase2,...
               'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
               'LineWidth', 1, 'CapSize', 10)
           xlabel('Trial Speed (\% of Preferred Walking Speed)')
           ylabel('Phase Error Magnitude [sec]')
           legend(dataStr)
           set(gca, 'FontSize', 14)
           
           for k = 1:8
               
               timeVecP = [cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).timeNormAvg)...
                   fliplr(cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).timeNormAvg))];
               timeVec = cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).timeNormAvg);

               dataVec1 = cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K1Avg(k));
               dataVec2 = cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K2Avg(k));

               patchVec1 = [(cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K1Avg(k))' +...
                   cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K1STD(k))')...
                   fliplr(cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K1Avg(k))' -...
                   cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K1STD(k))')];
               patchVec2 = [(cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K2Avg(k))' +...
                   cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K2STD(k))')...
                   fliplr(cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K2Avg(k))' -...
                   cell2mat(qS.(fieldQS1{i}).(fieldQS2{i}{j}).K2STD(k))')];
               
               isNANCheck1a = find(isnan(patchVec1), 1, 'first');
               isNANCheck1b = find(isnan(patchVec1), 1, 'last');
               isNANCheck2a = find(isnan(patchVec2), 1, 'first');
               isNANCheck2b = find(isnan(patchVec2), 1, 'last');
               
               timeVecP1a = timeVecP(1:find(isnan(patchVec1), 1, 'first')-1);
               timeVecP1b = timeVecP(find(isnan(patchVec1), 1, 'last')+1:end);
               
               timeVecP2a = timeVecP(1:find(isnan(patchVec2), 1, 'first')-1);
               timeVecP2b = timeVecP(find(isnan(patchVec2), 1, 'last')+1:end);
               
               patchVecP1a = patchVec1(1:find(isnan(patchVec1), 1, 'first')-1);
               patchVecP1b = patchVec1(find(isnan(patchVec1), 1, 'last')+1:end);
               
               patchVecP2a = patchVec2(1:find(isnan(patchVec2), 1, 'first')-1);
               patchVecP2b = patchVec2(find(isnan(patchVec2), 1, 'last')+1:end);
           
           %}
               
%                figure
%                hold on
%                patch(timeVecP1a, patchVecP1a, [0.0 0.8 0.9], 'FaceAlpha', 0.3)
%                patch(timeVecP1b, patchVecP1b, [0.0 0.8 0.9], 'FaceAlpha', 0.3)
%                patch(timeVecP2a, patchVecP2a, [1.0 0.8 0.9], 'FaceAlpha', 0.3)
%                patch(timeVecP2b, patchVecP2b, [1.0 0.8 0.9], 'FaceAlpha', 0.3)
%                plot(timeVec, dataVec1, 'linewidth', 2, 'color', [0.0 0.1 0.6])
%                plot(timeVec, dataVec2, 'linewidth', 2, 'color', [1.0 0.1 0.1])
%                xline(timeVec(101), 'k');
%                xline(timeVec(201), 'k');
%                xline(timeVec(301), 'k');
%                xlabel('Time []')
%                ylabel(['Spring Stiffness [Trial ' num2str(k) ']'])
%                legend('Mean Stiffness Leg 1', 'Mean Stiffness Leg 2');
%                legend([1 3 5 6], 'Standard Deviation Leg 1',...
%                    'Standard Deviation Leg 2',...
%                    'Mean Stiffness Leg 1',...
%                    'Mean Stiffness Leg 2', 'NumColumns', 2)
               
           %end
           
           %     timeVec = [dS.fit.time fliplr(dS.fit.time)];
           %     dataVec = [(dS.fit.comPosVert + dS.fit.comPosVertSD)...
           %         fliplr(dS.fit.comPosVert - dS.fit.comPosVertSD)];
           %     figure
           %     hold on
           %     patch(timeVec, dataVec, [0.1 0.9 1.0])
           %     plot(dS.fit.time, dS.fit.comPosVert,...
           %         'linewidth', 2, 'color', [0.1 0.0 0.4]);
           %     xlabel('Time [sec]')
           %     ylabel('Vertical COM Displacement [m]')
           %     legend('Standard Deviation', 'Mean Trajectory')
           %     hold off
       
       end
       
   end
   
   barYCOMerr = barYCOMFull(:);
   barYCOM2err = barYCOM2Full(:);
   barYGRFerr = barYGRFFull(:);
   barYGRF2err = barYGRF2Full(:);
   barYPhaseerr = barYPhaseFull(:);
   barYPhase2err = barYPhase2Full(:);
   
   barYCV = barYCOMFull([3 1],:);
   barYCC = barYCOMFull([4 2],:);
   barYGV = barYGRFFull([3 1],:);
   barYGC = barYGRFFull([4 2],:);
   barYPV = barYPhaseFull([3 1],:);
   barYPC = barYPhaseFull([4 2],:);
   
   barYCVS = barYCOM2Full([3 1],:);
   barYCCS = barYCOM2Full([4 2],:);
   barYGVS = barYGRF2Full([3 1],:);
   barYGCS = barYGRF2Full([4 2],:);
   barYPVS = barYPhase2Full([3 1],:);
   barYPCS = barYPhase2Full([4 2],:);
   
   barYCVP = barYCOMFull([2 1],:);
   barYCBS = barYCOMFull([4 3],:);
   barYGVP = barYGRFFull([2 1],:);
   barYGBS = barYGRFFull([4 3],:);
   barYPVP = barYPhaseFull([2 1],:);
   barYPBS = barYPhaseFull([4 3],:);
   
   barYCVPS = barYCOM2Full([2 1],:);
   barYCBSS = barYCOM2Full([4 3],:);
   barYGVPS = barYGRF2Full([2 1],:);
   barYGBSS = barYGRF2Full([4 3],:);
   barYPVPS = barYPhase2Full([2 1],:);
   barYPBSS = barYPhase2Full([4 3],:);
   
   figure
   hold on
   b1 = bar(barYCV');
   ngroups = size(barYCV, 2);
   nbars = size(barYCV, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYCV(i,:), NaN(1,8), barYCVS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [1 0 0];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (V)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYCC');
   ngroups = size(barYCC, 2);
   nbars = size(barYCC, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYCC(i,:), NaN(1,8), barYCCS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   b1(2).FaceColor = [0.4660 0.6740 0.1880];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   %b1(1).FaceColor = [0.9290 0.6940 0.1250];
   %b1(2).FaceColor = [1 0 0];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'VPP (C)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYCVP');
   ngroups = size(barYCVP, 2);
   nbars = size(barYCVP, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYCVP(i,:), NaN(1,8), barYCVPS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [0.4660 0.6740 0.1880];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('VPP (C)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYCBS');
   ngroups = size(barYCBS, 2);
   nbars = size(barYCBS, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYCBS(i,:), NaN(1,8), barYCBSS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [1 0 0];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'B-SLIP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYGV');
   ngroups = size(barYGV, 2);
   nbars = size(barYGV, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYGV(i,:), NaN(1,8), barYGVS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [1 0 0];
   ylim([0 0.4])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (V)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYGC');
   ngroups = size(barYGC, 2);
   nbars = size(barYGC, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYGC(i,:), NaN(1,8), barYGCS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   b1(2).FaceColor = [0.4660 0.6740 0.1880];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   ylim([0 0.4])
   %b1(1).FaceColor = [0.9290 0.6940 0.1250];
   %b1(2).FaceColor = [1 0 0];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'VPP (C)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYGVP');
   ngroups = size(barYGVP, 2);
   nbars = size(barYGVP, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYGVP(i,:), NaN(1,8), barYGVPS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [0.4660 0.6740 0.1880];
   ylim([0 0.4])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('VPP (C)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYGBS');
   ngroups = size(barYGBS, 2);
   nbars = size(barYGBS, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYGBS(i,:), NaN(1,8), barYGBSS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [1 0 0];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   ylim([0 0.4])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   %ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'B-SLIP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYPV');
   ngroups = size(barYPV, 2);
   nbars = size(barYPV, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYPV(i,:), NaN(1,8), barYPVS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [1 0 0];
   ylim([0 0.06])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (V)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYPC');
   ngroups = size(barYPC, 2);
   nbars = size(barYPC, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYPC(i,:), NaN(1,8), barYPCS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   b1(2).FaceColor = [0.4660 0.6740 0.1880];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   ylim([0 0.06])
   %b1(1).FaceColor = [0.9290 0.6940 0.1250];
   %b1(2).FaceColor = [1 0 0];
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'VPP (C)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYPVP');
   ngroups = size(barYPVP, 2);
   nbars = size(barYPVP, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYPVP(i,:), NaN(1,8), barYPVPS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [0.9290 0.6940 0.1250];
   b1(1).FaceColor = [0.4660 0.6740 0.1880];
   ylim([0 0.06])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   ylabel('Phase Error Magnitude [sec]')
   legend('VPP (C)', 'VPP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   figure
   hold on
   b1 = bar(barYPBS');
   ngroups = size(barYPBS, 2);
   nbars = size(barYPBS, 1);
   % Calculating the width for each bar group
   groupwidth = min(0.8, nbars/(nbars + 1.5));
   for i = 1:nbars
       x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
       errorbar(x, barYPBS(i,:), NaN(1,8), barYPBSS(i,:), 'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0], 'LineWidth', 1, 'CapSize', 10);
   end
   hold off
   %b1(1).FaceColor = [0.4660 0.6740 0.1880];
   %b1(2).FaceColor = [0.3010 0.7450 0.9330];
   b1(2).FaceColor = [1 0 0];
   b1(1).FaceColor = [0.3010 0.7450 0.9330];
   ylim([0 0.06])
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   %ylabel('COM Error Magnitude [$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]')
   %ylabel('GRF Error Magnitude [$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]')
   ylabel('Phase Error Magnitude [sec]')
   legend('B-SLIP (C)', 'B-SLIP (V)', 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   set(gca, 'XTickLabel', xTicks);
   
   %{
   figure
   hold on
   b1 = bar(barX, barYCOMFull);
%    errorbar(barX, barYCOMFull, barYCOM2Full,...
%        'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%        'LineWidth', 1, 'CapSize', 10)
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel(['COM Error Magnitude '...
       '[$\frac{\textrm{m RMSE}}{\textrm{m Leg Length}}$]'])
   legend(dataStrFull, 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   b1(1).FaceColor = [0.9290 0.6940 0.1250];
   b1(2).FaceColor = [0.4660 0.6740 0.1880];
   b1(3).FaceColor = [1 0 0];
   b1(4).FaceColor = [0.3010 0.7450 0.9330];
   
   figure
   hold on
   b2 = bar(barX, barYGRFFull);
%    errorbar(barX, barYGRFFull, barYGRF2Full,...
%        'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%        'LineWidth', 1, 'CapSize', 10)
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel(['GRF Error Magnitude '...
       '[$\frac{\textrm{N RMSE}}{\textrm{N Bodyweight}}$]'])
   legend(dataStrFull, 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   b2(1).FaceColor = [0.9290 0.6940 0.1250];
   b2(2).FaceColor = [0.4660 0.6740 0.1880];
   b2(3).FaceColor = [1 0 0];
   b2(4).FaceColor = [0.3010 0.7450 0.9330];
   
   figure
   hold on
   b3 = bar(barX, barYPhaseFull);
%    errorbar(barX, barYPhaseFull, barYPhase2Full,...
%        'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%        'LineWidth', 1, 'CapSize', 10)
   xlabel('Trial Speed (\% of Preferred Walking Speed)')
   ylabel('Phase Error Magnitude [sec]')
   legend(dataStrFull, 'NumColumns', 2)
   set(gca, 'FontSize', 14)
   b3(1).FaceColor = [0.9290 0.6940 0.1250];
   b3(2).FaceColor = [0.4660 0.6740 0.1880];
   b3(3).FaceColor = [1 0 0];
   b3(4).FaceColor = [0.3010 0.7450 0.9330];
   
   for i = 1:8
       
       timeVecV1 = cell2mat(qS.VPP.varyK.timeNormAvg);
       timeVecV2 = cell2mat(qS.VPP.constK.timeNormAvg);
       timeVecB1 = cell2mat(qS.BSLIP.varyK.timeNormAvg);
       timeVecB2 = cell2mat(qS.BSLIP.constK.timeNormAvg);

       dataVec1a = cell2mat(qS.VPP.varyK.K1Avg(i));
       dataVec1b = cell2mat(qS.VPP.varyK.K2Avg(i));
       
       dataVec2a = cell2mat(qS.VPP.constK.K1Avg(i));
       dataVec2b = cell2mat(qS.VPP.constK.K2Avg(i));
       
       dataVec3a = cell2mat(qS.BSLIP.varyK.K1Avg(i));
       dataVec3b = cell2mat(qS.BSLIP.varyK.K2Avg(i));
       
       dataVec4a = cell2mat(qS.BSLIP.constK.K1Avg(i));
       dataVec4b = cell2mat(qS.BSLIP.constK.K2Avg(i));
       
%        figure
%        hold on
%        plot(timeVecV1, dataVec1a, 'linewidth', 2, 'color', [0.0 0.7 1.0])
%        plot(timeVecV1, dataVec1b, 'linewidth', 2, 'color', [0.0 0.4 0.8])
%        plot(timeVecV2, dataVec2a, 'linewidth', 2, 'color', [0.0 0.9 0.4])
%        plot(timeVecV2, dataVec2b, 'linewidth', 2, 'color', [0.0 0.6 0.2])
%        plot(timeVecB1, dataVec3a, 'linewidth', 2, 'color', [1.0 0.5 0.7])
%        plot(timeVecB1, dataVec3b, 'linewidth', 2, 'color', [1.0 0.2 0.0])
%        plot(timeVecB2, dataVec4a, 'linewidth', 2, 'color', [1.0 0.8 0.0])
%        plot(timeVecB2, dataVec4b, 'linewidth', 2, 'color', [1.0 0.4 0.0])
%        xline(timeVec(101), 'k');
%        xline(timeVec(201), 'k');
%        xline(timeVec(301), 'k');
%        xlabel('Time []')
%        ylabel(['Spring Stiffness [Trial ' num2str(i) ']'])
%        legend('Leg 1 VPP (V)', 'Leg 2 VPP (V)',...
%            'Leg 1 VPP (C)', 'Leg 2 VPP (C)',...
%            'Leg 1 BSLIP (V)', 'Leg 2 BSLIP (V)',...
%            'Leg 1 BSLIP (C)', 'Leg 2 BSLIP (C)', 'NumColumns', 2);
%        hold off
       
   end
   
%    for i = 1:8
%        
%        figure
%        hold on
%        plot(timeVecV1, [qS.VPP.varyK.tau1Avg{i}], 'r', 'linewidth', 2)
%        plot(timeVecV1, [qS.VPP.varyK.tau2Avg{i}], 'b', 'linewidth', 2)
%        ylim([-1.5 1.5])
%        xlabel('Time []')
%        ylabel('Hip Torque [$\frac{\textrm{Nm}}{\textrm{kg}}$]')
%        legend('Lag Leg', 'Lead Leg', 'NumColumns', 2)
%        hold off
%        set(gca, 'FontSize', 14)
%        
%    end
   
%    for i = 1:8
%        
%        figure
%        hist(qS.VPP.varyK.errCOM(:, i))
%        title(['VPP Varying Trial ' num2str(i)])
%        
%        figure
%        hist(qS.VPP.constK.errCOM(:, i))
%        title(['VPP Constant Trial ' num2str(i)])
%        
%        figure
%        hist(qS.BSLIP.varyK.errCOM(:, i))
%        title(['BSLIP Varying Trial ' num2str(i)])
%        
%        figure
%        hist(qS.BSLIP.constK.errCOM(:, i))
%        title(['BSLIP Constant Trial ' num2str(i)])
%        
%    end
   %}
        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% Function Section %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function y = litcount(filename, literal)
% Search for number of string matches per line.

    fid = fopen(filename, 'rt');
    y = 0;

    while feof(fid) == 0

        tline = fgetl(fid);
        matches = contains(tline, literal);
        num = length(matches);

        if matches > 0

            y = y + num;
%             fprintf(1,'%d:%s\n',num,tline);

        end

    end

    fclose(fid);

end

function y = littime(filename, literal)
% Search for number of string matches per line.

    fid = fopen(filename, 'rt');
    y = 0;

    while feof(fid) == 0

        tline = fgetl(fid);
        
        if contains(tline, literal)
            
            y = y + str2double(tline(61:end));
            
        end

    end

    fclose(fid);

end
