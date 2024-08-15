%% %%%%%%%%%%%%%%%%%%%%%%%%%% postOptAnalysis %%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 27 June 2022
% Last Updated: 27 June 2022

% This code uses .mat files from optimized template models (see mainCRC.m)
% to analyze the resulting optimized models

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

close all
clearvars -except qS
clc

tempPath = fileparts(which('postOptAnalysis.m'));
addpath(genpath(tempPath));

tempPath = cd;
addpath([tempPath '\Data']);
addpath('C:\Users\dkell\Documents\MATLAB\CasADi');

import casadi.*

%% %%%%%%%%%%%%%%%%%%%%%%%%% Analysis Section %%%%%%%%%%%%%%%%%%%%%%%%%% %%

boolGraph = false;

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
subjIDEnd = 14;
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
cellK1 = cell(numSubjs, numTrials);
cellK2 = cell(numSubjs, numTrials);
cellT = cell(numSubjs, numTrials);
iter = 0;

for i = ((subjIDStart - 1)*8 + trialIDStart):((subjIDEnd - 1)*8 + trialIDEnd)
    
    load(matFiles{i});
    
    iter = iter + 1;
    indSubj = floor((iter-1)/numTrials) + 1;
    
    if ~mod(i,8)
        
        indTrial = 8;
        
    else
        
        indTrial = mod(i,8);
        
    end
    
    errCOM(indSubj, indTrial) = vS.quant.rmseCOMTimeND;
    errCOMH(indSubj, indTrial) = vS.quant.rmseCOMTimeHorND;
    errCOMV(indSubj, indTrial) = vS.quant.rmseCOMTimeVertND;
    errGRF(indSubj, indTrial) = mean(vS.quant.rmseGRFVertNDFull);
    errPhase(indSubj, indTrial) = vS.quant.rmsePhase;
    
    rangeTorso(indSubj, indTrial) = range(vS.optims.stateFull(:,3)*180/pi);
    
    springL1 = [];
    springL2 = springL1;

    for j = 1:(vS.params.steps/2)

       springL1 = [springL1; vS.optims.stateFull((4*(j-1)*100+1):(4*(j-1)*100+100), end-1);...
           zeros(100,1); vS.optims.stateFull((4*(j-1)*100+201):(4*(j-1)*100+400),end)];
       springL2 = [springL2; vS.optims.stateFull((4*(j-1)*100+1):(4*(j-1)*100+200), end);...
           vS.optims.stateFull((4*(j-1)*100+201):(4*(j-1)*100+300),end-1); zeros(100,1)];

    end
    
    cellK1(indSubj,indTrial) = {springL1*...
        (dS.subj.len0/dS.subj.bodyWeight)*1000};
    cellK2(indSubj,indTrial) = {springL2*...
        (dS.subj.len0/dS.subj.bodyWeight)*1000};
    
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
    
%     animTrajTime(vS, 0);
%     animForce(vS, dS, 0);

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
    
    if litcount(txtFiles{i}, 'Optimal Solution Found.')
        
        disp(['Subject ' num2str(floor((i-1)/8)+1) ' Trial '...
            num2str(mod(i,8)) ' Optimally Solved.'])
        count = count + 1;
        
    else
        
        disp(['Subject ' num2str(floor((i-1)/8)+1) ' Trial '...
            num2str(mod(i,8)) ' NOT OPTIMAL.'])
    
    end
    
end

disp(['Total Number of Trials Terminated at Optimal Solution: '...
    num2str(count)]);

errCOMAvg = mean(errCOM, 1);
errCOMSTD = std(errCOM, 0, 1);
rmTrialInd = zeros(size(errCOM));

for j = 1:length(errCOM(1,:))
    
    for i = 1:length(errCOM(:,1))
        
        if (errCOM(i, j) >= (errCOMAvg(j) + 2*errCOMSTD(j))) ||...
                (errCOM(i, j) <= (errCOMAvg(j) - 2*errCOMSTD(j)))
            
            errCOM(i, j) = NaN;
            rmTrialInd(i, j) = 1;
            
        end
            
    end
    
end

errCOMH(rmTrialInd == 1) = NaN;
errCOMV(rmTrialInd == 1) = NaN;
errGRF(rmTrialInd == 1) = NaN;
errPhase(rmTrialInd == 1) = NaN;
rangeTorso(rmTrialInd == 1) = NaN;
cellK1(rmTrialInd == 1) = {NaN(length(springL1), 1)};
cellK2(rmTrialInd == 1) = {NaN(length(springL1), 1)};
cellT(rmTrialInd == 1) = {NaN(1, 400)};

stiffAvg1 = cell(1, length(cellK1(1,:)));
stiffAvg2 = stiffAvg1;

stiffSTD1 = stiffAvg1;
stiffSTD2 = stiffAvg1;

timeAvg = stiffAvg1;
timeSTD = stiffAvg1;

for j = 1:length(cellK1(1,:))
    
    allK1 = reshape(cell2mat(cellK1(:,j)), 400, []);
    allK2 = reshape(cell2mat(cellK2(:,j)), 400, []);
    
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

if boolGraph
    
   xTicks = {'40\%', '55\%', '70\%', '85\%', '100\%', '115\%', '130\%', '145\%'};
   barX = categorical(xTicks);
   barX = reordercats(barX, xTicks);
   
   fieldQS1 = fieldnames(qS);
   numFields1 = numel(fieldQS1);
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
           
%            figure
%            hold on
%            bar(barX, barYCOM, 0.5, 'FaceColor', [0.6 0.9 0.1])
%            errorbar(barX, barYCOM, barYCOM2,...
%                'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%                'LineWidth', 1, 'CapSize', 10)
%            xlabel('Trial Speed (Based on Preferred Walking Speed)')
%            ylabel('COM Error Magnitude []')
%            legend(dataStr)
%            
%            figure
%            hold on
%            bar(barX, barYGRF, 0.5, 'FaceColor', [0.1 0.7 0.9])
%            errorbar(barX, barYGRF, barYGRF2,...
%                'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%                'LineWidth', 1, 'CapSize', 10)
%            xlabel('Trial Speed (Based on Preferred Walking Speed)')
%            ylabel('GRF Error Magnitude []')
%            legend(dataStr)
%            
%            figure
%            hold on
%            bar(barX, barYPhase, 0.5, 'FaceColor', [1.0 0.7 0.0])
%            errorbar(barX, barYPhase, barYPhase2,...
%                'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%                'LineWidth', 1, 'CapSize', 10)
%            xlabel('Trial Speed (Based on Preferred Walking Speed)')
%            ylabel('Phase Error Magnitude [sec]')
%            legend(dataStr)
           
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
               
           end
           
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
   
   figure
   hold on
   b1 = bar(barX, barYCOMFull);
%    errorbar(barX, barYCOMFull, barYCOM2Full,...
%        'LineStyle', 'none', 'Color', [0.0, 0.0, 0.0],...
%        'LineWidth', 1, 'CapSize', 10)
   xlabel('Trial Speed (Based on Preferred Walking Speed)')
   ylabel('COM Error Magnitude []')
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
   xlabel('Trial Speed (Based on Preferred Walking Speed)')
   ylabel('GRF Error Magnitude []')
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
   xlabel('Trial Speed (Based on Preferred Walking Speed)')
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

        if num > 0

            y = y + num;
%             fprintf(1,'%d:%s\n',num,tline);

        end

    end

    fclose(fid);

end
