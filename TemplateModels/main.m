%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 20 June 2022
% Last Updated: 15 November 2023

% This code uses the CasADi framework to model a walking gait based on
% the Virtual Pivot Point (VPP) or Bipedal Spring Loaded Inverted Pendulum
% (BSLIP) template.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

close all
%clear variables
clc

% Path for implementing CasADi framework
% NOTE: MUST HAVE CasADi PATH ESTABLISHED PRIOR TO RUNNING CODE
% EXAMPLES:
% addpath('C:\Users\username\Documents\MATLAB\CasADi');
import casadi.*

% Path for accessing all required files
% NOTE: THIS METHOD ADDS PATHS TO WHERE main.m AND THE HUMAN DATA RESIDES
% tempPath = fileparts(which('main.m'));
% addpath(genpath(tempPath))
% tempPath = cd;
% tempPath = fullfile(tempPath, '..');
% addpath([tempPath '/HumanData'])

% NOTE: THIS METHOD ADDS THE REQUIRED PATHS EXPLICITLY
% addpath(genpath(['C:\Users\username\Documents\MATLAB\templateModels']));
% addpath(['C:\Users\username\Documents\MATLAB\HumanData']);

% Initialize subjList struct
subjList.header = 'MoCap Human Data';

% Initiliaze subject-specific information from human data files
% Data is for the 24 young healthy subjects from WBDSinfo file
subjWeight = [74.3, 52.9, 75.85, 61.05, 77.55, 83.15, 71.75, 64, 61.25,...
    61.7, 77.55, 48.8, 95.4, 64.95, 89.3, 77.9, 79.05, 59.3, 66.25,...
    62.45, 61.5, 61.15, 44.9, 69.55]; % Subject weight
subjLeg = [0.89, 0.865, 0.94, 0.94, 0.86, 0.91, 0.84, 0.97, 0.84, 0.87,...
    0.98, 0.78, 0.85, 0.8, 0.92, 0.88, 0.98, 0.86, 0.89, 0.88, 0.91,...
    0.87, 0.81, 0.75]; % Subject leg length
subjTread = {[0.49, 0.67, 0.85, 1.03, 1.21, 1.4, 1.58, 1.76],...
    [0.5, 0.69, 0.88, 1.06, 1.25, 1.44, 1.63, 1.81],...
    [0.39, 0.54, 0.68, 0.83, 0.98, 1.12, 1.27, 1.42],...
    [0.52, 0.71, 0.91, 1.1, 1.3, 1.49, 1.69, 1.88],...
    [0.52, 0.71, 0.9, 1.1, 1.29, 1.48, 1.68, 1.87],...
    [0.51, 0.71, 0.9, 1.09, 1.28, 1.48, 1.67, 1.86],...
    [0.44, 0.61, 0.77, 0.94, 1.1, 1.27, 1.43, 1.6],...
    [0.55, 0.76, 0.97, 1.17, 1.38, 1.59, 1.79, 2],...
    [0.47, 0.64, 0.82, 0.99, 1.17, 1.35, 1.52, 1.7],...
    [0.55, 0.76, 0.96, 1.17, 1.37, 1.58, 1.79, 1.99],...
    [0.53, 0.72, 0.92, 1.12, 1.31, 1.51, 1.71, 1.9],...
    [0.42, 0.58, 0.74, 0.9, 1.05, 1.21, 1.37, 1.53],...
    [0.41, 0.56, 0.71, 0.87, 1.02, 1.17, 1.33, 1.48],...
    [0.45, 0.62, 0.78, 0.95, 1.12, 1.29, 1.46, 1.62],...
    [0.61, 0.85, 1.08, 1.31, 1.54, 1.77, 2, 2.23],...
    [0.44, 0.61, 0.77, 0.94, 1.1, 1.27, 1.43, 1.6],...
    [0.55, 0.75, 0.96, 1.16, 1.37, 1.57, 1.78, 1.98],...
    [0.5, 0.69, 0.88, 1.07, 1.26, 1.45, 1.64, 1.83],...
    [0.48, 0.66, 0.84, 1.03, 1.21, 1.39, 1.57, 1.75],...
    [0.51, 0.7, 0.9, 1.09, 1.28, 1.47, 1.66, 1.86],...
    [0.57, 0.78, 0.99, 1.21, 1.42, 1.63, 1.85, 2.06],...
    [0.53, 0.73, 0.93, 1.12, 1.32, 1.52, 1.72, 1.92],...
    [0.61, 0.84, 1.06, 1.29, 1.52, 1.75, 1.98, 2.2],...
    [0.41, 0.56, 0.71, 0.86, 1.02, 1.17, 1.2, 1.47]}; % Subject treadmill speeds

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

tempModel = 'VPP'; % Model type [VPP or BSLIP]
tempIntMethod = 'Collocation'; % Propagation method [Collocation or Multishooting]
tempStiffMethod = 'Varying'; % Leg stiffness [Varying or Constant]
tempVPPMethod = 'Constant'; % VPP position [Varying or Constant]
tempHum = 'MoCap'; % Human data source [MoCap]

tempN = 25; % Number of elements per phase (Shooting or Finite)
tempM = 4; % Number of points per element (RK4 or Collocation)

tempBoolGaitAvg = true; % Flag for using averaged gait data
tempBoolGaitAppend = false; % Flag for appending gait cycles (averaged only)

tempRuns = 1; % Number of times to run the optimization

tempQ = [1 0; 0 1]; % Object cost weights

tempVDiff = 0; % Legacy code variable [not used]
tempEps = 0.5; % Legacy code variable [still used]

% Create counter variables for input control variable creation
tempInputCountD = 1;
tempInputCountS = 1;
tempShootCount = 1;

% Create VPP-specific parameters
tempGamma = 0; % Angle offset from torso orientation [rad] [not in use]
tempRH = 0.1; % Offset of COM from hip along torso [m]
tempRVPP = [0.0, 0.3]; % VPP distance from COM [m] [DS,SS]
tempJ = 4.58; % Mass inertia [kg m^2]

% Boolean flags
tempRun1 = true; % Set if this is the first optimization run (leave true)
tempBoolPlot = true; % Plot COM/GRF data from humanDataX.m
tempBoolSave = false; % Clear tempX variables for saving human data to file
tempBoolSeg = false; % Legacy code (not in use)
tempBoolClear = false; % Clear ALL tempX variables after optimization
tempBoolSaveData = true; % Save data post-optimization

% Initialize time tracking variable
tempTimeTrack = 0;

% Run optimizer for the following \% of preferred walking speed (PWS)
for j = 1:8 % {1:40% 2:55% 3:70% 4:85% 5:100% 6:115% 7:130% 8:145%} PWS
    
    % Run optimizer for the following subjects
    for i = 1:24 % Vector of integers between 1-24
        
        % Store appropriate subject parameters in subjList struct
        subjList.(['Subj' num2str(i)]).weight = subjWeight(i);
        subjList.(['Subj' num2str(i)]).tread = subjTread{i};
        subjList.(['Subj' num2str(i)]).leg = subjLeg(i);
        subjList.(['Subj' num2str(i)]).trialNum = j;
        
        % Set subject ID parameter appropriately
        if (i < 10)
            
            tempSubjectID = ['0' num2str(i)]; % Subject ID [0X]
            
        else
            
            tempSubjectID = num2str(i);
            
        end
        
        % Set up trial-specific parameters
        tempSubjectTrial = ['0' num2str(subjList.(['Subj' num2str(i)]).trialNum)]; % Subject trial # [0X]
        tempIndStart = 2; % Trial data gait cycle starting index
        tempIndEnd = 25; % Trial data gait cycle ending index
        tempTreadmill = subjList.(['Subj' num2str(i)]).tread(...
            subjList.(['Subj' num2str(i)]).trialNum); % Trial treadmill speed [m/s]
        tempWeight = subjList.(['Subj' num2str(i)]).weight; % Subject mass [kg]
        tempLen0 = subjList.(['Subj' num2str(i)]).leg; % Subject leg length [m]
        tempFitType = 1; % Human data fitting method {1: fourier 2: polyval}
        
        % Stack parameters for human data retrieval
        tempParamsData = {tempSubjectID, tempSubjectTrial, tempIndStart,...
            tempIndEnd, tempTreadmill, tempWeight, tempLen0, tempBoolGaitAvg,...
            tempBoolGaitAppend, tempBoolPlot, tempBoolSeg, tempBoolSave};
        
        % Retrieve human data
        dataStruct = humanData(tempParamsData);
        
        % Save dataStruct struct to file
        % EXAMPLE:
        % save(['C:\Users\username\Documents\Data\S' tempSubjectID '_T' tempSubjectTrial '_DS.mat'], 'dataStruct');
        
        continue
        
        % Set the number of steps to consider [2x number of gait cycles]
        if tempBoolGaitAvg
            
            tempSteps = 2;
            
        else
            
            tempSteps = 2*dataStruct.inds.range;
            
        end
        
        % Stack parameters for varStruct struct initialization
        tempParamsVar = {tempSteps, tempGamma, tempRH, tempRVPP, tempJ,...
            tempVDiff, tempEps, tempN, tempM, tempRun1,...
            tempInputCountD, tempInputCountS, tempShootCount, tempTimeTrack,...
            tempModel, tempIntMethod, tempHum, tempStiffMethod, tempVPPMethod,...
            tempFitType, tempQ};
        
        % Stack parameters for optimization setup
        tempParamsSetup = {tempHum, tempRuns, tempBoolClear};
        
        % Call into optimization framework
        templateModelOpt(dataStruct, tempParamsSetup, tempParamsVar);
        
    end
    
end