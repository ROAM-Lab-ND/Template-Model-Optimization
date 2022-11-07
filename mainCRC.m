%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 3 June 2022
% Last Updated: 3 June 2022

% This code uses the CasADi framework to model a walking gait based on
% the Virtual Pivot Point (VPP) or Bipedal Spring Loaded Inverted Pendulum 
% (BSLIP) template.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Clear and close any existing data
close all
clear variables
clc

% Path for Implementing CasADi Framework
% NOTE: MUST HAVE CasADi PATH ESTABLISHED PRIOR TO RUNNING CODE
% addpath(['/afs/crc.nd.edu/user/d/dkelly7/Private/'...
%     'casadi-linux-matlabR2014b-v3.5.5'])

addpath('/afs/crc.nd.edu/user/d/dkelly7/Private/CasADi')
addpath(genpath('/afs/crc.nd.edu/user/d/dkelly7/Private/Matlab'))

% Path for Accessing All Required Files
% NOTE: THIS PATH SHOULD BE SET TO THE PARENT DIRECTORY OF THE HUMAN DATA
% AND MATLAB SCRIPTS
% tempPath = fileparts(which('mainCRC.m'));
% addpath(genpath(tempPath))
% tempPath = cd;
% tempPath = fullfile(tempPath, '..');
% addpath([tempPath '/HumanData'])

% Import CasADi toolbox
import casadi.*

% Young subject data values from mocap
subjList.header = 'MoCap Human Data';
subjWeight = [74.3, 52.9, 75.85, 61.05, 77.55, 83.15, 71.75, 64, 61.25,...
    61.7, 77.55, 48.8, 95.4, 64.95, 89.3, 77.9, 79.05, 59.3, 66.25,...
    62.45, 61.5, 61.15, 44.9, 69.55];
subjLeg = [0.89, 0.865, 0.94, 0.94, 0.86, 0.91, 0.84, 0.97, 0.84, 0.87,...
    0.98, 0.78, 0.85, 0.8, 0.92, 0.88, 0.98, 0.86, 0.89, 0.88, 0.91,...
    0.87, 0.81, 0.75];
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
    [0.41, 0.56, 0.71, 0.86, 1.02, 1.17, 1.2, 1.47]};

% trialNum = input(['Which trial number are you interested in?'...
%     ' \n1) 40% Preferred Walking Speed'...
%     ' \n2) 55% Preferred Walking Speed'...
%     ' \n3) 70% Below Preferred Walking Speed'...
%     ' \n4) 85% Preferred Walking Speed'...
%     ' \n5) Preferred Walking Speed'...
%     ' \n6) 115% Preferred Walking Speed'...
%     ' \n7) 130% Preferred Walking Speed'...
%     ' \n8) 145% Preferred Walking Speed'...
%     ' \n0) Exit Without Solving'...
%     ' \nChoice: ']);

trialNum = 5;

if (~trialNum)
    
    disp('Exit Status Called.')
    return
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Template model parameters to set
tempModel = 'BSLIP'; % Model Type [VPP or BSLIP]
tempIntMethod = 'Collocation'; % Int Method [Collocation or Multishooting]
tempStiffMethod = 'Constant'; % Leg Stiffness Option [Varying or Constant]
tempVPPMethod = 'Constant'; % VPP Location Option [Varying or Constant]
tempHum = 'MoCap'; % Human Data Source [MoCap or XSENS]

% Shooting and finite elements
tempN = 25; % Shooting
tempM = 4; % Finite

% Flags for gait cycle creation
tempBoolGaitAvg = true; % Flag for Using Averaged Gait Data
tempBoolGaitAppend = false; % Flag for Appending Gait Cycles (AVG ONLY)

% Preferred number of times to run optimization
tempRuns = 1; % Number of times to run the optimization

% Penalty gains for components of objective cost
tempQ = [1 0;
    0 1];

% Potentially deprecated parameters, don't touch for now
tempVDiff = 0;
tempEps = 0.5;

% Create iterators for each phase and overall cycle
tempInputCountD = 1;
tempInputCountS = 1;
tempShootCount = 1;

% VPP specific parameters
tempGamma = 0; % VP-Torso Orientation [rad]
tempRH = 0.1; % CoM-Hip offset [m]
tempRVPP = [0, 0.2]; % VP-CoM offset [m] [DS,SS]
tempJ = 4.58; % Trunk inertia [kg m^2]

% Flags for code options
tempRun1 = true; % First optimization flag
tempBoolPlot = false; % Plot subject data flag
tempBoolSave = false; % Clear temporary variables for saving purposes flag
tempBoolSeg = false; % Segment data collection flag [XSENS Only]
tempBoolClear = false; % Clear all temporary variables flag
tempBoolSaveData = true; % Save optimization result flag

% Initialize time tracking variable
tempTimeTrack = 0;

% For each trial speed
for j = 1:8
    
    % For each young healthy subject
    for i = 1:24

        % Create struct for current subject parameters
        subjList.(['Subj' num2str(i)]).weight = subjWeight(i); % Weight
        subjList.(['Subj' num2str(i)]).tread = subjTread{i}; % Trial speed
        subjList.(['Subj' num2str(i)]).leg = subjLeg(i); % Leg length
        subjList.(['Subj' num2str(i)]).trialNum = j; % Trial number
        
        % Check if subject number is less than 10 or not
        if (i < 10)
            
            % Create subject ID
            tempSubjectID = ['0' num2str(i)];
            
        else
            
            % Create subject ID
            tempSubjectID = num2str(i);
            
        end
        
        % Grab subject specific data
        tempSubjectTrial = ['0'...
            num2str(subjList.(['Subj' num2str(i)]).trialNum)];
        tempIndStart = 2; % Trial Gait Cycle Starting Index
        tempIndEnd = 25; % Trial Gait Cycle Ending Index
        tempTreadmill = subjList.(['Subj' num2str(i)]).tread(...
            subjList.(['Subj' num2str(i)]).trialNum);
        tempWeight = subjList.(['Subj' num2str(i)]).weight;
        tempLen0 = subjList.(['Subj' num2str(i)]).leg;
        
        % Choose Fourier fit or polyfit for subject data [1, 2]
        tempCheck = 1;
        
        % Stack parameters for subject data retrieval
        tempParamsData = {tempSubjectID, tempSubjectTrial, tempIndStart,...
            tempIndEnd, tempTreadmill, tempWeight, tempLen0,...
            tempBoolGaitAvg, tempBoolGaitAppend, tempBoolPlot,...
            tempBoolSeg, tempBoolSave};
        
        % Check source of human data
        switch tempHum
            
            case 'MoCap'
                
                % Call to mocap data script
                dataStruct = humanDataMoCap(tempParamsData);
                
            case 'XSENS'
                
                % Call to XSENS data script
                dataStruct = humanDataXSENS(tempParamsData);
                
            otherwise
                
                tempMsg = 'Invalid data source chosen';
                error(tempMsg);
                
        end
        
        % Check if gait cycle is averaged or not
        if tempBoolGaitAvg
            
            % Set number of steps to be taken by template model
            tempSteps = 2;
            
        else
            
            % Set number of steps to be taken by template model
            tempSteps = 2*dataStruct.inds.range;
            
        end
        
        % Stack parameters for structure initialization
        tempParamsVar = {tempSteps, tempGamma, tempRH, tempRVPP, tempJ,...
            tempVDiff, tempEps, tempN, tempM, tempRun1,...
            tempInputCountD, tempInputCountS, tempShootCount,...
            tempTimeTrack, tempModel, tempIntMethod, tempHum,...
            tempStiffMethod, tempVPPMethod, tempCheck, tempQ};
        
        % Stack parameters for optimization configuration
        tempParamsSetup = {tempHum, tempRuns, tempBoolClear};
        
        % Call into optimization framework
        templateModelOptCRC(dataStruct, tempParamsSetup, tempParamsVar);

    end
    
end