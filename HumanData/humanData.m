function mocap = humanDataMoCap(params)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% humanData %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    % Author: David Kelly
    % Created: 3 November 2021
    % Last Updated: 3 November 2021

    % This code uses the motion capture data from Roopak's data set to 
    % parameterize the trajectory of the COM over a gait cycle.  This will be 
    % determined by the contact events, starting at contact of the heel of one 
    % of the legs and concluding at the next instant of heel contact with the 
    % ground (2 steps).

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%% %%
    
    format short
    addpath(['C:\Users\dkell\Documents\MATLAB\Research\OSL_Modeling\'...
        'OSL_Templates\HumanData\WBDSascii']);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% Manual Info %%%%%%%%%%%%%%%%%%%%%%%%%% %%

    % Subject Information (based on excel sheet)
    mocap.subj.subjectID = num2str(params{1}); % Subject ID
    mocap.subj.testNum = num2str(params{2}); % Trial Number [0 to 8]
    mocap.subj.date = date; % Date String
    mocap.subj.fileName = ['subject' mocap.subj.subjectID 'test'...
        mocap.subj.testNum date '.mat']; % Filename for Savin to .mat
    mocap.subj.treadmill = params{5}; % Treadmill Speed [m/s]
    mocap.subj.weight = params{6}; % Subject Mass [kg]
    mocap.subj.len0 = params{7}; % Subject Measured Leg Length [m]
    mocap.subj.g = 9.81; % Gravity Constant [m/s^2]
    mocap.subj.bodyWeight = mocap.subj.weight*mocap.subj.g; % Subject Weight [N]
    mocap.subj.gaitAvg = params{8};
    mocap.subj.gaitAppend = params{9};

    % Indices for Gait Cycle Searching
    mocap.inds.start = params{3};
    mocap.inds.end = params{4};
    mocap.inds.range = params{4} - params{3};
    tempGaitNum = params{4} - params{3};

    % File Load For Markers and GRF text files
    markers = importdata(['WBDS' mocap.subj.subjectID 'walkT'...
        mocap.subj.testNum 'mkr.txt']);
    force = importdata(['WBDS' mocap.subj.subjectID 'walkT'...
        mocap.subj.testNum 'grf.txt']);

    % Polynomial Size
    mocap.fit.polyPosHor = 20;
    
    if params{8}
    
        mocap.fit.polyPosVert = 25;
        
    else
        
        mocap.fit.polyPosVert = 25;
        
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%% Data Parsing %%%%%%%%%%%%%%%%%%%%%%%%%% %%

    % Specify Labels to Look for COM and GRF data
    tempLabelsCOM = {'Time','R.ASISX','R.ASISY','R.ASISZ','L.ASISX','L.ASISY',...
        'L.ASISZ','R.PSISX','R.PSISY','R.PSISZ','L.PSISX','L.PSISY',...
        'L.PSISZ','R.HeelY','L.HeelY','R.MT1X','R.MT1Y','R.MT1Z','L.MT1X',...
        'L.MT1Y','L.MT1Z','R.MT5X','R.MT5Y','R.MT5Z','L.MT5X','L.MT5Y',...
        'L.MT5Z'};
    tempLabelsGRF = {'Time','Fy1','Fy2','Fx1','Fx2','Fz1','Fz2'};

    % Find Columns that Correspond to Above Labels
    tempindicesCOM = indexFind(markers.colheaders,tempLabelsCOM);
    tempindicesGRF = indexFind(force.colheaders,tempLabelsGRF);

    % Grab Time Data from MoCap
    mocap.data.timeCOM = markers.data(:,tempindicesCOM(1));

    % Grab COM Position Data
    mocap.data.COM = [centroid(markers.data,tempindicesCOM(1,[2 5 8 11])),...
        centroid(markers.data,tempindicesCOM(1,[3 6 9 12])),...
        centroid(markers.data,tempindicesCOM(1,[4 7 10 13]))];

    % Grab Heel Vertical Data, Apply Moving Average to Smooth Data
    mocap.data.heels = markers.data(:,tempindicesCOM([14,15]))/1000; %[right, left]
    mocap.data.heels = movmean(mocap.data.heels,10);

    % Grab Metatarsal Vertical Data
    mocap.data.metRight = [centroid(markers.data,tempindicesCOM(1,[16,22])),...
        centroid(markers.data,tempindicesCOM(1,[17,23])),...
        centroid(markers.data,tempindicesCOM(1,[18,24]))];

    mocap.data.metLeft = [centroid(markers.data,tempindicesCOM(1,[19,25])),...
        centroid(markers.data,tempindicesCOM(1,[20,26])),...
        centroid(markers.data,tempindicesCOM(1,[21,27]))];

    % Grab GRF Fore/Aft and Vertical Data, Apply Moving Average to Smooth
    mocap.data.grf.hor = [force.data(:,tempindicesGRF(4)),...
        force.data(:,tempindicesGRF(5))];
    mocap.data.grf.hor = [movmean(mocap.data.grf.hor(:,1),10),...
        movmean(mocap.data.grf.hor(:,2),10)];
    
    mocap.data.grf.vert = [force.data(:,tempindicesGRF(2)),...
        force.data(:,tempindicesGRF(3))];
    mocap.data.grf.vert = [movmean(mocap.data.grf.vert(:,1),20),...
        movmean(mocap.data.grf.vert(:,2),20)];
    
    % Set Negative Vertical GRF to Zero
    mocap.data.grf.vert(find(mocap.data.grf.vert(:,1)<0),1)=0;
    mocap.data.grf.vert(find(mocap.data.grf.vert(:,2)<0),2)=0;

    % Normalize GRF Fore/Aft and Vertical Data
    mocap.data.grf.horNorm = mocap.data.grf.vert/mocap.subj.bodyWeight;
    
    mocap.data.grf.vertNorm = mocap.data.grf.vert/mocap.subj.bodyWeight;

    % Create Time Vector Based on COM Time and GRF Length
    mocap.data.timeGRF = linspace(mocap.data.timeCOM(1),...
        mocap.data.timeCOM(end),length(mocap.data.grf.vert(:,1)));

    %% %%%%%%%%%%%%%%%%%%%%% Contact Events Search %%%%%%%%%%%%%%%%%%%%% %%
    
    % Liftoff and Touchdown Indices Based on GRF1Norm
%     mocap.inds.LO1 = find((mocap.data.grf.vertNorm(1:end-1,1)>0.01 &...
%         mocap.data.grf.vertNorm(2:end,1)<=0.01));
    
    mocap.inds.LO1 = find((mocap.data.grf.vertNorm(1:end-2,1)>0.01 &...
        mocap.data.grf.vertNorm(2:end-1,1)<=0.01 &...
        mocap.data.grf.vertNorm(3:end,1)<=0.01));
%     mocap.inds.LO1((find(diff(mocap.inds.LO1)<50))+1) = [];
    
%     mocap.inds.TD1 = find((mocap.data.grf.vertNorm(1:end-1,1)<=0.01 &...
%         mocap.data.grf.vertNorm(2:end,1)>0.01));
    
    mocap.inds.TD1 = find((mocap.data.grf.vertNorm(1:end-2,1)<=0.01 &...
        mocap.data.grf.vertNorm(2:end-1,1)>0.01 &...
        mocap.data.grf.vertNorm(3:end,1)>0.01));
%     mocap.inds.TD1((find(diff(mocap.inds.TD1)<50))+1) = [];

    % Liftoff and Touchdown Indices Based on GRF2Norm
%     mocap.inds.LO2 = find((mocap.data.grf.vertNorm(1:end-1,2)>0.01 &...
%         mocap.data.grf.vertNorm(2:end,2)<=0.01));
    
    mocap.inds.LO2 = find((mocap.data.grf.vertNorm(1:end-2,2)>0.01 &...
        mocap.data.grf.vertNorm(2:end-1,2)<=0.01 &...
        mocap.data.grf.vertNorm(3:end,2)<=0.01));
%     mocap.inds.LO2((find(diff(mocap.inds.LO2)<50))+1) = [];

%     mocap.inds.TD2 = find((mocap.data.grf.vertNorm(1:end-1,2)<=0.01 &...
%         mocap.data.grf.vertNorm(2:end,2)>0.01));
    
    mocap.inds.TD2 = find((mocap.data.grf.vertNorm(1:end-2,2)<=0.01 &...
        mocap.data.grf.vertNorm(2:end-1,2)>0.01 &...
        mocap.data.grf.vertNorm(3:end,2)>0.01));
%     mocap.inds.TD2((find(diff(mocap.inds.TD2)<50))+1) = [];

    % Shift Odd Indices to Even (GRF Data Freq = 2*COM Data Freq) 
    mocap.inds.LO1 = mocap.inds.LO1 + mod(mocap.inds.LO1,2);
    mocap.inds.LO2 = mocap.inds.LO2 + mod(mocap.inds.LO2,2);
    mocap.inds.TD1 = mocap.inds.TD1 - mod(mocap.inds.TD1,2);
    mocap.inds.TD2 = mocap.inds.TD2 - mod(mocap.inds.TD2,2);
    
    % Temporary Variable to Align Contact Events
    tempLenL01 = length(mocap.inds.LO1);
    tempLenL02 = length(mocap.inds.LO2);
    tempLenTD1 = length(mocap.inds.TD1);
    tempLenTD2 = length(mocap.inds.TD2);
    
    tempIndsSize = min([tempLenL01, tempLenL02, tempLenTD1, tempLenTD2]);
    
    % Resize Indices for Standardized Check
    mocap.inds.LO1 = mocap.inds.LO1(1:tempIndsSize);
    mocap.inds.LO2 = mocap.inds.LO2(1:tempIndsSize);
    mocap.inds.TD1 = mocap.inds.TD1(1:tempIndsSize);
    mocap.inds.TD2 = mocap.inds.TD2(1:tempIndsSize);
    
    mocap.inds.LO1 = mocap.inds.LO1(...
        find(mocap.inds.LO1>mocap.inds.TD1(1)):end);
    mocap.inds.LO2 = mocap.inds.LO2(...
        find(mocap.inds.LO2>mocap.inds.TD1(1)):end);
    mocap.inds.TD2 = mocap.inds.TD2(...
        find(mocap.inds.TD2>mocap.inds.TD1(1)):end);
    
    if mocap.inds.range > (length(mocap.inds.TD1) - mocap.inds.start - 1)
        
        mocap.inds.range = length(mocap.inds.TD1) - mocap.inds.start - 1;
        mocap.inds.end = mocap.inds.start + mocap.inds.range;
        
    end
    
    % Initialize Arrays to Store Phase Durations for Each Gait Cycle
    tempTimeDS1 = NaN(mocap.inds.range, 1);
    tempTimeSS1 = tempTimeDS1;
    tempTimeDS2 = tempTimeDS1;
    tempTimeSS2 = tempTimeDS1;
    
    % Calculate Phase Durations for Each Gait Cycle
    for i = 0:(mocap.inds.range - 1)
        
        tempTimeDS1(i+1) = mocap.data.timeCOM(...
            mocap.inds.LO2(mocap.inds.start + i)/2) -...
            mocap.data.timeCOM(mocap.inds.TD1(mocap.inds.start + i)/2);
        tempTimeSS1(i+1) = mocap.data.timeCOM(...
            mocap.inds.TD2(mocap.inds.start + i)/2) -...
            mocap.data.timeCOM(mocap.inds.LO2(mocap.inds.start + i)/2);
        tempTimeDS2(i+1) = mocap.data.timeCOM(...
            mocap.inds.LO1(mocap.inds.start + i)/2) -...
            mocap.data.timeCOM(mocap.inds.TD2(mocap.inds.start + i)/2);
        tempTimeSS2(i+1) = mocap.data.timeCOM(...
            mocap.inds.TD1(mocap.inds.start + i + 1)/2) -...
            mocap.data.timeCOM(mocap.inds.LO1(mocap.inds.start + i)/2);
        
    end
    
    % Check if Average Gait Data is Flagged
    if params{8}

        % Calculate Average Duration for Each Support Phase
        mocap.data.timeDS1 = mean(tempTimeDS1);
        mocap.data.timeSS1 = mean(tempTimeSS1);
        mocap.data.timeDS2 = mean(tempTimeDS2);
        mocap.data.timeSS2 = mean(tempTimeSS2);
        
    else
        
        % Store Duration for Each Support Phase
        mocap.data.timeDS1 = tempTimeDS1;
        mocap.data.timeSS1 = tempTimeSS1;
        mocap.data.timeDS2 = tempTimeDS2;
        mocap.data.timeSS2 = tempTimeSS2;
        
    end
    
    % Store 
    mocap.data.timeFull = [mocap.data.timeDS1'; mocap.data.timeSS1';...
        mocap.data.timeDS2'; mocap.data.timeSS2'];
    mocap.data.timeFull = mocap.data.timeFull(:);

    mocap.data.timeGait = sum(mocap.data.timeDS1 + mocap.data.timeSS1 +...
        mocap.data.timeDS2 + mocap.data.timeSS2);
    
    tempSkipCount = 0;
    tempSkipInd = [];
    
    % Check if Average Gait Data is Flagged
    if params{8}
        
        % Initialize Cells to Store Gait Times for COM/GRF (Full + Norm)
        tempTimeCOM = cell(1, mocap.inds.range);
        tempTimeCOM2 = tempTimeCOM;
        tempTimeCOMNorm = tempTimeCOM;
        tempTimeGRF = tempTimeCOM;
        tempTimeGRFNorm = tempTimeCOM;
        
        % Initialize Cells to Store COM/GRF Fore/Aft and Vertical Gait Data
        tempTrajCOMHor = tempTimeCOM;
        tempTrajCOMHorF1 = tempTimeCOM;
        tempTrajCOMHorF2 = tempTimeCOM;
        tempTrajCOMMedLat = tempTimeCOM;
        tempTrajCOMVert = tempTimeCOM;
        tempTrajGRFHor = tempTimeCOM;
        tempTrajGRFVert = tempTimeCOM;
        
        % Initialize Arrays to Store Interpolation for Time, COM, and GRF
        tempTime = 0:0.001:1;
        tempTrajCOMHor2 = NaN(mocap.inds.range, length(tempTime));
        tempTrajCOMHor2F1 = NaN(mocap.inds.range, length(tempTime));
        tempTrajCOMHor2F2 = NaN(mocap.inds.range, length(tempTime));
        tempTrajCOMMedLat2 = NaN(mocap.inds.range, length(tempTime));
        tempTrajCOMVert2 = NaN(mocap.inds.range, length(tempTime));
        tempTrajGRFHor1 = NaN(mocap.inds.range, length(tempTime));
        tempTrajGRFVert1 = NaN(mocap.inds.range, length(tempTime));
        tempTrajGRFHor2 = NaN(mocap.inds.range, length(tempTime));
        tempTrajGRFVert2 = NaN(mocap.inds.range, length(tempTime));
        
        for i = 0:(mocap.inds.range-1)
            
            if isnan(mocap.data.COM(...
                    mocap.inds.TD1(mocap.inds.start + i)/2, 1))
               
                tempSkipCount = tempSkipCount + 1;
                tempSkipInd = [tempSkipInd, i+1];
                
                continue
                
            end
            
            % Grab Time Arrays for Each Gait Cycle, Calculate Norm for Each
                % COM
            tempTimeCOM{i+1} = mocap.data.timeCOM(mocap.inds.TD1(...
                (mocap.inds.start + i))/2:mocap.inds.TD1(...
                (mocap.inds.start + i + 1))/2) -...
                mocap.data.timeCOM(mocap.inds.TD1(...
                (mocap.inds.start + i))/2);
            
            tempTimeCOM2{i+1} = mocap.data.timeCOM(mocap.inds.TD1(...
                (mocap.inds.start + i))/2:mocap.inds.TD1(...
                (mocap.inds.start + i + 1))/2);
            
            tempTimeCOMNorm{i+1} = tempTimeCOM{i+1}/...
                (mocap.data.timeCOM(mocap.inds.TD1(...
                (mocap.inds.start + i + 1))/2) -...
                mocap.data.timeCOM(mocap.inds.TD1(...
                (mocap.inds.start + i))/2));
                
                % GRF
            tempTimeGRF{i+1} = mocap.data.timeGRF(mocap.inds.TD1(...
                (mocap.inds.start + i)):mocap.inds.TD1(...
                (mocap.inds.start + i + 1))) -...
                mocap.data.timeGRF(mocap.inds.TD1((mocap.inds.start + i)));
            
            tempTimeGRFNorm{i+1} = tempTimeGRF{i+1}/...
                (mocap.data.timeGRF(mocap.inds.TD1(...
                (mocap.inds.start + i + 1))) -...
                mocap.data.timeGRF(mocap.inds.TD1(...
                (mocap.inds.start + i))));        
            
            % Grab Fore/Aft and Vertical COM Trajectory for Each Gait Cycle
            tempTrajCOMHor{i+1} = (mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2:...
                mocap.inds.TD1(mocap.inds.start + i + 1)/2, 1) -...
                mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2, 1))/1000 +...
                mocap.subj.treadmill*tempTimeCOM{i+1};
            
            tempTrajCOMHorF1{i+1} = (mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2:...
                mocap.inds.TD1(mocap.inds.start + i + 1)/2, 1) -...
                mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2, 1))/1000;
            
            tempTrajCOMHorF2{i+1} = mocap.subj.treadmill*tempTimeCOM{i+1};
            
            %{
            tempTrajCOMMedLat{i+1} = (mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2:...
                mocap.inds.TD1(mocap.inds.start + i + 1)/2, 3) -...
                mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2, 3))/1000;
            %}
            tempTrajCOMMedLat{i+1} = (mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2:...
                mocap.inds.TD1(mocap.inds.start + i + 1)/2, 3))/1000;
            
            tempTrajCOMVert{i+1} = mocap.data.COM(mocap.inds.TD1(...
                mocap.inds.start + i)/2:...
                mocap.inds.TD1(mocap.inds.start + i + 1)/2, 2)/1000;            
            
            % Grab Fore/Aft and Vertical GRF Profile for Each Gait Cycle
            tempTrajGRFHor{i+1} = mocap.data.grf.hor(...
                mocap.inds.TD1(mocap.inds.start + i):...
                mocap.inds.TD1(mocap.inds.start + i + 1),:);
            
            tempTrajGRFVert{i+1} = mocap.data.grf.vert(...
                mocap.inds.TD1(mocap.inds.start + i):...
                mocap.inds.TD1(mocap.inds.start + i + 1),:);         
            
            % Use Spline Interpolation to Align COM Data Points
            tempTrajCOMHor2(i+1,:) = interp1(tempTimeCOMNorm{i+1},...
                tempTrajCOMHor{i+1}, tempTime, 'spline');
            tempTrajCOMMedLat2(i+1,:) = interp1(tempTimeCOMNorm{i+1},...
                tempTrajCOMMedLat{i+1}, tempTime, 'spline');
            tempTrajCOMVert2(i+1,:) = interp1(tempTimeCOMNorm{i+1},...
                tempTrajCOMVert{i+1}, tempTime, 'spline');
            
            tempTrajCOMHor2F1(i+1,:) = interp1(tempTimeCOMNorm{i+1},...
                tempTrajCOMHorF1{i+1}, tempTime, 'spline');
            tempTrajCOMHor2F2(i+1,:) = interp1(tempTimeCOMNorm{i+1},...
                tempTrajCOMHorF2{i+1}, tempTime, 'spline');
            
            % Use Spline Interpolation to Align GRF Data Points
            tempTrajGRFHor1(i+1,:) = interp1(tempTimeGRFNorm{i+1},...
                tempTrajGRFHor{i+1}(:,1)', tempTime, 'spline');
            tempTrajGRFVert1(i+1,:) = interp1(tempTimeGRFNorm{i+1},...
                tempTrajGRFVert{i+1}(:,1)', tempTime, 'spline');
            tempTrajGRFHor2(i+1,:) = interp1(tempTimeGRFNorm{i+1},...
                tempTrajGRFHor{i+1}(:,2)', tempTime, 'spline');
            tempTrajGRFVert2(i+1,:) = interp1(tempTimeGRFNorm{i+1},...
                tempTrajGRFVert{i+1}(:,2)', tempTime, 'spline');
            
        end
        
        tempTimeCOM = tempTimeCOM(~cellfun('isempty', tempTimeCOM));
        tempTimeCOM2 = tempTimeCOM2(~cellfun('isempty', tempTimeCOM2));
        tempTimeCOMNorm = tempTimeCOMNorm(~cellfun('isempty', tempTimeCOMNorm));
        tempTimeGRF = tempTimeGRF(~cellfun('isempty', tempTimeGRF));
        tempTimeGRFNorm = tempTimeGRFNorm(~cellfun('isempty', tempTimeGRFNorm));
        
        tempTrajCOMHor = tempTrajCOMHor(~cellfun('isempty', tempTrajCOMHor));
        tempTrajCOMHorF1 = tempTrajCOMHorF1(~cellfun('isempty', tempTrajCOMHorF1));
        tempTrajCOMHorF2 = tempTrajCOMHorF2(~cellfun('isempty', tempTrajCOMHorF2));
        tempTrajCOMMedLat = tempTrajCOMMedLat(~cellfun('isempty', tempTrajCOMMedLat));
        tempTrajCOMVert = tempTrajCOMVert(~cellfun('isempty', tempTrajCOMVert));
        
        tempTrajGRFHor = tempTrajGRFHor(~cellfun('isempty', tempTrajGRFHor));
        tempTrajGRFVert = tempTrajGRFVert(~cellfun('isempty', tempTrajGRFVert));
        
        tempTrajCOMHor2(tempSkipInd,:) = [];
        tempTrajCOMMedLat2(tempSkipInd,:) = [];
        tempTrajCOMVert2(tempSkipInd,:) = [];
        tempTrajCOMHor2F1(tempSkipInd,:) = [];
        tempTrajCOMHor2F2(tempSkipInd,:) = [];
        
        tempTrajGRFHor1(tempSkipInd,:) = [];
        tempTrajGRFHor2(tempSkipInd,:) = [];
        tempTrajGRFVert1(tempSkipInd,:) = [];
        tempTrajGRFVert2(tempSkipInd,:) = [];
        
        % Calculate Average Gait Cycle Duration (GRF and COM Compare)
        tempTimeDurGRF = mean(cellfun(@(x) x(end), tempTimeGRF), 'omitnan');
        tempTimeDurCOM = mean(cellfun(@(x) x(end), tempTimeCOM), 'omitnan');
        
        % Store Individual Gait Cycles for COM
        mocap.gait.comTimeC = tempTimeCOM;
        mocap.gait.comTimeC2 = tempTimeCOM2;
        mocap.gait.comPosHorC = tempTrajCOMHor;
        mocap.gait.comPosMedLatC = tempTrajCOMMedLat;
        mocap.gait.comPosVertC = tempTrajCOMVert;
        mocap.gait.comPosHor = tempTrajCOMHor2;
        mocap.gait.comPosMedLat = tempTrajCOMMedLat2;
        mocap.gait.comPosVert = tempTrajCOMVert2;
        mocap.gait.comGait = tempTime;
        
        % Calculate Average Fore/Aft and Vertical COM Trajectory
        mocap.fit.comPosHor = mean(tempTrajCOMHor2, 1, 'omitnan');
        mocap.fit.comPosMedLat = mean(tempTrajCOMMedLat2, 1, 'omitnan');
        mocap.fit.comPosVert = mean(tempTrajCOMVert2, 1, 'omitnan');
        
        mocap.fit.comPosHorF1 = mean(tempTrajCOMHor2F1, 1, 'omitnan');
        mocap.fit.comPosHorF2 = mean(tempTrajCOMHor2F2, 1, 'omitnan');
        
        % Calculate Standard Deviation of Fore/Aft and Vertical COM
        mocap.fit.comPosHorSD = std(tempTrajCOMHor2, 0, 1, 'omitnan');
        mocap.fit.comPosVertSD = std(tempTrajCOMVert2, 0, 1, 'omitnan');
        
        mocap.fit.comPosHorFSD = mocap.fit.comPosHorSD;
        mocap.fit.comPosVertFSD = mocap.fit.comPosVertSD;
        
        % Calculate Average Fore/Aft and Vertical GRF Profile
        mocap.fit.grfHor = [mean(tempTrajGRFHor1, 1, 'omitnan')',...
            mean(tempTrajGRFHor2, 1, 'omitnan')'];
        mocap.fit.grfVert = [mean(tempTrajGRFVert1, 1, 'omitnan')',...
            mean(tempTrajGRFVert2, 1, 'omitnan')'];
        
        % Calculate Average Fore/Aft and Vertical GRF Profile
        mocap.fit.grfHorSD = [std(tempTrajGRFHor1, 0, 1, 'omitnan')',...
            std(tempTrajGRFHor2, 0, 1, 'omitnan')'];
        mocap.fit.grfVertSD = [std(tempTrajGRFVert1, 0, 1, 'omitnan')',...
            std(tempTrajGRFVert2, 0, 1, 'omitnan')'];
        
        % Calculate Average Gait Cycle Time Array
        mocap.fit.time = tempTime*tempTimeDurCOM;
        
        % Store Average Gait Cycle Normalized Time Array (COM and GRF)
        mocap.data.timeCOMNorm = tempTime;
        mocap.data.timeGRFNorm = tempTime;
        
        if params{9}
            
            mocap.fit.comPosHor = [mocap.fit.comPosHor...
                (mocap.fit.comPosHor(end) + mocap.fit.comPosHor(2:end))...
                (2*mocap.fit.comPosHor(end) + mocap.fit.comPosHor(2:end))];
            
            mocap.fit.comPosVert = [mocap.fit.comPosVert...
                repmat(mocap.fit.comPosVert(2:end), 1, 2)];
            
            mocap.fit.comPosHorF1 = [mocap.fit.comPosHorF1...
                repmat(mocap.fit.comPosHorF1(2:end), 1, 2)];
            mocap.fit.comPosHorF2 = [mocap.fit.comPosHorF2...
                (mocap.fit.comPosHorF2(end) +...
                mocap.fit.comPosHorF2(2:end))...
                (2*mocap.fit.comPosHorF2(end) +...
                mocap.fit.comPosHorF2(2:end))];
            
            mocap.fit.comPosHorSD = [mocap.fit.comPosHorSD...
                repmat(mocap.fit.comPosHorSD(2:end), 1, 2)];
            mocap.fit.comPosHorFSD = mocap.fit.comPosHorFSD;
            
            mocap.fit.comPosVertSD = [mocap.fit.comPosVertSD...
                repmat(mocap.fit.comPosVertSD(2:end), 1, 2)];
            mocap.fit.comPosVertFSD = mocap.fit.comPosVertFSD;
            
            mocap.fit.grfHor = [mocap.fit.grfHor;...
                repmat(mocap.fit.grfHor(2:end,:), 2, 1)];
            
            mocap.fit.grfVert = [mocap.fit.grfVert;...
                repmat(mocap.fit.grfVert(2:end,:), 2, 1)];
            
            mocap.fit.time = [mocap.fit.time...
                (mocap.fit.time(end) + mocap.fit.time(2:end))...
                (2*mocap.fit.time(end) + mocap.fit.time(2:end))];
            
            mocap.data.timeCOMNorm = [mocap.data.timeCOMNorm...
                (mocap.data.timeCOMNorm(end) +...
                mocap.data.timeCOMNorm(2:end))...
                (2*mocap.data.timeCOMNorm(end) +...
                mocap.data.timeCOMNorm(2:end))];
            mocap.data.timeGRFNorm = mocap.data.timeCOMNorm;
            
            mocap.data.timeDS1 = repmat(mocap.data.timeDS1, 1, 3);
            mocap.data.timeSS1 = repmat(mocap.data.timeSS1, 1, 3);
            mocap.data.timeDS2 = repmat(mocap.data.timeDS2, 1, 3);
            mocap.data.timeSS2 = repmat(mocap.data.timeSS2, 1, 3);
            
            mocap.data.timeFull = repmat(mocap.data.timeFull, 3, 1);

            mocap.data.timeGait = 3*mocap.data.timeGait;
            
        end
        
    else
        
        % Grab Time Array for Entire Data, Calculate Norm
            % COM
        tempTimeCOM = mocap.data.timeCOM(mocap.inds.TD1(...
            mocap.inds.start)/2:mocap.inds.TD1(mocap.inds.end)/2) -...
            mocap.data.timeCOM(mocap.inds.TD1(mocap.inds.start)/2);
        
        mocap.data.timeCOMNorm = tempTimeCOM/...
            (mocap.data.timeCOM(mocap.inds.TD1(mocap.inds.end)/2) -...
            mocap.data.timeCOM(mocap.inds.TD1(mocap.inds.start)/2));
        
            % GRF     
        tempTimeGRF = mocap.data.timeGRF(mocap.inds.TD1(...
            mocap.inds.start):mocap.inds.TD1(mocap.inds.end)) -...
            mocap.data.timeGRF(mocap.inds.TD1(mocap.inds.start));
        
        mocap.data.timeGRFNorm = tempTimeGRF/...
            (mocap.data.timeGRF(mocap.inds.TD1(mocap.inds.end)) -...
            mocap.data.timeGRF(mocap.inds.TD1(mocap.inds.start)))';
        
        % Grab Fore/Aft and Vertical COM Data along with Time Array        
        mocap.fit.comPosHor = ((mocap.data.COM(mocap.inds.TD1(...
            mocap.inds.start)/2:mocap.inds.TD1(mocap.inds.end)/2,1) -...
            mocap.data.COM(mocap.inds.TD1(mocap.inds.start)/2,1))/1000 +...
            mocap.subj.treadmill*tempTimeCOM)';
        
        mocap.fit.comPosHorF1 = ((mocap.data.COM(mocap.inds.TD1(...
            mocap.inds.start)/2:mocap.inds.TD1(mocap.inds.end)/2,1) -...
            mocap.data.COM(mocap.inds.TD1(mocap.inds.start)/2,1))/1000)';
        
        mocap.fit.comPosHorF2 = mocap.subj.treadmill*tempTimeCOM';

        mocap.fit.comPosVert = (mocap.data.COM(mocap.inds.TD1(...
            mocap.inds.start)/2:mocap.inds.TD1(mocap.inds.end)/2,2)/1000)';

        mocap.fit.time = (mocap.data.timeCOM(mocap.inds.TD1(...
            mocap.inds.start)/2:mocap.inds.TD1(mocap.inds.end)/2) -...
            mocap.data.timeCOM(mocap.inds.TD1(mocap.inds.start)/2))';
        
        % Grab Fore/Aft and Vertical GRF Profile
        mocap.fit.grfHor = mocap.data.grf.hor(...
            mocap.inds.TD1(mocap.inds.start):...
            mocap.inds.TD1(mocap.inds.end),:);
        mocap.fit.grfVert = mocap.data.grf.vert(...
            mocap.inds.TD1(mocap.inds.start):...
            mocap.inds.TD1(mocap.inds.end),:); 
    
    end

    %% %%%%%%%%%%%%%%%%%%%% COM Trajectory Polyfit %%%%%%%%%%%%%%%%%%%%% %%

    % Calculate polynomial coefficients for vertical and fore/aft COM
    mocap.fit.coefPosHor = polyfit(mocap.fit.time,...
        mocap.fit.comPosHor, mocap.fit.polyPosHor);
    mocap.fit.coefPosVert = polyfit(mocap.fit.time,...
        mocap.fit.comPosVert, mocap.fit.polyPosVert);
    
    mocap.fit.coefPosHorF1 = coeffvalues(fit(mocap.fit.time',...
        mocap.fit.comPosHorF1', 'fourier8'));
    mocap.fit.coefPosHorF2 = polyfit(mocap.fit.time,...
        mocap.fit.comPosHorF2, 1);
    mocap.fit.coefPosHorF = [mocap.fit.coefPosHorF2...
        mocap.fit.coefPosHorF1];
    mocap.fit.coefPosVertF = coeffvalues(fit(mocap.fit.time',...
        mocap.fit.comPosVert', 'fourier8'));
    
    % Polynomial coefficients for vertical and fore/aft COM velocity
    mocap.fit.coefVelHor = polyder(mocap.fit.coefPosHor);
    mocap.fit.coefVelVert = polyder(mocap.fit.coefPosVert);

    % Generate polyfit vertical and fore/aft COM trajectory
    mocap.fit.evalPosHor = polyval(mocap.fit.coefPosHor, mocap.fit.time);
    mocap.fit.evalPosVert = polyval(mocap.fit.coefPosVert, mocap.fit.time);
    
    mocap.fit.evalPosHorF = fitFourier(mocap.fit.coefPosHorF,...
        mocap.fit.time, 'Hor');
    mocap.fit.evalPosVertF = fitFourier(mocap.fit.coefPosVertF,...
        mocap.fit.time, 'Vert');
    
    mocap.fit.evalVelHor = polyval(mocap.fit.coefVelHor,mocap.fit.time);
    mocap.fit.evalVelVert = polyval(mocap.fit.coefVelVert,mocap.fit.time);
    
    mocap.fit.evalVelHorF = fitFourierDer(mocap.fit.coefPosHorF,...
        mocap.fit.time, 'Hor');
    mocap.fit.evalVelVertF = fitFourierDer(mocap.fit.coefPosVertF,...
        mocap.fit.time, 'Vert');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    if params{10}

%         figure
%         plot(mocap.fit.time, mocap.fit.comPosHor,'k')
%         hold on
%         plot(mocap.fit.time, mocap.fit.evalPosHor,'r--')
%         xlabel('Time [s]')
%         ylabel('Fore/Aft COM Trajectory [m]')
%         legend('MoCap COM','Polyfit COM')
        
%         figure
%         plot(mocap.fit.time, mocap.fit.comPosVert,'k')
%         hold on
%         plot(mocap.fit.time, mocap.fit.evalPosVert,'r--')
%         xlabel('Time [s]')
%         ylabel('Vertical COM Trajectory [m]')
%         legend('MoCap COM','Polyfit COM')
        
        figure
        plot(mocap.fit.time, mocap.fit.comPosHor,'k')
        hold on
        plot(mocap.fit.time, mocap.fit.evalPosHorF,'r--')
        xlabel('Time [s]')
        ylabel('Fore/Aft COM Trajectory [m]')
        legend('MoCap COM','Harmonic Fit COM')
        
        figure
        plot(mocap.fit.time, mocap.fit.comPosVert,'k')
        hold on
        plot(mocap.fit.time, mocap.fit.evalPosVertF,'r--')
        xlabel('Time [s]')
        ylabel('Vertical COM Trajectory [m]')
        legend('MoCap COM','Harmonic Fit COM')

        figure
        hold on
        plot(mocap.data.timeGRFNorm, mocap.fit.grfHor(:,2),...
            'g-.', 'Linewidth', 2)
        plot(mocap.data.timeGRFNorm, mocap.fit.grfVert(:,2),...
            'g-', 'Linewidth', 2)
        plot(mocap.data.timeGRFNorm, mocap.fit.grfHor(:,1),...
            'm-.', 'Linewidth', 2)
        plot(mocap.data.timeGRFNorm, mocap.fit.grfVert(:,1),...
            'm-', 'Linewidth', 2)
        plot(mocap.data.timeGRFNorm,...
            mocap.subj.bodyWeight*...
            ones(length(mocap.data.timeGRFNorm),1), 'k--')
        axis([0 mocap.data.timeGRFNorm(end)...
            -mocap.subj.bodyWeight*0.5 mocap.subj.bodyWeight*2])
        xlabel('Percentage of Gait Cycle')
        ylabel('Force [N]')
        legend('Leg 1 Fore/Aft GRF', 'Leg 1 Vertical GRF',...
            'Leg 2 Fore/Aft GRF', 'Leg 2 Vertical GRF',...
            'Bodyweight')

    end
    
    mocap.inds.range = mocap.inds.range - tempSkipCount;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Data Save %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    if params{12}

        clear temp* boolPlot boolSave markers force 

        save(mocap.subj.fileName)

    end

    %% %%%%%%%%%%%%%%%%%%%% User-Defined Functions %%%%%%%%%%%%%%%%%%%%% %%

    function cent = centroid(data,indices)

        size = length(indices);
        temp = zeros(length(data), 1);

        for ct = 1:size

            temp = temp + data(:, indices(ct));

        end

        cent = temp/size;

    end

end