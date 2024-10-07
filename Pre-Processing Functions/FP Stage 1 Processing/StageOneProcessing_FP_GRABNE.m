%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%o
%   Purpose: Data acquired during trials must be in a form that Matlab can work with easily. This code converts the
%            various forms of data listed below into MATLAB structures that can be easily manipulated.
%
%            .bin - Cameras
%            .tdms - Digital and Analog Data
%            .tdms_index - Index for the LabVIEW data in the .tdms file
%________________________________________________________________________________________________________________________
% function StageOneProcessing_FP_GRABNE
%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window
zap;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Asks the user to load all files with a '_WhiskerCam.bin' extension
fileNames = uigetfile('*_WhiskerCam.bin','MultiSelect','on');   % CTL-A to select all files
% find the filepath to correct for naming error
filepath = pwd;
AnimalIndX = strfind(filepath(1,:),'\');
AnimalIDX = filepath(AnimalIndX(end-1)+1:AnimalIndX(end)-1);
AnimalIndx = strfind(AnimalIDX(1,:),'-');
AnimalID = AnimalIDX(1:AnimalIndx(1)-1);
%% BLOCK PURPOSE: [1] Preparing to create RawData files.
disp('Analyzing Block [1] Preparing to create RawData file(s).'); disp(' ')
% Load in each file one at a time, looping through the list
for a = 1:length(fileNames)
    disp(['Analyzing WhiskerCam file (' num2str(a) ' of ' num2str(length(fileNames)) ')']); disp(' ')
    % Adapt to list or single file. The purpose of this is control the way uigetfile handles an instance of a
    % single file input (character string) vs. multiple files, which it puts in cells
    if iscell(fileNames) == true
        indFile = fileNames{a};
    else
        indFile = fileNames;
    end
    % Pull out the file ID for the file - this is the numerical string after the animal name/hemisphere
    [~,~,fileID] = GetFileInfo_FP(indFile);
    % Determine if a RawData file has already been created for this file. If it has, skip it
%     fileExist = ls(['*' fileID '_RawData.mat']);
%     if isempty(fileExist)
        %% BLOCK PURPOSE: [2] Import .tdms data (All channels).
        disp('Analyzing Block [2] Importing .tdms data from all channels.'); disp(' ')
        trialData = ReadInTDMSWhiskerTrials_FP([fileID '.tdms']);
        % Left, Right, and hippocampal electrodes
        dataRow = strcmp(trialData.data.names,'cortLH');
        cortical_LH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'cortRH');
        cortical_RH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        dataRow = strcmp(trialData.data.names,'hipp');
        hippocampus = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        
        % Force sensor and EMG
        dataRow = strcmp(trialData.data.names,'forceSensor');
        forceSensor = trialData.data.vals(dataRow,:);
        dataRow = strcmp(trialData.data.names,'EMG');
        EMG = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
        % Left, Right, Auditory solenoids. Combine the arrays together.
        dataRow = strcmp(trialData.data.names,'LPadSol');
        LPadSol = gt(trialData.data.vals(dataRow,:),0.5)*1; % ID amplitude is 1
        dataRow = strcmp(trialData.data.names,'RPadSol');
        RPadSol = gt(trialData.data.vals(dataRow,:),0.5)*2; % ID amplitude is 2
        dataRow = strcmp(trialData.data.names,'AudSol');
        AudSol = gt(trialData.data.vals(dataRow,:),0.5)*3;  % ID amplitude is 3
        stimulations = LPadSol + RPadSol + AudSol;
        if sum(strcmp(trialData.data.names,'OptoStim'))
            dataRow = strcmp(trialData.data.names,'OptoStim');
            OptoStim = gt(trialData.data.vals(dataRow,:),0.5)*4;  % ID amplitude is 4
            stimulations = LPadSol + RPadSol + AudSol + OptoStim;
        end
        %% BLOCK PURPOSE: [3] Save the notes and data.
        disp('Analyzing Block [3] Evaluating data to save to RawData file.'); disp(' ')
        % notes - all variables are descriptive
        RawFileId = [AnimalID '_' fileID '_RawData.mat'];
        try % try to see if a rawData file already exist
            load(RawFileId);
            RawData.notes.experimenter = trialData.experimenter;
            RawData.notes.hemisphere = trialData.hemisphere;
            RawData.notes.solenoidPSI = str2double(trialData.solenoidPSI);
            RawData.notes.isofluraneTime = str2double(trialData.isofluraneTime);
            RawData.notes.sessionID = trialData.sessionID;
            RawData.notes.amplifierGain = str2double(trialData.amplifierGain);        
            RawData.notes.whiskCamSamplingRate = str2double(trialData.whiskCamSamplingRate);
            RawData.notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
            RawData.notes.trialDuration_long = str2double(trialData.trialDuration_sec);
            RawData.notes.whiskCamPixelHeight = str2double(trialData.whiskCamPixelHeight);
            RawData.notes.whiskCamPixelWidth = str2double(trialData.whiskCamPixelWidth);
            RawData.notes.pupilCamPixelHeight = str2double(trialData.pupilCamPixelHeight);
            RawData.notes.pupilCamPixelWidth = str2double(trialData.pupilCamPixelWidth);
            RawData.notes.pupilCamSamplingRate = str2double(trialData.pupilCamSamplingRate);
            RawData.notes.droppedPupilCamFrameIndex = str2double(trialData.droppedPupilCamFrameIndex);
            RawData.notes.droppedwhiskCamFrameIndex = str2double(trialData.droppedwhiskCamFrameIndex);
            % stimulus parameters
            RawData.notes.Interstim = str2double(trialData.Interstim);
            RawData.notes.TrialOffset= str2double(trialData.TrialOffset);
            RawData.notes.SolenoidDuration = str2double(trialData.SolenoidDuration);
            RawData.notes.SolenoidDutycyle = str2double(trialData.SolenoidDutycyle);
            RawData.notes.SolenoidFreq = str2double(trialData.SolenoidFreq);
        catch % raw file doesn't exist proceed with regular task
            RawData.notes.experimenter = trialData.experimenter;
            RawData.notes.animalID = AnimalID;%trialData.animalID;
            RawData.notes.hemisphere = trialData.hemisphere;
            RawData.notes.solenoidPSI = str2double(trialData.solenoidPSI);
            RawData.notes.isofluraneTime = str2double(trialData.isofluraneTime);
            RawData.notes.sessionID = trialData.sessionID;
            RawData.notes.amplifierGain = str2double(trialData.amplifierGain);        
            RawData.notes.whiskCamSamplingRate = str2double(trialData.whiskCamSamplingRate);
            RawData.notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
            RawData.notes.trialDuration_long = str2double(trialData.trialDuration_sec);
            RawData.notes.whiskCamPixelHeight = str2double(trialData.whiskCamPixelHeight);
            RawData.notes.whiskCamPixelWidth = str2double(trialData.whiskCamPixelWidth);
            RawData.notes.pupilCamPixelHeight = str2double(trialData.pupilCamPixelHeight);
            RawData.notes.pupilCamPixelWidth = str2double(trialData.pupilCamPixelWidth);
            RawData.notes.pupilCamSamplingRate = str2double(trialData.pupilCamSamplingRate);
            RawData.notes.droppedPupilCamFrameIndex = str2double(trialData.droppedPupilCamFrameIndex);
            RawData.notes.droppedwhiskCamFrameIndex = str2double(trialData.droppedwhiskCamFrameIndex);
            % stimulus parameters
            RawData.notes.Interstim = str2double(trialData.Interstim);
            RawData.notes.TrialOffset= str2double(trialData.TrialOffset);
            RawData.notes.SolenoidDuration = str2double(trialData.SolenoidDuration);
            RawData.notes.SolenoidDutycyle = str2double(trialData.SolenoidDutycyle);
            RawData.notes.SolenoidFreq = str2double(trialData.SolenoidFreq);
        end
        % Data
        RawData.data.cortical_LH_long = cortical_LH;
        RawData.data.cortical_RH_long = cortical_RH;
        RawData.data.hippocampus_long = hippocampus;
        RawData.data.forceSensor_long = forceSensor;
        RawData.data.EMG_long = EMG;
        RawData.data.stimulations_long = stimulations;
        %% BLOCK PURPOSE: [4] Start Whisker tracker.
        disp('Analyzing Block [4] Starting whisker tracking.'); disp(' ')
        if isfield(RawData.data,'whiskerAngle_long') == false || isfield(RawData.data,'whiskerAngle_long') == true
            [whiskerAngle] = WhiskerTrackerParallel_FP(fileID,RawData.notes);
            indsA = isnan(whiskerAngle); 
            indsB = indsA == 1; 
            whiskerAngle(indsB) = [];
            RawData.data.whiskerAngle_long = whiskerAngle;
        else
            disp('Whisker tracking already exist.Continuing...'); disp(' ')
        end
        disp(['File Created. Saving RawData File ' num2str(a) '...']); disp(' ')
        save([AnimalID '_' fileID '_RawData'],'RawData','-v7.3')
        clear RawData
end
disp('Fiber Photometry One Processing - Complete.'); disp(' ')
