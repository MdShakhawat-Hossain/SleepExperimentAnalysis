function StageThreeProcessing_FP_EventRest()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 1) Categorize behavioral (rest,whisk,stim) data using previously processed data structures, add 'flags'
%            2) Create a temporary RestData structure that contains periods of rest - use this for initial figures
%            3) Analyze neural data and create different spectrograms for each file's electrodes
%            4) Uses periods when animal is not being stimulated or moving to establish an initial baseline
%            5) Manually select awake files for a slightly different baseline not based on hard time vals
%            6) Use the best baseline to convert reflectance changes to total hemoglobin
%            7) Re-create the RestData structure now that we can deltaHbT
%            8) Create an EventData structure looking at the different data types after whisking or stimulation
%            9) Apply the resting baseline to each data type to create a percentage change
%            10) Use the time indeces of the resting baseline file to apply a percentage change to the spectrograms
%            11) Use the time indeces of the resting baseline file to create a reflectance pixel-based baseline
%            12) Generate a summary figure for all of the analyzed and processed data
%________________________________________________________________________________________________________________________
% addpath(genpath('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Documents\Research_Codes\FiberPhotometry\Data-Analysis-master'))
%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command windyow.
% zap;
% clear all; close all; clc
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
[animalID,~,~] = GetFileInfo_FP(procDataFileIDs(1,:));
% parameters used for various animal analysis
curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
stimulationType = input('Input stimulation type (single or pulse): ','s'); disp(' ')
dataTypes = {'TRITC','GCaMP7s','cortical_LH','cortical_RH','hippocampus','EMG'};
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
basefile = ([animalID '_RestingBaselines.mat']);
%% BLOCK PURPOSE: [1] Categorize data
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Analyzing file ' num2str(a) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
    CategorizeData_FP(procDataFileID,stimulationType)
end
%% BLOCK PURPOSE: [2] Create RestData data structure
disp('Analyzing Block [2] Create RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData_FP(procDataFileIDs,dataTypes);
%% BLOCK PURPOSE: [8] Create the EventData structure for CBV and neural data
disp('Analyzing Block [8] Create EventData struct for CBV and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_FP(procDataFileIDs,dataTypes);
%% BLOCK PURPOSE: [3] Analyze the spectrogram for each session.
disp('Analyzing Block [3] Analyzing the spectrogram for each file.'); disp(' ')
CreateTrialSpectrograms_FP(rawDataFileIDs,neuralDataTypes);
%% BLOCK PURPOSE: [4] Create Baselines data structure
disp('Analyzing Block [4] Create baselines structure for CBV and neural data.'); disp(' ')
baselineType = 'setDuration';
% calculate trial duration
load(procDataFileIDs(1,:))
trialDuration_sec = length(ProcData.data.EMG.emg)/ProcData.notes.dsFs;
targetMinutes = trialDuration_sec/60;
[RestingBaselines] = CalculateRestingBaselines_FP(animalID,targetMinutes,trialDuration_sec,RestData);
% Find spectrogram baselines for each day
[RestingBaselines] = CalculateSpectrogramBaselines_FP(animalID,neuralDataTypes,trialDuration_sec,RestingBaselines,baselineType);
% Normalize spectrogram by baseline
NormalizeSpectrograms_FP(neuralDataTypes,RestingBaselines);
%% BLOCK PURPOSE: [5] Manually select files for custom baseline calculation
disp('Analyzing Block [5] Manually select files for custom baseline calculation.'); disp(' ')
% [RestingBaselines] = CalculateManualRestingBaselinesTimeIndeces_FP();
[RestingBaselines] = CalculateManualRestingBaselinesTimeIndeces_FP_No_Redo(); % Don't need to redo  the baseline selection
%% BLOCK PURPOSE: [9] Normalize RestData and EventData structures by the resting baseline
updatedBaselineType = 'manualSelection';
% Character list of all ProcData files
restDataFileStruct = dir('*_RestData.mat');
restDataFiles = {restDataFileStruct.name}';
restDataFileIDs = char(restDataFiles);
load(restDataFileIDs)
% Character list of all ProcData files
eventDataFileStruct = dir('*_EventData.mat');
eventDataFiles = {eventDataFileStruct.name}';
eventDataFileIDs = char(eventDataFiles);
load(eventDataFileIDs)
% Character list of all ProcData files
baseDataFileStruct = dir('*_RestingBaselines.mat');
baseDataFiles = {baseDataFileStruct.name}';
baseDataFileIDs = char(baseDataFiles);
load(baseDataFileIDs)
disp('Analyzing Block [9] Normalizing RestData and EventData structures by the resting baseline.'); disp(' ')
[RestData] = NormBehavioralDataStruct_FP(RestData,RestingBaselines,updatedBaselineType);
save([animalID '_RestData.mat'],'RestData','-v7.3')
[EventData] = NormBehavioralDataStruct_FP(EventData,RestingBaselines,updatedBaselineType);
save([animalID '_EventData.mat'],'EventData','-v7.3')
%% BLOCK PURPOSE: [10] Analyze the spectrogram baseline for each session.
disp('Analyzing Block [10] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
% Find spectrogram baselines for each day
% [RestingBaselines] = CalculateSpectrogramBaselines_FP(animalID,neuralDataTypes,trialDuration_sec,RestingBaselines,updatedBaselineType);
% Normalize spectrogram by baseline
% NormalizeSpectrograms_FP(neuralDataTypes,RestingBaselines);
% Create a structure with all spectrograms for convenient analysis further downstream
% CreateAllSpecDataStruct_FP(animalID,neuralDataTypes)
%% BLOCK PURPOSE [11] Generate single trial figures
% disp('Analyzing Block [11] Generating single trial summary figures'); disp(' ')
% saveFigs = 'y';
% for bb = 1:size(procDataFileIDs,1)
%     procDataFileID = procDataFileIDs(bb,:);
%     [figHandle] = GenerateSingleFigures_FP(procDataFileID,saveFigs);
%     close(figHandle)
% end

disp('Stage Three Processing - Complete.'); disp(' ')
