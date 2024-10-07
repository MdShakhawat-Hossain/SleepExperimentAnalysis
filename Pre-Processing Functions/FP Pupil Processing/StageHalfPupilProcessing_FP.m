    
% function StageHalfPupilProcessing_FP()

clc; clear all; close all
disp('Analyzing Block [0] reading the raw files to process the pupil data.'); disp(' ')
% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
FileID = rawDataFileIDs(1,:);
[animalID,~,~] = GetFileInfo_FP_Shak(FileID);
%% BLOCK PURPOSE: [1] Pupil Data
disp('Analyzing Block [1] loading pupil camera data and process pupil changes.'); disp(' ')
% run pupil tracker
% RunPupilTracker_JNeurosci2022_FP(rawDataFileIDs)
RunPupilTracker_FP(rawDataFileIDs)
% patch any missing frames.
disp('Patch any missing frames for pupil measurements'); disp(' ')
PatchPupilArea_FP(rawDataFileIDs) 

 % update the centroid and major minor axises for missing frames
disp('Update the centroid and major minor axises for missing frames'); disp(' ')
AccessoryPupilUpdate_FP(rawDataFileIDs)

% convert area to diameter. Convert pixels to mm
disp('Convert area to diameter. Convert pixels to mm'); disp(' ')
ConvertPupilAreaToDiameter_FP(rawDataFileIDs) 

