
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Create ProcData structure using threshholds for the observed data
%            2) Extract the average pixel reflectance changes within those ROIs and save to RawData/ProcData files
%            3) Use spectral analysis of the reflectance dataper to pull out the animal's heart rate
%            4) Remove pixel-drift trends from each file, if necessary
%________________________________________________________________________________________________________________________   

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
% zap;
function StageTwoProcessing_FP_GRABNE()
clc; clear all; close all
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
[animalID,~,~] = GetFileInfo_FP(rawDataFileIDs(1,:));
% Identify whether this analysis is for bilateral or single hemisphere imaging
curDir = cd;
dirBreaks = strfind(curDir,'\');
curFolder = curDir(dirBreaks(end) + 1:end);
%% BLOCK PURPOSE: [1] Fiber Photometry data
disp('Analyzing Block [1] Processing TDT fiber photometry data.'); disp(' ')
% generate locations for fiberData
filelocation = pwd;
% calculate baseline data for fiber zscore
%%
% FiberBaseline = TDTFiberPhotometry_bilateralGRABNE_Fiber_Baseline(filelocation,rawDataFiles);
%%
% Final Fiber Analysis
% TDTFiberPhotometry_bilateralGRABNE_BaselineUpdate(filelocation,rawDataFiles,FiberBaseline);

%% run the fiber data analysis
codeVersion = input('Input stimulation type (z or perc or raw): ','s'); disp(' '); % select which analysis to run

    if strcmp(codeVersion,'z') == true
        % Original version
        TDTFiberPhotometry_bilateralGRABNE(filelocation,rawDataFiles);
    elseif strcmp(codeVersion,'perc') == true
        % percentage changes
        time_Trim_Start = input('How many seconds of data to be removed from start: '); disp(' ')
        time_Trim_End = input('How many seconds of data to be removed from end: '); disp(' ')
        TDTFiberPhotometry_bilateralGRABNE_Percent(filelocation,rawDataFiles,time_Trim_Start,time_Trim_End);
    elseif strcmp(codeVersion,'raw') == true
        % Raw changes
        TDTFiberPhotometry_bilateralGRABNE_RawChange(filelocation,rawDataFiles);
    elseif strcmp(codeVersion,'none') == true
        time_Trim_Start = input('How many seconds of data to be removed from start: '); disp(' ')        
    else
        error('wrong input'); 
        dbquit; % quit the debug mode
    end
%% Block Purpose: [2] Process the Deeplabcut tracked pupil diameter from CSV files
disp('Analyzing Block [2] processing DLC tracked pupil diameter.'); disp(' ')
DLCPupilTrack_Processing_8point(animalID);
%% BLOCK PURPOSE: [3] Correct the offset between the MScan and LabVIEW acquisiton.
disp('Analyzing Block [3] Correcting LabVIEW time offset.'); disp(' ')
fiberDataFileStruct = dir('*_FiberData.mat');
fiberDataFiles = {fiberDataFileStruct.name}';
fiberDataFileID = char(fiberDataFiles);
trimTime = time_Trim_Start; %seconds
TemplateMatchFiberData_FP_GRABNE(fiberDataFileID,rawDataFileIDs,trimTime)
%% BLOCK PURPOSE: [4] Process the RawData structure -> Create Threshold data structure and ProcData structure.
disp('Analyzing Block [4] Creating ProcData files and processing analog data.'); disp(' ')
ProcessRawDataFiles_FP_GRABNE(rawDataFileIDs)
%% fin.
disp('Fiber Photometry Stage Two Processing - Complete.'); disp(' ')
