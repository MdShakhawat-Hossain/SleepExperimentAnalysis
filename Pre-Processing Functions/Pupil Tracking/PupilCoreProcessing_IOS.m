function [] = PupilCoreProcessing_IOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 1) Track pupil diameter and detect periods of blinking
%            2) Patch any NaN values or droppped camera frames via interpolation
%            3) Manually check the first 5 and last 5 frames of each session to verify eye integrity/discharge
%            4) Manually check each blink for false positives
%            5) Manually check each file's pupil diameter
%            6) Extract resting pupil area and add it to RestData.mat structure
%            7) Determine baseline pupil area during rest and add it to RestingBaselines.mat
%            8) Extract whisking/stimulus triggered pupil area and add it to EventData.mat
%            9) Normalize RestData.mat and EventData.mat structures using resting baseline
%           10) Update pupil data in SleepData.mat
%________________________________________________________________________________________________________________________

%% load the script's necessary variables and data structures.
% clear the workspace variables and command windyow.
zap;
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
[animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
%% track pupil area and blink detetction
RunPupilTracker_IOS(procDataFileIDs)
%% patch pupil area
for aa = 1:size(procDataFileIDs,1)
    disp(['Patching pupil area of file ' num2str(aa) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    PatchPupilArea_IOS(procDataFileIDs(aa,:))
end
%% check eye quality
for bb = 1:size(procDataFileIDs,1)
    disp(['Manually checking eye of file ' num2str(bb) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilVideoFrames_IOS(procDataFileIDs(bb,:))
end
%% verify blinks
for cc = 1:size(procDataFileIDs,1)
    disp(['Manually checking blinks of file ' num2str(cc) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilBlinks_IOS(procDataFileIDs(cc,:))
end
%% verify pupil area
for dd = 1:size(procDataFileIDs,1)
    disp(['Manually checking pupil area of file ' num2str(dd) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilDiameter_IOS(procDataFileIDs(dd,:))
end
%% convert area to diameter
for ee = 1:size(procDataFileIDs,1)
    disp(['Converting pupil area to pupil diameter of file ' num2str(ee) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    ConvertPupilAreaToDiameter_IOS(procDataFileIDs(ee,:))
end
%% add pupil area to RestData.mat
dataTypes = {'pupilArea','diameter','mmArea','mmDiameter'};
ExtractPupilRestingData_IOS(procDataFileIDs,dataTypes);
%% add pupil baseline to Restingbaselines.mat
[RestingBaselines] = AddPupilRestingBaseline_IOS();
%% zScore pupil data
for ff = 1:size(procDataFileIDs,1)
    disp(['Z-scoring pupil data of file ' num2str(ff) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    procDataFileID = procDataFileIDs(ff,:);
    zScorePupilData_IOS(procDataFileID,RestingBaselines)
end
%% add pupil area to RestData.mat
dataTypes = {'pupilArea','diameter','mmArea','mmDiameter','zArea','zDiameter','LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
[RestData] = ExtractPupilRestingData_IOS(procDataFileIDs,dataTypes);
%% add pupil area to EventData.mat
[EventData] = ExtractPupilEventTriggeredData_IOS(procDataFileIDs);
%% normalize Rest/Event data structures
[RestData] = NormBehavioralDataStruct_IOS(RestData,RestingBaselines,'manualSelection');
save([animalID '_RestData.mat'],'RestData','-v7.3')
[EventData] = NormBehavioralDataStruct_IOS(EventData,RestingBaselines,'manualSelection');
save([animalID '_EventData.mat'],'EventData','-v7.3')
%% add pupil data to SleepData.mat
AddPupilSleepParameters_IOS(procDataFileIDs)
UpdatePupilSleepData_IOS(procDataFileIDs)
