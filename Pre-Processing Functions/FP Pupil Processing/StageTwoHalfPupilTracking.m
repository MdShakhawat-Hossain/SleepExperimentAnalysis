
function StageTwoHalfPupilTracking

clc; clear all; close all
disp('Analyzing Block [0] reading the proc files to check the quality of pupil analysis.'); disp(' ')
% Character list of all RawData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
%% check the quality of pupil analysis
% check eye quality
for bb = 1:size(procDataFileIDs,1)
    disp(['Manually checking eye of file ' num2str(bb) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilVideoFrames_FP(procDataFileIDs(bb,:))
end
% verify blinks
for cc = 1:size(procDataFileIDs,1)
    disp(['Manually checking blinks of file ' num2str(cc) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilBlinks_FP(procDataFileIDs(cc,:))
end
% verify pupil area
for dd = 1:size(procDataFileIDs,1)
    disp(['Manually checking pupil area of file ' num2str(dd) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilDiameter_FP(procDataFileIDs(dd,:))
end


