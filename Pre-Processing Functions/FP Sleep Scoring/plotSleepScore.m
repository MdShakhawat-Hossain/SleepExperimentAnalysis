% function plotSleepScore()
clc;
clear all; close all;
% FileLocation = 'H:\Sleep_GCaMP7s_ChATCre\T281\CombinedImaging';
% cd(FileLocation);
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    [figHandle,ax1,ax2,ax3,ax5,ax6] = GenerateSingleFigures_Sleep_FP_Plot(procDataFileID,'y');
%    wait
    pause
end
