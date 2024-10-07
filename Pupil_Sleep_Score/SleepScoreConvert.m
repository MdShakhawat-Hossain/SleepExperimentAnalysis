clc; clear all; close all

% character list of all training data files
procDataFileStruct = dir('*_TrainingData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    trainingDataSetID = [procDataFileID(1:end-16) 'TrainingData.mat'];
    load(trainingDataSetID)

    disp(['Converting to sleep score:  ' trainingDataSetID '...' ]); disp(' ')

    SleepScore.Label =   trainingTable.behavState;
    Score = zeros(length(SleepScore.Label),1);
    Score(SleepScore.Label == "Not Sleep") = 1;
    Score(SleepScore.Label == "NREM Sleep") = 2;
    Score(SleepScore.Label == "REM Sleep") = 3;
    Score = categorical(Score);

    SleepScore.Score = Score;
    SleepScoreID = [procDataFileID(1:end-16) 'SleepScore.mat'];
    save(SleepScoreID,'SleepScore')

    clear SleepScore Score trainingTable
end