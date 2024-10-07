function [animalID] = GetAnimalID()

startingDirectory = cd;
trainingDirectory = [startingDirectory];% '\Training Data\'];
cd(trainingDirectory)
% character list of all training files
trainingDataFileStruct = dir('*_ModelData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
trainingTableFileID = trainingDataFileIDs(1,:);
% pull animal ID
[animalID,~,~] = GetFileInfo_FP(trainingTableFileID);
end