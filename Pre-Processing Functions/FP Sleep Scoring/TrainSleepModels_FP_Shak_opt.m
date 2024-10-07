% function [animalID] = TrainSleepModels_FP_Shak_opt()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Train several machine learning techniques on manually scored sleep data, and evaluate each model's accuracy
%________________________________________________________________________________________________________________________

%% load in all the data to create a table of values
startingDirectory = cd;
trainingDirectory = [startingDirectory];% '\Training Data\'];
cd(trainingDirectory)
% character list of all training files
trainingDataFileStruct = dir('*_ModelData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
% Load each updated training set and concatenate the data into table
for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    dataLength = size(ModelData.TrainLabels,1);
    for lp = 1:1:dataLength
        AllFeatures{((bb-1)*dataLength)+lp } = ModelData.Features{lp}';
    end
    
    AllLabels(((bb-1)*dataLength)+1 : (bb*dataLength)) = ModelData.TrainLabels;
end

% pull animal ID
[animalID,~,~] = GetFileInfo_FP(trainingTableFileID);
% directory path for saving data
dirpath = [startingDirectory '\Figures\Sleep Models\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
%% prepare data
crossvaridation = cvpartition(AllLabels,'HoldOut',0.3,'Stratify',true);
test_idx = crossvaridation.test;
train_idx = crossvaridation.training;
trainData  = AllFeatures(train_idx);
trainLabel  = AllLabels(train_idx)';
testData  = AllFeatures(test_idx);
testLabel  = AllLabels(test_idx)';
%% Network architecture
 [bestIdx,BestfileName] = CNN_Opti(animalID,trainData,trainLabel,testData,testLabel);
%% train
savedStruct = load(BestfileName);  %load('H:\Sleep_GCaMP7s_ChATCre\T281\CombinedImaging\T281_0.11004.mat');
valError = savedStruct.valError
%% predict
YPredCNN = classify(savedStruct.trainedNet,testData, ...
    MiniBatchSize=64, ...
    SequencePaddingDirection="left");
%% calculate accuracy
CNN_confMat = figure;plotconfusion(testLabel,YPredCNN)
%% save model in desired location
CNN_MDL = savedStruct.trainedNet;
save([dirpath animalID '_FP_CNN_SleepScoringModel.mat'],'CNN_MDL')
% save model and figure
savefig(CNN_confMat,[dirpath animalID '_FP_CNN_ConfusionMatrix']);
close(CNN_confMat)
cd(startingDirectory)
%% save confusion matrix results
% save([dirpath animalID '_ConfusionData.mat'],'ConfusionData')
% end
