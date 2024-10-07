function [animalID] = TrainSleepModels_FP_GRABNE()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Train several machine learning techniques on manually scored sleep data, and evaluate each model's accuracy
%________________________________________________________________________________________________________________________
clear all
close all
clc
%% load in all the data to create a table of values
startingDirectory = cd;
trainingDirectory = [startingDirectory];% '\Training Data\'];
cd(trainingDirectory)
% character list of all training files
trainingDataFileStruct = dir('*_TrainingData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
% Load each updated training set and concatenate the data into table
for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    if bb == 1
        joinedTable = trainingTable;
    end
    joinedTable = vertcat(joinedTable,trainingTable);
    dataLength = size(trainingTable,1);
end
%% separate the data for test and train
% prepare data
AllLabels = joinedTable.behavState;
% AllLabels = trainingTable.behavState;
crossvaridation = cvpartition(AllLabels,'HoldOut',0.2,'Stratify',true);
test_idx = crossvaridation.test;
train_idx = crossvaridation.training;

XTrain  = joinedTable(train_idx,1:end-1);
YTrain  = joinedTable(train_idx,end);
XTest  = joinedTable(test_idx,1:end-1);
YTest  = joinedTable(test_idx,end);
%% pull animal ID
[animalID,~,~] = GetFileInfo_FP(trainingTableFileID);
% directory path for saving data
dirpath = [startingDirectory '\Figures\Sleep Models\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
%% Train Support Vector Machine (SVM) classifier
t = templateSVM('Standardize',true,'KernelFunction','gaussian');
disp('Training Support Vector Machine...'); disp(' ')
SVM_MDL = fitcecoc(XTrain,YTrain,'Learners',t,'FitPosterior',true,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Verbose',2);
% save model in desired location
save([dirpath animalID '_FP_SVM_SleepScoringModel.mat'],'SVM_MDL')
% determine k-fold loss of the model
disp('Cross-validating the support vector machine classifier...'); disp(' ')
CV_SVM_MDL = crossval(SVM_MDL,'kfold',10);
loss = kfoldLoss(CV_SVM_MDL);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
% use the model to generate a set of scores for the Test set of data
XTrainLabels = predict(SVM_MDL,XTrain);
YTrains = table2array(YTrain);
XTestLabels = predict(SVM_MDL,XTest);
YTests = table2array(YTest);
%preallocate
XTrainLabelK = zeros(3,length(XTrainLabels));
YTrainK = zeros(3,length( YTrains));
XTestLabelK = zeros(3,length(XTestLabels));
YTestK = zeros(3,length( YTests));
for LXT = 1:length(XTrainLabels)
    if XTrainLabels(LXT) == "Not Sleep"
        XTrainLabelK(1,LXT) = 1;
    elseif XTrainLabels(LXT) == "NREM Sleep"
        XTrainLabelK(2,LXT) = 1;
    elseif XTrainLabels(LXT) == "REM Sleep"
        XTrainLabelK(3,LXT) = 1;
    end
    if YTrains(LXT) == "Not Sleep"
        YTrainK(1,LXT) = 1;
    elseif YTrains(LXT) == "NREM Sleep"
        YTrainK(2,LXT) = 1;
    elseif YTrains(LXT) == "REM Sleep"
        YTrainK(3,LXT) = 1;
    end
end

for LXT = 1:length(XTestLabels)
    if XTestLabels(LXT) == "Not Sleep"
        XTestLabelK(1,LXT) = 1;
    elseif XTestLabels(LXT) == "NREM Sleep"
        XTestLabelK(2,LXT) = 1;
    elseif XTestLabels(LXT) == "REM Sleep"
        XTestLabelK(3,LXT) = 1;
    end
    if YTests(LXT) == "Not Sleep"
        YTestK(1,LXT) = 1;
    elseif YTests(LXT) == "NREM Sleep"
        YTestK(2,LXT) = 1;
    elseif YTests(LXT) == "REM Sleep"
        YTestK(3,LXT) = 1;
    end
end

SVM_confMat = figure; 
subplot(1,2,1)
plotconfusion(YTrainK,XTrainLabelK);
title('Train Confusion')
savefig(SVM_confMat,[dirpath animalID '_FP_SVM_Train_ConfusionMatrix']);
subplot(1,2,2)
plotconfusion(YTestK,XTestLabelK);
title('Test Confusion')
savefig(SVM_confMat,[dirpath animalID '_FP_SVM_Test_ConfusionMatrix']);
close(SVM_confMat)
ConfusionData.SVM.trainYlabels = YTrain.behavState;
ConfusionData.SVM.trainXlabels = XTrainLabels;
ConfusionData.SVM.testYlabels = YTest.behavState;
ConfusionData.SVM.testXlabels = XTestLabels;
%% Ensemble classification - AdaBoostM2, Subspace, Bag, LPBoost,RUSBoost, TotalBoost
disp('Training Ensemble Classifier...'); disp(' ')
t = templateTree('Reproducible',true);
EC_MDL = fitcensemble(XTrain,YTrain,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% save model in desired location
save([dirpath animalID '_FP_EC_SleepScoringModel.mat'],'EC_MDL')
% determine k-fold loss of the model
disp('Cross-validating (3-fold) the ensemble classifier...'); disp(' ')
CV_EC_MDL = crossval(EC_MDL,'kfold',3);
loss = kfoldLoss(CV_EC_MDL);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
% use the model to generate a set of scores for the Test set of data
XTrainLabels = predict(EC_MDL,XTrain);
YTrains = table2array(YTrain);
XTestLabels = predict(EC_MDL,XTest);
YTests = table2array(YTest);
%preallocate
XTrainLabelK = zeros(3,length(XTrainLabels));
YTrainK = zeros(3,length( YTrains));
XTestLabelK = zeros(3,length(XTestLabels));
YTestK = zeros(3,length( YTests));

for LXT = 1:length(XTrainLabels)
    if XTrainLabels(LXT) == "Not Sleep"
        XTrainLabelK(1,LXT) = 1;
    elseif XTrainLabels(LXT) == "NREM Sleep"
        XTrainLabelK(2,LXT) = 1;
    elseif XTrainLabels(LXT) == "REM Sleep"
        XTrainLabelK(3,LXT) = 1;
    end
    if YTrains(LXT) == "Not Sleep"
        YTrainK(1,LXT) = 1;
    elseif YTrains(LXT) == "NREM Sleep"
        YTrainK(2,LXT) = 1;
    elseif YTrains(LXT) == "REM Sleep"
        YTrainK(3,LXT) = 1;
    end
end

for LXT = 1:length(XTestLabels)
    if XTestLabels(LXT) == "Not Sleep"
        XTestLabelK(1,LXT) = 1;
    elseif XTestLabels(LXT) == "NREM Sleep"
        XTestLabelK(2,LXT) = 1;
    elseif XTestLabels(LXT) == "REM Sleep"
        XTestLabelK(3,LXT) = 1;
    end
    if YTests(LXT) == "Not Sleep"
        YTestK(1,LXT) = 1;
    elseif YTests(LXT) == "NREM Sleep"
        YTestK(2,LXT) = 1;
    elseif YTests(LXT) == "REM Sleep"
        YTestK(3,LXT) = 1;
    end
end

EC_confMat = figure; 
subplot(1,2,1)
plotconfusion(YTrainK,XTrainLabelK);
title('Train Confusion')
savefig(EC_confMat,[dirpath animalID '_FP_EC_Train_ConfusionMatrix']);
subplot(1,2,2)
plotconfusion(YTestK,XTestLabelK);
title('Test Confusion')
savefig(EC_confMat,[dirpath animalID '_FP_EC_Test_ConfusionMatrix']);
close(EC_confMat)
% save labels for later confusion matrix
ConfusionData.EC.trainYlabels= YTrain.behavState;
ConfusionData.EC.trainXlabels = XTrainLabels;
ConfusionData.EC.testYlabels = YTest.behavState;
ConfusionData.EC.testXlabels = XTestLabels;
%% Random forest
disp('Training Random Forest Classifier...'); disp(' ')
numTrees = 128;
RF_MDL = TreeBagger(numTrees,XTrain,YTrain,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% save model in desired location
save([dirpath animalID '_FP_RF_SleepScoringModel.mat'],'RF_MDL')
% determine the misclassification probability (for classification trees) for out-of-bag observations in the training data
RF_OOBerror = oobError(RF_MDL,'Mode','Ensemble');
disp(['Random Forest out-of-bag error: ' num2str(RF_OOBerror*100) '%']); disp(' ')
% use the model to generate a set of scores for the Test set of data
XTrainLabels = predict(RF_MDL,XTrain);
YTrains = table2array(YTrain);
XTestLabels = predict(RF_MDL,XTest);
YTests = table2array(YTest);
%preallocate
XTrainLabelK = zeros(3,length(XTrainLabels));
YTrainK = zeros(3,length( YTrains));
XTestLabelK = zeros(3,length(XTestLabels));
YTestK = zeros(3,length( YTests));
for LXT = 1:length(XTrainLabels)
    if XTrainLabels(LXT) == "Not Sleep"
        XTrainLabelK(1,LXT) = 1;
    elseif XTrainLabels(LXT) == "NREM Sleep"
        XTrainLabelK(2,LXT) = 1;
    elseif XTrainLabels(LXT) == "REM Sleep"
        XTrainLabelK(3,LXT) = 1;
    end
    if YTrains(LXT) == "Not Sleep"
        YTrainK(1,LXT) = 1;
    elseif YTrains(LXT) == "NREM Sleep"
        YTrainK(2,LXT) = 1;
    elseif YTrains(LXT) == "REM Sleep"
        YTrainK(3,LXT) = 1;
    end
end

for LXT = 1:length(XTestLabels)
    if XTestLabels(LXT) == "Not Sleep"
        XTestLabelK(1,LXT) = 1;
    elseif XTestLabels(LXT) == "NREM Sleep"
        XTestLabelK(2,LXT) = 1;
    elseif XTestLabels(LXT) == "REM Sleep"
        XTestLabelK(3,LXT) = 1;
    end
    if YTests(LXT) == "Not Sleep"
        YTestK(1,LXT) = 1;
    elseif YTests(LXT) == "NREM Sleep"
        YTestK(2,LXT) = 1;
    elseif YTests(LXT) == "REM Sleep"
        YTestK(3,LXT) = 1;
    end
end

RF_confMat = figure; 
subplot(1,2,1)
plotconfusion(YTrainK,XTrainLabelK);
title('Train Confusion')
savefig(RF_confMat,[dirpath animalID '_FP_RF_Train_ConfusionMatrix']);
subplot(1,2,2)
plotconfusion(YTestK,XTestLabelK);
title('Test Confusion')
savefig(RF_confMat,[dirpath animalID '_FP_RF_Test_ConfusionMatrix']);
close(RF_confMat)
% save labels for later confusion matrix
ConfusionData.RF.trainYlabels = YTrain.behavState;
ConfusionData.RF.trainXlabels = XTrainLabels;
ConfusionData.RF.testYlabels = YTest.behavState;
ConfusionData.RF.testXlabels = XTestLabels;

cd(startingDirectory)

% save confusion matrix results
save([dirpath animalID '_ConfusionData.mat'],'ConfusionData')
% end
