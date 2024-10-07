function [ScoringResults] = PredictBehaviorEvents_FP_S(animalID,startingDirectory,baselineDirectory,modelDataFileIDs,modelName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

disp(['Predicting behavior events using ' modelName ' model']); disp(' ')
% load appropriate model
modelDirectory = [startingDirectory '\Figures\Sleep Models\'];
cd(modelDirectory)
notManual = 'y';
if strcmp(modelName,'CNN') == true
    modelName = [animalID '_FP_CNN_SleepScoringModel.mat'];
    load(modelName)
    modelType = 'CNN';
    MDL = CNN_MDL;
elseif strcmp(modelName,'Manual') == true
    notManual = 'n';
    modelType = 'Manual';
end
cd(baselineDirectory)
% go through each manual file and sleep score it using the chosen model
if strcmp(notManual,'y') == true
    for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    dataLength = size(ModelData.TrainLabels,1);
    for lp = 1:1:dataLength
        AllFeatures{((bb-1)*dataLength)+lp } = ModelData.Features{lp}';
    end
    
    AllLabels(((bb-1)*dataLength)+1 : (bb*dataLength)) = ModelData.TrainLabels;
    end
    scoringTable = AllFeatures;
    [labels,~] = predict(MDL,scoringTable);
    % apply a logical patch on the REM events
    REMindex = strcmp(labels,'REM Sleep');
    numFiles = length(labels)/dataLength;
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for c = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,c);
        patchedREMarray = LinkBinaryEvents_FP(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for d = 1:length(labels)
        if patchedREMindex(d,1) == 1
            labels{d,1} = 'REM Sleep';
        end
    end
    % export results
    reshapedLabels = reshape(labels,dataLength,numFiles);
    for e = 1:size(modelDataFileIDs,1)
        modelDataFileID = modelDataFileIDs(e,:);
        [~,~,fileID] = GetFileInfo_FP(modelDataFileID);
        fileIDs{e,1} = fileID;
        labelArrays{e,1} = reshapedLabels(:,e);
    end
    ScoringResults.fileIDs = fileIDs;
    ScoringResults.labels = labelArrays;
    ScoringResults.allfileIDs = joinedFileList;
    ScoringResults.alllabels = labels;
else
    ScoringResults = [];
end
save([animalID '_' modelType '_ScoringResults'],'ScoringResults')
cd(startingDirectory)

end
