function [] = UpdateTrainingDataSets_FP(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Add existing labels to (potentially) re-processed model tables
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
    trainingDataSetID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    load(modelDataSetID)
    load(trainingDataSetID)
    disp(['Updating training data set for ' trainingDataSetID '...' ]); disp(' ')
    behaveState =   trainingTable.behavState; %traininglables;%
   
    % I have removed extra 3 minutes of data from each side after I did the
    % sleep scoring last time.
%     if length(behaveState) ~= size(trainingTable,1)
%         removedIDx = (3*60)/5;
%         behaveState = behaveState(removedIDx+1:end-removedIDx);
%     end
    paramsTable.behavState = behaveState;
   
    trainingTable = paramsTable;
    save(trainingDataSetID,'trainingTable')
end
