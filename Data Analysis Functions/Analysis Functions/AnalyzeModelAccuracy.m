function [AnalysisResults] = AnalyzeModelAccuracy(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the out-of-bag error (model accuracy) of each random forest classification model
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = { 'T281','T282','T283','T284','T285'};
%% only run analysis for valid animal IDs
modelLocation = [rootFolder '\' animalID '\Figures\Sleep Models\'];
cd(modelLocation)
% load the random forest model for evaluating the cross-validation
modelName = [animalID '_FP_RF_SleepScoringModel.mat'];
load(modelName,'-mat')
iterations = 100;
X = RF_MDL.X;
Y = RF_MDL.Y;
% determine the misclassification probability (for classification trees) for out-of-bag observations in the training data
AnalysisResults.(animalID).ModelAccuracy.oobErr = oobError(RF_MDL,'Mode','Ensemble')*100;
% re-create the model 100 times with shuffled data and determine the oobError distribution
for aa = 1:iterations
    shuffYIdx = randperm(numel(Y));
    shuffY = Y(shuffYIdx);
    numTrees = 128;
    shuffRF_MDL = TreeBagger(numTrees,X,shuffY,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    AnalysisResults.(animalID).ModelAccuracy.shuff_oobErr(aa,1) = oobError(shuffRF_MDL,'Mode','Ensemble')*100;
end
% save figure if desired
if strcmp(saveFigs,'y') == true
    distributionFig = figure;
    histogram(AnalysisResults.(animalID).ModelAccuracy.shuff_oobErr,'Normalization','probability','FaceColor','k')
    title('OOB error for shuffled data')
    ylabel('Probability')
    xlabel('error (%)')
    axis square
    set(gca,'box','off')
    % Save the figure to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Cross Validation/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(distributionFig,[dirpath animalID '_OOBErrorDistribution']);
    close(distributionFig)
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
