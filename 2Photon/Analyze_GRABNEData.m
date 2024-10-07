
function [] = Analyze_GRABNEData()
    currentFolder = pwd;
    addpath(genpath(currentFolder));
    fileparts = strsplit(currentFolder,filesep);
    
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
  
    % add root folder to Matlab's working directory
    addpath(genpath(rootFolder))
    % FP animal IDs
    FP_animalIDs = {'NEACh008'};
    saveFigs = 'y';
    if exist('AnalysisResults.mat','file') == 2
        load('AnalysisResults.mat')
    else
        AnalysisResults = [];
    end
    %% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
    runFromStart = 'y';
    for pp = 1:length(FP_animalIDs)
        if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,FP_animalIDs);
        end
        multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
    end
    %% Plot the response
    
    [AnalysisResults] = Whisk_GRABNE_Response_2Photon(rootFolder,saveFigs,delim,AnalysisResults,FP_animalIDs);
    [AnalysisResults] = Stim_GRABNE_Response_2Photon(rootFolder,saveFigs,delim,AnalysisResults,FPanimalIDs);

end
