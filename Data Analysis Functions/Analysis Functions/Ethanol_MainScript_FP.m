function [] = Ethanol_MainScript_FP()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. eLife2020
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________
% addpath(genpath('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Documents\Research_Codes\FiberPhotometry\Data-Analysis-master'))
clear; clc; close all;
%% Animal IDs
FP_animalIDs =    {'NEACh006'};%,'NEACh006'};% 
%% make sure the code repository and data are present in the current directory
firstHrs = "false";
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'n';
saveFigs = 'y';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
%     multiWaitbar('Analyzing sleep probability',0,'Color','B'); pause(0.25);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults] = AnalyzeData(rootFolder,FP_animalIDs);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end

%% main figure panels
% [AnalysisResults] = Fig1_S4_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
% [AnalysisResults] = Fig1_S3_Whisk_GRABNE_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
[AnalysisResults] = Fig1_S2_Stim_GRABNE_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);

% [AnalysisResults] = Fig1_S5_Movement_GRABNE_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
% [AnalysisResults] = Fig4_FP_Transition_GRABNE_SingleMouse_consolidated_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs);
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData(rootFolder,FP_animalIDs)

saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end
% these data are not from first hours
firstHrs = "false";
%% Block [1] Analyze the arousal-state probability of trial duration and resting events (IOS)
runFromStart = 'n';
for aa = 1:length(FP_animalIDs)
    if isfield(AnalysisResults,(FP_animalIDs{1,aa})) == false || isfield(AnalysisResults.(FP_animalIDs{1,aa}),'SleepProbability') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeAwakeProbability_GRBANE(FP_animalIDs{1,aa},saveFigs,rootFolder,AnalysisResults,firstHrs);
    end
    multiWaitbar('Analyzing sleep probability','Value',aa/length(FP_animalIDs));
end
%% Block [2] Analyze the arousal-state distribution of different behavioral measurements (IOS)
runFromStart = 'n';
for bb = 1:length(FP_animalIDs)
    if isfield(AnalysisResults,(FP_animalIDs{1,bb})) == false || isfield(AnalysisResults.(FP_animalIDs{1,bb}),'BehaviorDistributions') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeBehavioralDistributions_GRABNE(FP_animalIDs{1,bb},rootFolder,AnalysisResults,firstHrs);
    end
    multiWaitbar('Analyzing behavioral distributions','Value',bb/length(FP_animalIDs));
end
%% Block [3] Analyze the transitions between different arousal-states (IOS)
runFromStart = 'n';
for dd = 1:length(FP_animalIDs)
    if isfield(AnalysisResults,(FP_animalIDs{1,dd})) == false || isfield(AnalysisResults.(FP_animalIDs{1,dd}),'Transitions') == false || strcmp(runFromStart,'y') == true
          [AnalysisResults] = AnalyzeTransitionalAverages_FP_Movements_GRABNE(FP_animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults,firstHrs);

    end
    multiWaitbar('Analyzing behavioral transitions triggered changes','Value',dd/length(FP_animalIDs));
end
%% Block [4] Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
% runFromStart = 'n';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanRhodamine') == false || strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeMeanRhodamine_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%                 [AnalysisResults] = AnalyzeMeanRhodamine_FP_GRABNE_Sleep(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing behavioral Rhodamine','Value',ff/length(FP_animalIDs));
% end
%% Block [5] Analyze the hemodynamic signal GRAB during different arousal states (IOS)
% runFromStart = 'n';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanGFP') == false || strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%                 [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE_Sleep(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing behavioral GFP','Value',ff/length(FP_animalIDs));
% end
%% Block [5] Analyze the Pupil diameter during different arousal states (IOS)
% runFromStart = 'y';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanGFP') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeMeanPupilDiameter_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing behavioral GFP','Value',ff/length(FP_animalIDs));
% end

%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% runFromStart = 'y';
% for pp = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Movement') == false ||isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
% end
%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% runFromStart = 'y';
% for pp = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeOptoStimEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs);
%         [AnalysisResults] = AnalyzeOptoStimEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
% end
%% Block [7] Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'y';
% for jj = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,jj})) == false || isfield(AnalysisResults.(FP_animalIDs{1,jj}),'Coherence') == false || strcmp(runFromStart,'y') == true
%         % [AnalysisResults] = AnalyzeCoherence_Pupil_GRABNE(FP_animalIDs{1,jj},rootFolder,AnalysisResults,firstHrs);
%         [AnalysisResults] = AnalyzeCoherence_DataDivision(FP_animalIDs{1,jj},rootFolder,AnalysisResults,firstHrs);
% 
%     end
%     multiWaitbar('Analyzing coherence','Value',jj/length(FP_animalIDs));
% end
%% Block [8] Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
% runFromStart = 'n';
% for jj = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,jj})) == false || isfield(AnalysisResults.(FP_animalIDs{1,jj}),'NeuralHemoCoherence') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeNeuralHemoCoherence(FP_animalIDs{1,jj},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing neural-hemo coherence','Value',jj/length(FP_animalIDs));
% end
%% Block [9] Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'y';
% for kk = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,kk})) == false || isfield(AnalysisResults.(FP_animalIDs{1,kk}),'PowerSpectrum') == false || strcmp(runFromStart,'y') == true
% %         [AnalysisResults] = AnalyzePowerSpectrum(FP_animalIDs{1,kk},rootFolder,AnalysisResults);
%           % [AnalysisResults] = AnalyzePowerSpectrum_Pupil_GRABNE(FP_animalIDs{1,kk},rootFolder,AnalysisResults);
%           [AnalysisResults] = AnalyzePowerSpectrum_DataDivision(FP_animalIDs{1,kk},rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing power spectra','Value',kk/length(FP_animalIDs));
% end
%% Block [10] Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'y';
% for mm = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,mm})) == false || isfield(AnalysisResults.(FP_animalIDs{1,mm}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeCorrCoeffs(FP_animalIDs{1,mm},rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing Pearson''s correlation coefficients','Value',mm/length(FP_animalIDs));
% end

 %% Block [11] Analyze the cross-correlation between neural activity and hemodynamics [HbT] 
% runFromStart = 'y';
% for nn = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,nn})) == false || isfield(AnalysisResults.(FP_animalIDs{1,nn}),'CrossCorrelation') == false || strcmp(runFromStart,'y') == true
% %         [AnalysisResults] = AnalyzeCrossCorrelation_GRABNE(FP_animalIDs{1,nn},rootFolder,AnalysisResults,firstHrs);
%         % [AnalysisResults] = AnalyzeCrossCorrelation_GRABNE_CBV(FP_animalIDs{1,nn},rootFolder,AnalysisResults,firstHrs);
%         [AnalysisResults] = AnalyzeCrossCorrelation_DataDivision(FP_animalIDs{1,nn},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing cross correlation','Value',nn/length(FP_animalIDs));
% end
%% Block [14] Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults.(FP_animalIDs{1,qq}),'TRITCvsGamma') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeCBVGammaRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing CBV-Gamma relationship','Value',qq/length(FP_animalIDs));
% end
 %% Block [15] Analyze the relationship between gamma-band power and GCaMP7s
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults.(FP_animalIDs{1,qq}),'GCaMP7svsGamma') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeGCaMP7sGammaRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing GCaMP7s-Gamma relationship','Value',qq/length(FP_animalIDs));
% end
%% Block [16] Analyze the relationship between HBT and GCaMP7s
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults.(FP_animalIDs{1,qq}),'TRITCvsGCaMP7s') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeHbTGCaMP7sRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults);
%     end
%     multiWaitbar('Analyzing CBV-GCaMP7s relationship','Value',qq/length(FP_animalIDs));
% end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
