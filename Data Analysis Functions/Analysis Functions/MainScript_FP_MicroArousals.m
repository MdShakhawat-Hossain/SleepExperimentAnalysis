function [] = MainScript_FP_MicroArousals()
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
%% make sure the code repository and data are present in the current directory
% currentFolder = 'H:\Sleep_GCaMP7s_ChATCre\';
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
    [AnalysisResults] = AnalyzeData(rootFolder);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end
%% supplemental figure panels
% [AnalysisResults] = Fig8_S1(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig7_S3(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig7_S2(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig7_S1(rootFolder,saveFigs,delim,AnalysisResults);
%% [AnalysisResults] = Fig6_S1(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig5_S1(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig3_S5(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig3_S4(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig3_S3(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig3_S2(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig3_S1(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig2_S2(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig2_S1(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S9(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S8(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S7(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S6(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S5(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig4_S1_GRABNE_Rhodamine(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig4_S1_GCaMP7s(rootFolder,saveFigs,delim,AnalysisResults);

% [AnalysisResults] = Fig4_FP_Bilateral_No_REM_TRITCChAT_movement(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig4_FP_Bilateral_No_REM_TRITCChAT_movement_compare(rootFolder,saveFigs,delim,AnalysisResults);


% [AnalysisResults] = Fig1_S1(rootFolder,saveFigs,delim,AnalysisResults);
%% main figure panels
% [AnalysisResults] = Fig1_S4_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig1_S3_Whisk_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig1_S2_Stim_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);

% [AnalysisResults] = PlotCOherence_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = PlotpowerSpectrum_GRABNE(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = PlotCrossCorrelation_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig4_FP_Transition_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig4_FP_Transition_GRABNE_SingleMouse(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);
% [AnalysisResults] = Fig4_FP_Transition_GRABNE_SingleMouse_consolidated(rootFolder,saveFigs,delim,AnalysisResults,firstHrs);

% [AnalysisResults] = Fig1_S5_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S6_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults);
% [AnalysisResults] = Fig1_S7_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults);
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData(rootFolder)
% FP animal IDs
FP_animalIDs = {'NEACh001'};

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
runFromStart = 'y';
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
%         [AnalysisResults] = AnalyzeMeanRhodamine_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral Rhodamine','Value',ff/length(FP_animalIDs));
% end
%% Block [5] Analyze the hemodynamic signal [GCaMP7s] during different arousal states (IOS)
% runFromStart = 'n';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults.(FP_animalIDs{1,ff}),'MeanGFP') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral GFP','Value',ff/length(FP_animalIDs));
% end
%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% runFromStart = 'y';
% for pp = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults,firstHrs);
%     end
%     multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
% end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
